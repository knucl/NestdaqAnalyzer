#include <fairmq/Device.h>
#include <fairmq/runDevice.h>

#include <cctype>
#include <filesystem>
#include <iomanip>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "TSystem.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THttpServer.h"

#include "NodeUnpacker.hh"
#include "LeafProcessor.hh"
#include "LeafSchema.hh"
#include "NodeHeaderSchema.hh"
#include "SubTimeFrameHeader.h"

// Helper timer
#include "KTimer.cxx"

#include "detector_map.hpp"

namespace bpo = boost::program_options;

struct OnlineAnalysisNode : fair::mq::Device
{
    struct OptionKey {
        static constexpr std::string_view InputChannelName {"in-chan-name"};
        static constexpr std::string_view SamplingMode     {"sampling-mode"};
        static constexpr std::string_view Prescale         {"prescale"};
        static constexpr std::string_view ChannelMapName   {"channel-map"};
        static constexpr std::string_view DefinitionName   {"detector-definition"};
    };

    OnlineAnalysisNode()
        : fLeafProcessor(new nestdaq::unpacker::LeafProcessor())
    {
    }

    void InitTask() override
    {
        using opt = OptionKey;
        fInputChannelName = fConfig->GetValue<std::string>(opt::InputChannelName.data());
        fSamplingMode     = fConfig->GetValue<std::string>(opt::SamplingMode.data());
        fPrescale         = fConfig->GetValue<int>(opt::Prescale.data());
        fChannelMapName   = fConfig->GetValue<std::string>(opt::ChannelMapName.data());
        fDefinitionName   = fConfig->GetValue<std::string>(opt::DefinitionName.data());
//        LOG(info) << "DEBUG: InitTask Started. Channel: " << fInputChannelName << " Mode: " << fSamplingMode << " Prescale: " << fPrescale;
        LOG(info) << "DEBUG: InitTask Started. Channel: " << fInputChannelName << " Mode: " << fSamplingMode << " Prescale: " << fPrescale << " detector definition: " << fDefinitionName << " channel map: " << fChannelMapName;

        gROOT->SetBatch(kTRUE);

        if (!gApplication) {
            static int argc = 1;
            static char* argv[] = {(char*)"OnlineAnalysisNode", nullptr};
            new TApplication("OnlineAnalysisNode", &argc, argv);
        }

        if (!fServer) {
            fServer = new THttpServer("http:8888");
            fServer->SetReadOnly(kTRUE);
            LOG(info) << "DEBUG: THttpServer started on port 8888";
        }

        fDrawTimer.SetDuration(100); 

        // Global Histograms (All FEMs aggregated) - Guard against duplicates
        if (!fH2HitPattern) {
             fH2HitPattern = new TH2F("h2_hitpat", "Hit Pattern;Channel;FEM ID", 128, 0, 128, 10, 0, 10);
             if(fServer) {
                 fServer->Register("/Summary", fH2HitPattern);
                 
                 // Update monitor settings
                 fServer->SetItemField("/", "_monitoring", "1000");
                 fServer->SetItemField("/Summary", "_monitoring", "1000");
             }
        } else {
            fH2HitPattern->Reset();
        }
        
        // Hits on FE channels that the channel map does not know (FE channel x FE board)
        if (!fH2Unmapped) {
            fH2Unmapped = new TH2F("h2_unmapped", "Hits on FE channels not in the channel map;FE channel;FE board",
                                   128, 0, 128, 1, 0, 1);
            fH2Unmapped->SetCanExtend(TH1::kYaxis);
            if (fServer) fServer->Register("/Detector", fH2Unmapped);
        }
        ReloadDetectorMap();

        LOG(info) << "DEBUG: InitTask Finished.";
    }

    void PreRun() override
    {
        LOG(info) << "OnlineAnalysisNode: PreRun - Resetting histograms for new run.";
        ReloadDetectorMap();
        ResetHistograms();
    }

    // ---- channel map -------------------------------------------------------
    // One set of histograms per detector, plane and data kind. The x axis is
    // segment x ch (bin = segment index * number of ch + ch index).
    struct DetHists {
        TH1F* hitpat = nullptr;
        TH2F* tdc    = nullptr;  // made at the first hit, with the TDC range of its FE type
        TH2F* tot    = nullptr;
        std::string dir;         // THttpServer folder
        std::string name;        // ROOT name stem
        std::string title;
        std::string xtitle;
        int  nbins   = 0;
        bool use_seg = false;
        bool use_ch  = false;
    };
    using MapStamp = std::pair<std::filesystem::file_time_type, std::filesystem::file_time_type>;

    static std::string Sanitize(const std::string& s)
    {
        std::string r;
        for (char c : s) r += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
        return r;
    }

    MapStamp CurrentMapStamp() const
    {
        std::error_code ec1, ec2;
        return {std::filesystem::last_write_time(fDefinitionName, ec1),
                std::filesystem::last_write_time(fChannelMapName, ec2)};
    }

    // Reads the detector definition and the channel map when they changed.
    // A map with errors is not used: the map of the previous run stays, and
    // every problem is logged with its file and line.
    void ReloadDetectorMap()
    {
        if (fDefinitionName.empty() || fChannelMapName.empty()) {
            LOG(warn) << "Channel map is not set (detector-definition='" << fDefinitionName
                      << "', channel-map='" << fChannelMapName << "'): no detector histograms";
            return;
        }
        const auto stamp = CurrentMapStamp();
        if (fDetMapLoaded && stamp == fDetMapStamp) {
            LOG(info) << "Channel map unchanged since the last run";
            return;
        }
        chmap::DetectorMap next;
        const bool ok = next.load(fDefinitionName, fChannelMapName);
        for (const auto& d : next.diagnostics()) {
            const std::string where = d.file + (d.line > 0 ? ":" + std::to_string(d.line) : "");
            if (d.level == chmap::Diagnostic::Level::kError) {
                LOG(error) << where << ": " << d.message;
            } else {
                LOG(warn) << where << ": " << d.message;
            }
        }
        if (!ok) {
            LOG(error) << "Channel map has errors: "
                       << (fDetMapLoaded ? "keeping the map of the previous run"
                                         : "no detector histograms until it is fixed");
            return;
        }
        fDetMap = std::move(next);
        fDetMapStamp = stamp;
        fDetMapLoaded = true;
        RebuildDetectorHistograms();
        LOG(info) << "Channel map loaded: " << fDetMap.entries().size() << " FE channels, "
                  << fDetMap.devices().size() << " detectors";
    }

    void DeleteDetectorHistograms()
    {
        for (auto& h : fDetHists) {
            for (TH1* p : {static_cast<TH1*>(h.hitpat), static_cast<TH1*>(h.tdc), static_cast<TH1*>(h.tot)}) {
                if (!p) continue;
                if (fServer) fServer->Unregister(p);
                delete p;
            }
        }
        fDetHists.clear();
        fDetOffset.clear();
    }

    // Bin labels such as "0U", "0D", "1U", or "3-12" for numeric ch, when they are few enough to read.
    static void LabelAxis(TAxis* axis, const chmap::DeviceDef& d, const DetHists& h)
    {
        if (h.nbins > 64) return;
        const int nch = d.chs.size();
        for (int b = 0; b < h.nbins; ++b) {
            const std::string seg = h.use_seg ? d.segments[b / nch] : "";
            const std::string ch  = h.use_ch  ? d.chs[b % nch] : "";
            const bool numeric_ch = !ch.empty() && std::isdigit(static_cast<unsigned char>(ch[0]));
            const std::string label = (h.use_seg && numeric_ch) ? seg + "-" + ch : seg + ch;
            if (!label.empty()) axis->SetBinLabel(b + 1, label.c_str());
        }
    }

    void RebuildDetectorHistograms()
    {
        DeleteDetectorHistograms();
        for (const auto& d : fDetMap.devices()) {
            fDetOffset.push_back(fDetHists.size());
            const bool use_plane = !(d.planes.size() == 1 && d.planes[0] == "-");
            for (std::size_t p = 0; p < d.planes.size(); ++p) {
                for (std::size_t x = 0; x < d.datas.size(); ++x) {
                    DetHists h;
                    h.use_seg = !(d.segments.size() == 1 && d.segments[0] == "-");
                    h.use_ch  = !(d.chs.size() == 1 && d.chs[0] == "-");
                    h.nbins   = d.segments.size() * d.chs.size();
                    h.dir     = "/Detector/" + d.name + (use_plane ? "/" + d.planes[p] : "");
                    h.name    = Sanitize(d.name) + (use_plane ? "_" + Sanitize(d.planes[p]) : "")
                              + (d.datas.size() > 1 ? "_" + Sanitize(d.datas[x]) : "");
                    h.title   = d.name + (use_plane ? " plane " + d.planes[p] : "")
                              + (d.datas.size() > 1 ? " " + d.datas[x] : "");
                    h.xtitle  = h.use_seg && h.use_ch ? "segment, ch" : h.use_seg ? "segment" : h.use_ch ? "ch" : "";
                    h.hitpat  = new TH1F(("h_" + h.name + "_hitpat").c_str(),
                                         (h.title + " hit pattern;" + h.xtitle + ";counts").c_str(),
                                         h.nbins, 0, h.nbins);
                    LabelAxis(h.hitpat->GetXaxis(), d, h);
                    if (fServer) fServer->Register(h.dir.c_str(), h.hitpat);
                    fDetHists.push_back(std::move(h));
                }
            }
        }
        if (fServer) fServer->SetItemField("/Detector", "_monitoring", "1000");
    }

    void CreateDetectorTimeHistograms(DetHists& h, const chmap::DeviceDef& d, bool is_hr)
    {
        const double tdc_max = is_hr ? 600000000.0 : 250000.0;
        const double tot_max = is_hr ? 100000.0 : 200.0;
        h.tdc = new TH2F(("h_" + h.name + "_tdc").c_str(), (h.title + " TDC;" + h.xtitle + ";TDC").c_str(),
                         h.nbins, 0, h.nbins, 200, 0, tdc_max);
        h.tot = new TH2F(("h_" + h.name + "_tot").c_str(), (h.title + " TOT;" + h.xtitle + ";TOT").c_str(),
                         h.nbins, 0, h.nbins, 200, 0, tot_max);
        LabelAxis(h.tdc->GetXaxis(), d, h);
        LabelAxis(h.tot->GetXaxis(), d, h);
        if (fServer) {
            fServer->Register(h.dir.c_str(), h.tdc);
            fServer->Register(h.dir.c_str(), h.tot);
        }
    }

    void FillDetector(const chmap::Address& a, uint32_t tdc, uint32_t tot, bool is_hr)
    {
        const auto& d = fDetMap.devices()[a.device];
        auto& h = fDetHists[fDetOffset[a.device] + a.plane * d.datas.size() + a.data];
        const int bin = a.segment * d.chs.size() + a.ch;
        h.hitpat->Fill(bin);
        if (!h.tdc) CreateDetectorTimeHistograms(h, d, is_hr);
        h.tdc->Fill(bin, tdc);
        h.tot->Fill(bin, tot);
    }

    void FillUnmapped(uint32_t fem_id, uint32_t ch)
    {
        const std::string ip = ToIPAddress(fem_id);
        if (fH2Unmapped) fH2Unmapped->Fill(ch, ip.c_str(), 1.);
        if (fUnmappedSeen.insert((static_cast<uint64_t>(fem_id) << 32) | ch).second) {
            LOG(warn) << "FE " << ip << " ch " << ch << " is not in the channel map (reported once per run)";
        }
    }

    void ResetHistograms()
    {
        if (fH2HitPattern) fH2HitPattern->Reset();
        for (auto& [id, h] : fMapHitPattern) {
            if (h) h->Reset();
        }
        for (auto& [id, vec] : fMapTDC) {
            for (auto* h : vec) if (h) h->Reset();
        }
        for (auto& [id, vec] : fMapTOT) {
            for (auto* h : vec) if (h) h->Reset();
        }
        
        for (auto& h : fDetHists) {
            if (h.hitpat) h.hitpat->Reset();
            if (h.tdc) h.tdc->Reset();
            if (h.tot) h.tot->Reset();
        }
        if (fH2Unmapped) fH2Unmapped->Reset();
        fUnmappedSeen.clear();

        fProcessCount = 0;
        LOG(info) << "OnlineAnalysisNode: Histograms reset completed.";
    }

    bool ConditionalRun() override
    {
        try {
            LOG(info) << "ConditionalRun: Start. Checking channel: " << fInputChannelName;
            
            // Check if the channel is registered
            bool channel_exists = false;
            for (auto const& [name, channel] : fChannels) {
                if (name == fInputChannelName) {
                    channel_exists = true;
                    break;
                }
            }

            if (!channel_exists) {
                std::stringstream ss;
                ss << "Channel " << fInputChannelName << " is not registered! (Available: ";
                for (auto const& [name, channel] : fChannels) ss << name << " ";
                ss << ")";
                LOG(error) << ss.str();
                std::this_thread::sleep_for(std::chrono::seconds(1));
                return true; 
            }

            FairMQParts parts;
            // timeout 100ms
            if (Receive(parts, fInputChannelName, 0, 100) > 0) {
                LOG(info) << "ConditionalRun: Received parts: " << parts.Size(); 
                
                // Conditionally drain queue based on sampling mode
                if (fSamplingMode == "latest") {
                    // Drain queue: only keep the most recent message
                    while(true) {
                        FairMQParts newer;
                        if (Receive(newer, fInputChannelName, 0, 0) > 0) {
                            parts = std::move(newer);
                        } else {
                            break;
                        }
                    }
                }
                
                fProcessCount++;
                if (fProcessCount % fPrescale != 0) {
                    if (fDrawTimer.Check()) {
                        UpdateDisplay();
                    }
                    return true; // Skip processing to save CPU
                }

                // Prepare Processor
                fLeafProcessor->clear();
                LOG(info) << "LeafProcessor cleared. Processing parts...";

                std::unordered_map<uint32_t, uint32_t> current_fem_types;

                // Merge all parts into a single aligned buffer
                size_t total_size = 0;
                for (const auto& msg : parts) {
                    total_size += msg->GetSize();
                }

                if (total_size == 0) return true;

                size_t n_words = (total_size + 7) / 8;
                std::vector<uint64_t> merged_buffer(n_words, 0);
                uint8_t* dest_ptr = reinterpret_cast<uint8_t*>(merged_buffer.data());
                size_t offset = 0;
                for (const auto& msg : parts) {
                    std::memcpy(dest_ptr + offset, msg->GetData(), msg->GetSize());
                    offset += msg->GetSize();
                }

                uint64_t* data = merged_buffer.data();
                uint64_t* end_ptr = data + n_words;
                uint64_t* current_ptr = data;

                while (current_ptr < end_ptr) {
                    uint64_t word = *current_ptr;
                    // Mask top 8 bits (Version) to match registry: 0x00...
                    uint64_t magic = word & 0x00FFFFFFFFFFFFFF; 
                    
                    // 1. Handle Known Containers that we want to scan inside manually
                    // TimeFrame-v1 (TIMEFRAM) and @TF-HEAD (v0)
                    if (magic == 0x004d5246454d4954 || magic == 0x004145482d465440) {
                         LOG(info) << "Found Container Header (skipping 24-byte header): 0x" << std::hex << magic << std::dec;
                         current_ptr += 3; 
                         continue;
                    }

                    if (magic == SubTimeFrame::MAGIC) {
                        auto* pstf = reinterpret_cast<SubTimeFrame::Header*>(current_ptr);
                        current_fem_types[pstf->femId] = pstf->femType;
                    }

                    // 2. Handle Recognized Headers (like SubTimeFrame)
                    if (nestdaq::g_header_schema.count(magic) > 0) {
                        LOG(info) << "Found Recognized Header: 0x" << std::hex << magic << std::dec;
                        auto unpacker = std::make_shared<nestdaq::unpacker::NodeUnpacker>(nestdaq::g_header_schema.at(magic));
                        unpacker->set_data(current_ptr);
                        unpacker->unpack();
                        
                        auto header_data = unpacker->get_header();

                        // Extract leaf nodes (HBFs/TDC blocks)
                        auto leaves = unpacker->extract_leafnodes();
                        for (auto& leaf : leaves) {
                            fLeafProcessor->set_leaf_node(std::get<0>(leaf), std::get<1>(leaf), std::get<2>(leaf));
                        }

                        // Determine the total length to skip for this block
                        uint64_t length_bytes = header_data.count("Length") ? header_data["Length"] : 8;
                        if (length_bytes == 0) length_bytes = 8;
                        
                        current_ptr += (length_bytes / 8);
                        continue;
                    }
                    
                    // 3. Skip unrecognized words
                    current_ptr++;
                }

                // Decode and Process Data
                if (fLeafProcessor->get_num_frame() > 0) {
                    fLeafProcessor->decode_node_body();
                    fLeafProcessor->decode_heartbeat_delimiter(); // Required to avoid crash in get_leafnode_data
                }
                
                // Iterate processed FEMs
                if (fLeafProcessor->get_num_frame() > 0) {
                    auto fem_ids = fLeafProcessor->get_node_ids();
                    
                    for (auto fem_id : fem_ids) {
                        // get_leafnode_data returns {body_vec, hbd_vec}
                        auto [body_vec, hbd_vec] = fLeafProcessor->get_leafnode_data(fem_id);
                    
                        LOG(info) << "FEM " << fem_id << ": " << body_vec.size() << " body datasets, " << hbd_vec.size() << " hbd datasets";

                        int hits_count = 0;
                        for (const auto& dataset : body_vec) {
                            const std::string& type_name = std::get<0>(dataset);
                            const auto& leaf_data = std::get<1>(dataset);
                            
                            if (type_name == "L-TDC" || type_name == "T-TDC") {
                                bool is_hr = false;
                                if (current_fem_types.count(fem_id) > 0) {
                                    uint32_t fem_type = current_fem_types.at(fem_id);
                                    is_hr = (fem_type == SubTimeFrame::TDC64H || fem_type == SubTimeFrame::TDC64H_V3);
                                } else {
                                    is_hr = (type_name == "T-TDC" || type_name == "H-TDC");
                                }

                                CheckAndCreateHistograms(static_cast<uint32_t>(fem_id), is_hr);

                                if (leaf_data.count("Ch") && leaf_data.count("TDC") && leaf_data.count("TOT")) {
                                    uint32_t ch = leaf_data.at("Ch");
                                    uint32_t tdc = leaf_data.at("TDC");
                                    uint32_t tot = leaf_data.at("TOT");
                                    
                                    if (fH2HitPattern) fH2HitPattern->Fill(ch, fem_id);
                                    if (fMapHitPattern.count(fem_id)) fMapHitPattern[fem_id]->Fill(ch);

                                    int max_ch = is_hr ? 64 : 128;
                                    if (ch < static_cast<uint32_t>(max_ch) && fMapTDC.count(fem_id) && fMapTDC[fem_id].size() > ch) {
                                        fMapTDC[fem_id][ch]->Fill(tdc);
                                        fMapTOT[fem_id][ch]->Fill(tot);
                                    }
                                    if (fDetMapLoaded) {
                                        if (const auto* a = fDetMap.find(static_cast<uint32_t>(fem_id), ch)) {
                                            FillDetector(*a, tdc, tot, is_hr);
                                        } else {
                                            FillUnmapped(static_cast<uint32_t>(fem_id), ch);
                                        }
                                    }

                                    hits_count++;
                                }
                            }
                        }
                        if (hits_count > 0) {
                            LOG(info) << "FEM " << fem_id << ": Processed " << hits_count << " hits";
                        }
                    }
                }
                }
                
            if (fDrawTimer.Check()) {
                UpdateDisplay();
            }

        } catch (const std::exception& e) {
            LOG(warn) << "Exception in ConditionalRun: " << e.what();
            LOG(warn) << "Assuming End-of-Stream. Entering Keep-Alive mode to maintain HTTP Server.";
            
            // Keep Alive Loop
            while (true) {
                if (fDrawTimer.Check()) {
                    UpdateDisplay();
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                
                // Allow State Machine transitions/interrupts
                if (GetCurrentState() != fair::mq::State::Running) {
                    break;
                }
            }
        } catch (...) {
            LOG(error) << "Unknown exception in ConditionalRun. Retrying...";
        }

        return true;
    }

    std::string ToIPAddress(uint32_t fem_id) {
        uint8_t* ip = reinterpret_cast<uint8_t*>(&fem_id);
        return std::to_string(ip[3]) + "." + std::to_string(ip[2]) + "." + std::to_string(ip[1]) + "." + std::to_string(ip[0]);
    }

    void CheckAndCreateHistograms(uint32_t fem_id, bool is_hr) {
        if (fMapTDC.count(fem_id) > 0) return; // Already exists

        std::string ip_str = ToIPAddress(fem_id);
        int max_ch = is_hr ? 64 : 128;
        LOG(info) << "Registering New FEM ID: " << fem_id << " (" << ip_str << ") - " << (is_hr ? "HR" : "LR");

        fMapTDC[fem_id].resize(max_ch, nullptr);
        // Create Hit Pattern histogram for this FEM
        fMapHitPattern[fem_id] = new TH1F(Form("h1_hitpat_%s", ip_str.c_str()), 
                                          Form("FEM %s Hit Pattern;Channel;Counts", ip_str.c_str()), 
                                          max_ch, 0, max_ch);

        fMapTOT[fem_id].resize(max_ch, nullptr);
        // Create Histograms
        for (int i = 0; i < max_ch; ++i) {
            double tdc_max = is_hr ? 600000000.0 : 250000.0;
            double tot_max = is_hr ? 100000.0 : 200.0;

            fMapTDC[fem_id][i] = new TH1F(Form("h1_tdc_%s_ch%d", ip_str.c_str(), i), 
                                          Form("FEM %s TDC Ch%d;TDC;Counts", ip_str.c_str(), i), 
                                          1000, 0, tdc_max);
            fMapTOT[fem_id][i] = new TH1F(Form("h1_tot_%s_ch%d", ip_str.c_str(), i), 
                                          Form("FEM %s TOT Ch%d;TOT;Counts", ip_str.c_str(), i), 
                                          1000, 0, tot_max);
            
            if (fServer) {
                fServer->Register(Form("/FEM_%s/TDC", ip_str.c_str()), fMapTDC[fem_id][i]);
                fServer->Register(Form("/FEM_%s/TOT", ip_str.c_str()), fMapTOT[fem_id][i]);
            }
        }
        
        if (fServer) {
            fServer->Register(Form("/FEM_%s", ip_str.c_str()), fMapHitPattern[fem_id]);
        }

        // Create Canvases for this FEM
        TCanvas* c_tdc_1 = new TCanvas(Form("c_%s_tdc_0_31", ip_str.c_str()), Form("FEM %s TDC 0-31", ip_str.c_str()), 1600, 1000);
        c_tdc_1->Divide(8, 4);
        TCanvas* c_tdc_2 = new TCanvas(Form("c_%s_tdc_32_63", ip_str.c_str()), Form("FEM %s TDC 32-63", ip_str.c_str()), 1600, 1000);
        c_tdc_2->Divide(8, 4);
        TCanvas* c_tdc_3 = nullptr;
        TCanvas* c_tdc_4 = nullptr;
        
        if (!is_hr) {
            c_tdc_3 = new TCanvas(Form("c_%s_tdc_64_95", ip_str.c_str()), Form("FEM %s TDC 64-95", ip_str.c_str()), 1600, 1000);
            c_tdc_3->Divide(8, 4);
            c_tdc_4 = new TCanvas(Form("c_%s_tdc_96_127", ip_str.c_str()), Form("FEM %s TDC 96-127", ip_str.c_str()), 1600, 1000);
            c_tdc_4->Divide(8, 4);
        }

        TCanvas* c_tot_1 = new TCanvas(Form("c_%s_tot_0_31", ip_str.c_str()), Form("FEM %s TOT 0-31", ip_str.c_str()), 1600, 1000);
        c_tot_1->Divide(8, 4);
        TCanvas* c_tot_2 = new TCanvas(Form("c_%s_tot_32_63", ip_str.c_str()), Form("FEM %s TOT 32-63", ip_str.c_str()), 1600, 1000);
        c_tot_2->Divide(8, 4);
        TCanvas* c_tot_3 = nullptr;
        TCanvas* c_tot_4 = nullptr;

        if (!is_hr) {
            c_tot_3 = new TCanvas(Form("c_%s_tot_64_95", ip_str.c_str()), Form("FEM %s TOT 64-95", ip_str.c_str()), 1600, 1000);
            c_tot_3->Divide(8, 4);
            c_tot_4 = new TCanvas(Form("c_%s_tot_96_127", ip_str.c_str()), Form("FEM %s TOT 96-127", ip_str.c_str()), 1600, 1000);
            c_tot_4->Divide(8, 4);
        }

        for (int i = 1; i <= 32; ++i) {
            if (auto* p = c_tot_1->GetPad(i)) p->SetLogy(1);
            if (auto* p = c_tot_2->GetPad(i)) p->SetLogy(1);
            if (!is_hr) {
                if (auto* p = c_tot_3->GetPad(i)) p->SetLogy(1);
                if (auto* p = c_tot_4->GetPad(i)) p->SetLogy(1);
            }
        }

        for (int i = 0; i < 32; ++i) {
            c_tdc_1->cd(i+1); fMapTDC[fem_id][i]->Draw();
            c_tot_1->cd(i+1); fMapTOT[fem_id][i]->Draw();
            c_tdc_2->cd(i+1); fMapTDC[fem_id][i+32]->Draw();
            c_tot_2->cd(i+1); fMapTOT[fem_id][i+32]->Draw();
            
            if (!is_hr) {
                c_tdc_3->cd(i+1); fMapTDC[fem_id][i+64]->Draw();
                c_tot_3->cd(i+1); fMapTOT[fem_id][i+64]->Draw();
                c_tdc_4->cd(i+1); fMapTDC[fem_id][i+96]->Draw();
                c_tot_4->cd(i+1); fMapTOT[fem_id][i+96]->Draw();
            }
        }

        fCanvases.push_back(c_tdc_1);
        fCanvases.push_back(c_tdc_2);
        if (!is_hr) {
            fCanvases.push_back(c_tdc_3);
            fCanvases.push_back(c_tdc_4);
        }
        
        fCanvases.push_back(c_tot_1);
        fCanvases.push_back(c_tot_2);
        if (!is_hr) {
            fCanvases.push_back(c_tot_3);
            fCanvases.push_back(c_tot_4);
        }

        if (fServer) {
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_1);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_2);
            if (!is_hr) {
                fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_3);
                fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_4);
            }
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_1);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_2);
            if (!is_hr) {
                fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_3);
                fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_4);
            }

            fServer->SetItemField(Form("/FEM_%s", ip_str.c_str()), "_monitoring", "1000");
        }
    }

    void UpdateDisplay()
    {
        for (auto* c : fCanvases) {
            if (c) {
                for(int i=0; i<32; ++i) { 
                   auto* pad = c->GetPad(i+1);
                   if(pad) pad->Modified(); 
                }
                c->Update();
            }
        }
        if(gSystem) gSystem->ProcessEvents(); 
    }

private:
   std::string fInputChannelName;
   std::string fSamplingMode; // "all" or "latest"
   int fPrescale = 1;
   std::string fChannelMapName;
   int fProcessCount = 0;
   std::unique_ptr<nestdaq::unpacker::LeafProcessor> fLeafProcessor;
   
   KTimer fDrawTimer;

   TH2F* fH2HitPattern = nullptr;
   
   // Map: FEM ID -> Vector of Channel Histograms
   std::map<uint32_t, TH1F*> fMapHitPattern;
   std::map<uint32_t, std::vector<TH1F*>> fMapTDC;
   std::map<uint32_t, std::vector<TH1F*>> fMapTOT;

   // channel map
   std::string fDefinitionName;
   chmap::DetectorMap fDetMap;
   bool fDetMapLoaded = false;
   MapStamp fDetMapStamp;
   std::vector<DetHists> fDetHists;
   std::vector<std::size_t> fDetOffset;   // device index -> first entry in fDetHists
   TH2F* fH2Unmapped = nullptr;
   std::set<uint64_t> fUnmappedSeen;

   std::vector<TCanvas*> fCanvases;
   
   THttpServer* fServer = nullptr;
};

void addCustomOptions(bpo::options_description& options)
{
    using opt = OnlineAnalysisNode::OptionKey;
    options.add_options()
        (opt::InputChannelName.data(), bpo::value<std::string>()->default_value("in"), "Name of the input channel")
        (opt::SamplingMode.data(), bpo::value<std::string>()->default_value("latest"), "Sampling mode: 'all' (process all) or 'latest' (drop old if slow)")
        (opt::Prescale.data(), bpo::value<int>()->default_value(1), "Process only 1 in every N messages (default: 1, i.e., process all)")
        (opt::DefinitionName.data(), bpo::value<std::string>()->default_value(""), "Detector definition CSV (device, plane, segment, ch[, data])")
        (opt::ChannelMapName.data(), bpo::value<std::string>()->default_value(""), "Channel map CSV (fe, fe.ch, device, plane, segment, ch[, data]); re-read at each run start when changed");
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<OnlineAnalysisNode>();
}
