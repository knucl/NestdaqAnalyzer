// DetectorView: detector-ordered histograms for the online analyzers.
//
// It owns the channel map (format v0: detector definition + channel map, see
// channel-map/README.md) and, for each detector, plane and data kind:
//   /Detector/<name>[/<plane>]/h_<name>_hitpat     hit pattern over segment x ch
//   /Detector/<name>[/<plane>]/h_<name>_multi      hits per event (only when events are used)
//   /Detector/<name>[/<plane>]/TDC/, TOT/          one histogram per channel
//   /Detector/<name>[/<plane>]/Canvases/           32 channels per 8x4 canvas
//   /Detector/h2_unmapped                          hits on FE channels not in the map
//
// Usage:
//   InitTask: view.Configure(server, definition_path, map_path, hr_ranges, lr_ranges, use_events);
//             view.Reload();
//   PreRun:   view.Reload(); view.Reset();
//   per hit:  view.Fill(fem_id, ch, tdc, tot, is_hr);
//   events:   view.BeginEvent(); ...Fill()...; view.EndEvent();   (optional)
//   display:  view.UpdateCanvases();
//
// Reload() re-reads the files only when one of them changed. A map with
// errors is not used: the previous map stays and every problem is logged
// with its file and line.

#ifndef NESTDAQ_ANALYZER_DETECTOR_VIEW_H_
#define NESTDAQ_ANALYZER_DETECTOR_VIEW_H_

// Default struct packing for everything below (see detector_map.hpp).
#pragma pack(push)
#pragma pack()

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <fairlogger/Logger.h>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "THttpServer.h"

#include "detector_map.hpp"

class DetectorView {
public:
    // Histogram ranges for one FE type (HR-TDC or LR-TDC).
    struct Ranges {
        int         tdc_bins = 1000;
        double      tdc_min  = 0;
        double      tdc_max  = 1;
        double      tot_max  = 1;
        std::string tdc_axis = "TDC";
    };

    void Configure(THttpServer* server, const std::string& definition_path, const std::string& map_path,
                   const Ranges& hr, const Ranges& lr, bool use_events)
    {
        fServer         = server;
        fDefinitionPath = definition_path;
        fMapPath        = map_path;
        fRangesHR       = hr;
        fRangesLR       = lr;
        fUseEvents      = use_events;
        if (!fUnmapped) {
            fUnmapped = new TH2F("h2_unmapped", "Hits on FE channels not in the channel map;FE channel;FE board",
                                 128, 0, 128, 1, 0, 1);
            fUnmapped->SetCanExtend(TH1::kYaxis);
            if (fServer) fServer->Register("/Detector", fUnmapped);
        }
    }

    bool Loaded() const { return fLoaded; }

    void Reload()
    {
        if (fDefinitionPath.empty() || fMapPath.empty()) {
            LOG(warn) << "Channel map is not set (detector-definition='" << fDefinitionPath
                      << "', channel-map='" << fMapPath << "'): no detector histograms";
            return;
        }
        const auto stamp = CurrentStamp();
        if (fLoaded && stamp == fStamp) {
            LOG(info) << "Channel map unchanged since the last run";
            return;
        }
        chmap::DetectorMap next;
        const bool ok = next.load(fDefinitionPath, fMapPath);
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
                       << (fLoaded ? "keeping the map of the previous run" : "no detector histograms until it is fixed");
            return;
        }
        fMap    = std::move(next);
        fStamp  = stamp;
        fLoaded = true;
        Rebuild();
        LOG(info) << "Channel map loaded: " << fMap.entries().size() << " FE channels, "
                  << fMap.devices().size() << " detectors";
    }

    void Reset()
    {
        for (auto& s : fSets) {
            if (s.hitpat) s.hitpat->Reset();
            if (s.multi)  s.multi->Reset();
            for (auto* p : s.tdc) if (p) p->Reset();
            for (auto* p : s.tot) if (p) p->Reset();
        }
        if (fUnmapped) fUnmapped->Reset();
        fUnmappedSeen.clear();
    }

    void BeginEvent()
    {
        for (auto& s : fSets) s.event_hits = 0;
    }

    void EndEvent()
    {
        for (auto& s : fSets) if (s.multi) s.multi->Fill(s.event_hits);
    }

    // tdc and tot in the units the caller chose for the Ranges of this FE type.
    void Fill(uint32_t fem_id, uint32_t ch, double tdc, double tot, bool is_hr)
    {
        if (!fLoaded) return;
        const auto* a = fMap.find(fem_id, ch);
        if (!a) {
            FillUnmapped(fem_id, ch);
            return;
        }
        const auto& d = fMap.devices()[a->device];
        auto& s = fSets[fOffset[a->device] + a->plane * d.datas.size() + a->data];
        const int bin = a->segment * d.chs.size() + a->ch;
        s.hitpat->Fill(bin);
        ++s.event_hits;
        if (s.tdc.empty()) CreateTimeHistograms(s, d, is_hr ? fRangesHR : fRangesLR);
        s.tdc[bin]->Fill(tdc);
        s.tot[bin]->Fill(tot);
    }

    void UpdateCanvases()
    {
        for (auto& s : fSets) {
            for (auto* c : s.canvases) {
                for (int i = 1; i <= 32; ++i) if (auto* pad = c->GetPad(i)) pad->Modified();
                c->Update();
            }
        }
    }

    const chmap::DetectorMap& Map() const { return fMap; }

private:
    // One set per detector, plane and data kind; x = segment index * number of ch + ch index.
    struct Set {
        TH1F* hitpat = nullptr;
        TH1F* multi  = nullptr;
        std::vector<TH1F*> tdc;   // per channel; made at the first hit, with the ranges of its FE type
        std::vector<TH1F*> tot;
        std::vector<TCanvas*> canvases;
        std::string dir, name, title, xtitle;
        int  nbins   = 0;
        bool use_seg = false;
        bool use_ch  = false;
        int  event_hits = 0;
    };
    using Stamp = std::pair<std::filesystem::file_time_type, std::filesystem::file_time_type>;

    Stamp CurrentStamp() const
    {
        std::error_code ec1, ec2;
        return {std::filesystem::last_write_time(fDefinitionPath, ec1), std::filesystem::last_write_time(fMapPath, ec2)};
    }

    static std::string Sanitize(const std::string& s)
    {
        std::string r;
        for (char c : s) r += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
        return r;
    }

    // "0U", "0D", "1U", or "3-12" when ch is a number.
    static std::string ChannelLabel(const chmap::DeviceDef& d, const Set& s, int bin)
    {
        const int nch = d.chs.size();
        const std::string seg = s.use_seg ? d.segments[bin / nch] : "";
        const std::string ch  = s.use_ch  ? d.chs[bin % nch] : "";
        const bool numeric_ch = !ch.empty() && std::isdigit(static_cast<unsigned char>(ch[0]));
        const std::string label = (s.use_seg && numeric_ch) ? seg + "-" + ch : seg + ch;
        return label.empty() ? std::to_string(bin) : label;
    }

    void Delete()
    {
        for (auto& s : fSets) {
            for (auto* c : s.canvases) {   // canvases first: they draw the histograms below
                if (fServer) fServer->Unregister(c);
                delete c;
            }
            std::vector<TH1*> all{s.hitpat, s.multi};
            all.insert(all.end(), s.tdc.begin(), s.tdc.end());
            all.insert(all.end(), s.tot.begin(), s.tot.end());
            for (TH1* p : all) {
                if (!p) continue;
                if (fServer) fServer->Unregister(p);
                delete p;
            }
        }
        fSets.clear();
        fOffset.clear();
    }

    void Rebuild()
    {
        Delete();
        for (const auto& d : fMap.devices()) {
            fOffset.push_back(fSets.size());
            const bool use_plane = !(d.planes.size() == 1 && d.planes[0] == "-");
            for (std::size_t p = 0; p < d.planes.size(); ++p) {
                for (std::size_t x = 0; x < d.datas.size(); ++x) {
                    Set s;
                    s.use_seg = !(d.segments.size() == 1 && d.segments[0] == "-");
                    s.use_ch  = !(d.chs.size() == 1 && d.chs[0] == "-");
                    s.nbins   = d.segments.size() * d.chs.size();
                    s.dir     = "/Detector/" + d.name + (use_plane ? "/" + d.planes[p] : "");
                    s.name    = Sanitize(d.name) + (use_plane ? "_" + Sanitize(d.planes[p]) : "")
                              + (d.datas.size() > 1 ? "_" + Sanitize(d.datas[x]) : "");
                    s.title   = d.name + (use_plane ? " plane " + d.planes[p] : "")
                              + (d.datas.size() > 1 ? " " + d.datas[x] : "");
                    s.xtitle  = s.use_seg && s.use_ch ? "segment, ch" : s.use_seg ? "segment" : s.use_ch ? "ch" : "";
                    s.hitpat  = new TH1F(("h_" + s.name + "_hitpat").c_str(),
                                         (s.title + " hit pattern;" + s.xtitle + ";counts").c_str(), s.nbins, 0, s.nbins);
                    if (s.nbins <= 64) {
                        for (int b = 0; b < s.nbins; ++b)
                            s.hitpat->GetXaxis()->SetBinLabel(b + 1, ChannelLabel(d, s, b).c_str());
                    }
                    if (fServer) fServer->Register(s.dir.c_str(), s.hitpat);
                    if (fUseEvents) {
                        s.multi = new TH1F(("h_" + s.name + "_multi").c_str(),
                                           (s.title + " multiplicity;hits/event;counts").c_str(),
                                           s.nbins + 1, -0.5, s.nbins + 0.5);
                        if (fServer) fServer->Register(s.dir.c_str(), s.multi);
                    }
                    fSets.push_back(std::move(s));
                }
            }
        }
        if (fServer) fServer->SetItemField("/Detector", "_monitoring", "1000");
    }

    void CreateTimeHistograms(Set& s, const chmap::DeviceDef& d, const Ranges& r)
    {
        s.tdc.assign(s.nbins, nullptr);
        s.tot.assign(s.nbins, nullptr);
        for (int b = 0; b < s.nbins; ++b) {
            const std::string label = ChannelLabel(d, s, b);
            const std::string tag   = Sanitize(label);
            s.tdc[b] = new TH1F(("h_" + s.name + "_tdc_" + tag).c_str(),
                                (s.title + " " + label + " TDC;" + r.tdc_axis + ";counts").c_str(),
                                r.tdc_bins, r.tdc_min, r.tdc_max);
            s.tot[b] = new TH1F(("h_" + s.name + "_tot_" + tag).c_str(),
                                (s.title + " " + label + " TOT;TOT;counts").c_str(), 1000, 0, r.tot_max);
            if (fServer) {
                fServer->Register((s.dir + "/TDC").c_str(), s.tdc[b]);
                fServer->Register((s.dir + "/TOT").c_str(), s.tot[b]);
            }
        }
        for (int first = 0; first < s.nbins; first += 32) {
            const int last = std::min(first + 32, s.nbins) - 1;
            const std::string range = ChannelLabel(d, s, first) + " - " + ChannelLabel(d, s, last);
            for (const bool is_tot : {false, true}) {
                auto* c = new TCanvas(("c_" + s.name + (is_tot ? "_tot_" : "_tdc_") + std::to_string(first / 32)).c_str(),
                                      (s.title + (is_tot ? " TOT " : " TDC ") + range).c_str(), 1600, 1000);
                c->Divide(8, 4);
                for (int b = first; b <= last; ++b) {
                    c->cd(b - first + 1);
                    if (is_tot) {
                        if (auto* pad = c->GetPad(b - first + 1)) pad->SetLogy(1);
                        s.tot[b]->Draw();
                    } else {
                        s.tdc[b]->Draw();
                    }
                }
                s.canvases.push_back(c);
                if (fServer) fServer->Register((s.dir + "/Canvases").c_str(), c);
            }
        }
    }

    void FillUnmapped(uint32_t fem_id, uint32_t ch)
    {
        const std::string ip = chmap::DetectorMap::ip_to_string(fem_id);
        if (fUnmapped) fUnmapped->Fill(ch, ip.c_str(), 1.);
        if (fUnmappedSeen.insert((static_cast<uint64_t>(fem_id) << 32) | ch).second) {
            LOG(warn) << "FE " << ip << " ch " << ch << " is not in the channel map (reported once per run)";
        }
    }

    THttpServer* fServer = nullptr;
    std::string  fDefinitionPath, fMapPath;
    Ranges       fRangesHR, fRangesLR;
    bool         fUseEvents = false;

    chmap::DetectorMap fMap;
    bool  fLoaded = false;
    Stamp fStamp;

    std::vector<Set>         fSets;
    std::vector<std::size_t> fOffset;   // device index -> first set
    TH2F*                    fUnmapped = nullptr;
    std::set<uint64_t>       fUnmappedSeen;
};

#pragma pack(pop)

#endif
