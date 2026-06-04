/**
 * @file         gaugwriter.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        AUG file writer for PPP-AR augmentation
 * @version      1.0
 * @date         2024-08-29
 * 
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 * 
 */
#include "gproc/gaugwriter.h"
#include "gutils/gsysconv.h"
#include "gutils/gcommon.h"
#include "gutils/gfileconv.h"

using namespace great;

t_gaugwriter::t_gaugwriter(t_gsetbase *gset, t_spdlog spdlog)
    : _fout(nullptr), _header_written(false), _spdlog(spdlog)
{
    t_gsetout *gsetout = dynamic_cast<t_gsetout *>(gset);
    if (gsetout)
    {
        // Prefer explicit aug output path; fall back to ppp path derivation
        string aug_out = gsetout->outputs("aug");
        if (!aug_out.empty())
        {
            _path = aug_out;
            substitute(_path, GFILE_PREFIX, "");
        }
        else
        {
            string ppp_out = gsetout->outputs("ppp");
            if (!ppp_out.empty())
            {
                string tmp = ppp_out;
                substitute(tmp, GFILE_PREFIX, "");
                size_t last_sep = tmp.find_last_of("/\\");
                size_t last_dot = tmp.find_last_of('.');
                if (last_dot != string::npos && (last_sep == string::npos || last_dot > last_sep))
                {
                    _path = tmp.substr(0, last_dot) + ".aug";
                }
                else
                {
                    _path = tmp + ".aug";
                }
            }
        }
    }
}

t_gaugwriter::~t_gaugwriter()
{
    close();
}

void t_gaugwriter::setSite(const string& site)
{
    _site = site;
    if (!_path.empty())
    {
        substitute(_path, "$(rec)", _site, false);
    }
}

void t_gaugwriter::setPath(const string& path)
{
    if (!path.empty())
    {
        _path = path;
        substitute(_path, GFILE_PREFIX, "");
    }
}

bool t_gaugwriter::setHeader(const t_gtriple& xyz, const t_gtriple& blh,
                             const string& rcv_type, const string& ant_type, const string& ant_radome,
                             const t_gtime& beg, const t_gtime& end,
                             const set<string>& sys)
{
    _xyz = xyz;
    _blh = blh;
    _rcv = rcv_type;
    _ant = ant_type;
    _radome = ant_radome;
    _beg = beg;
    _end = end;

    ostringstream oss;
    for (auto it = sys.begin(); it != sys.end(); ++it)
    {
        if (it != sys.begin()) oss << "  ";
        oss << *it;
    }
    _sys_str = oss.str();

    return true;
}

bool t_gaugwriter::_writeHeader()
{
    if (_header_written) return true;
    if (_path.empty()) return false;

    if (!_fout)
    {
        _fout = new t_giof(_path);
    }

    auto fmt_comment = [](const string& content, const string& comment) -> string {
        string line = content;
        size_t end = line.find_last_not_of(' ');
        if (end != string::npos) line.erase(end + 1);
        int pad = 60 - (int)line.size();
        if (pad < 1) pad = 1;
        line.append(pad, ' ');
        line += "# " + comment + "\n";
        return line;
    };

    ostringstream os;
    os << fmt_comment("SAT  IPPLAT  IPPLON  AZIM  ELEV  STEC(m) SIGMA(m)", "DATA TYPES");
    os << fmt_comment("     NLOCK  EL_Max  StecP1P2  FIX_TAG", "DATA TYPES");
    os << fmt_comment("SIon_given = SIon_true - 1/((F1/F2)**2-1) * (DCB_r + DCB_s)", "COMMENT");
    os << fmt_comment("63781363", "R OF EARTH (cm)");
    os << fmt_comment("450000", "HEIGHT OF SHELL (cm)");

    // TIME START
    {
        ostringstream tmp;
        int wk = _beg.gwk();
        double sow = _beg.sow() + _beg.dsec();
        tmp << " (" << setw(4) << setfill('0') << wk << " " << fixed << setprecision(1) << setw(10) << sow << "s) ";
        tmp << setfill(' ');
        tmp << _beg.str_ymd() << " " << _beg.str_hms() << ".0 GPST";
        os << fmt_comment(tmp.str(), "TIME START");
    }

    // TIME END
    {
        ostringstream tmp;
        int wk = _end.gwk();
        double sow = _end.sow() + _end.dsec();
        tmp << " (" << setw(4) << setfill('0') << wk << " " << fixed << setprecision(1) << setw(10) << sow << "s) ";
        tmp << setfill(' ');
        tmp << _end.str_ymd() << " " << _end.str_hms() << ".0 GPST";
        os << fmt_comment(tmp.str(), "TIME END");
    }

    // RCV/ANT TYPE
    {
        ostringstream tmp;
        tmp << left << setw(30) << _rcv << " " << setw(15) << _ant << " " << setw(6) << _radome;
        os << fmt_comment(tmp.str(), "RCV/ANT TYPE");
    }

    // ACCURATE POSITION XYZ (retain 4 decimals => 0.1 mm)
    {
        ostringstream tmp;
        tmp << fixed << setprecision(4)
            << setw(16) << _xyz[0] << " " << setw(16) << _xyz[1] << " " << setw(16) << _xyz[2];
        os << fmt_comment(tmp.str(), "ACCURATE POSITION XYZ");
    }

    // ACCURATE POSITION BLH (lat/lon in degrees, height in m)
    {
        ostringstream tmp;
        tmp << fixed << setprecision(10) << setw(18) << _blh[0] * R2D << " ";
        tmp << fixed << setprecision(10) << setw(18) << _blh[1] * R2D << " ";
        tmp << fixed << setprecision(4) << setw(16) << _blh[2];
        os << fmt_comment(tmp.str(), "ACCURATE POSITION BLH");
    }

    // SAT SYSTEM
    {
        ostringstream tmp;
        tmp << " " << left << setw(59) << _sys_str;
        os << fmt_comment(tmp.str(), "SAT SYSTEM");
    }

    // EXECUTABLE NAME
    os << fmt_comment(" GREAT-PVT", "EXECUTABLE NAME");

    // END OF HEADER
    os << fmt_comment("", "END OF HEADER");

    _fout->write(os.str().c_str(), os.str().size());
    _fout->flush();
    _header_written = true;

    return true;
}

bool t_gaugwriter::_writeObsCodePairs(const vector<t_gsatdata>& data)
{
    if (!_fout) return false;

    // Collect one obs-code pair per system from the first available satellite
    map<char, pair<string, string>> sys_pairs;

    for (const auto& satdata : data)
    {
        char sys_char = satdata.sat()[0];
        if (sys_pairs.find(sys_char) != sys_pairs.end()) continue;

        set<GOBSBAND> bands = satdata.band_avail(true);
        if (bands.size() < 2) continue;

        // Find first two valid phase observations
        GOBS obs1 = X;
        GOBS obs2 = X;
        auto it = bands.begin();
        while (it != bands.end() && obs1 == X)
        {
            obs1 = satdata.id_phase(*it);
            ++it;
        }
        while (it != bands.end() && obs2 == X)
        {
            obs2 = satdata.id_phase(*it);
            ++it;
        }

        if (obs1 != X && obs2 != X)
        {
            sys_pairs[sys_char] = make_pair(gobs2str(obs1), gobs2str(obs2));
        }
    }

    if (sys_pairs.empty()) return true;

    ostringstream os;
    auto sys_name = [](char c) -> string {
        switch (c) {
            case 'G': return "GPS";
            case 'R': return "GLO";
            case 'E': return "GAL";
            case 'C': return "BDS";
            case 'J': return "QZS";
            case 'I': return "IRN";
            case 'S': return "SBAS";
            default:  return "UNK";
        }
    };

    for (const auto& sp : sys_pairs)
    {
        const string& s1 = sp.second.first;
        const string& s2 = sp.second.second;
        string pair_str = s1 + "-" + s2;
        string name = sys_name(sp.first);

        os << "  " << name << ": " << pair_str;
        int prefix_len = 2 + (int)name.size() + 2 + (int)pair_str.size();
        int padding = 60 - prefix_len;
        if (padding < 1) padding = 1;
        os << string(padding, ' ') << "# OBS CODE PAIR\n";
    }

    _fout->write(os.str().c_str(), os.str().size());
    _fout->flush();
    return true;
}

bool t_gaugwriter::writeEpoch(const t_gtime& epoch,
                              const t_gallpar& X_fix,
                              const t_gallpar& X_float,
                              const SymmetricMatrix& Q_fix,
                              const vector<t_gsatdata>& data,
                              const map<string, double>& ele,
                              t_gambiguity* ambfix,
                              const map<string, int>& lock_epo,
                              const map<string, double>& el_max,
                              const string& site)
{
    if (!_header_written)
    {
        if (!_writeHeader())
            return false;
        _writeObsCodePairs(data);
    }

    if (!_fout) return false;

    return _writeEpochBody(epoch, X_fix, X_float, Q_fix, data, ele, ambfix, lock_epo, el_max, site);
}

bool t_gaugwriter::_writeEpochBody(const t_gtime& epoch,
                                   const t_gallpar& X_fix,
                                   const t_gallpar& X_float,
                                   const SymmetricMatrix& Q_fix,
                                   const vector<t_gsatdata>& data,
                                   const map<string, double>& ele,
                                   t_gambiguity* ambfix,
                                   const map<string, int>& lock_epo,
                                   const map<string, double>& el_max,
                                   const string& site)
{
    // 1. Get coordinates
    t_gtriple xyz, blh;
    X_fix.getCrdParam(site, xyz);
    xyz2ell(xyz, blh, false);

    // 2. Get ZWD and gradients
    int itrp = X_fix.getParam(site, par_type::TRP, "");
    int igrn = X_fix.getParam(site, par_type::GRD_N, "");
    int igre = X_fix.getParam(site, par_type::GRD_E, "");
    double zwd = (itrp >= 0) ? X_fix.getPar(itrp).value() : 0.0;
    double zwd_std = (itrp >= 0) ? sqrt(Q_fix(itrp + 1, itrp + 1)) : 0.5;
    double grd_n = (igrn >= 0) ? X_fix.getPar(igrn).value() : 0.0;
    double grd_e = (igre >= 0) ? X_fix.getPar(igre).value() : 0.0;

    // 3. Collect per-satellite records
    struct AugSatRec
    {
        string prn;
        double ipp_lat;
        double ipp_lon;
        double azim;
        double elev;
        double stec;
        double sigma;
        int lck;
        double el_max;
        double stec_p1p2;
        int fix_EWL;
        int fix_WL;
        int fix_NL;
    };

    vector<AugSatRec> recs;

    // Pre-collect DD ambiguity fix states
    map<string, tuple<int, int, int>> sat_fix; // prn -> (EWL, WL, NL)
    if (ambfix)
    {
        const t_DD_amb& dd = ambfix->getDD();
        for (const auto& itdd : dd)
        {
            for (const auto& s : itdd.ddSats)
            {
                string prn = get<0>(s);
                auto it = sat_fix.find(prn);
                if (it == sat_fix.end())
                {
                    sat_fix[prn] = make_tuple(0, 0, 0);
                    it = sat_fix.find(prn);
                }
                if (itdd.isEwlFixed) get<0>(it->second) = 1;
                if (itdd.isWlFixed)  get<1>(it->second) = 1;
                if (itdd.isNlFixed)  get<2>(it->second) = 1;
            }
        }
    }

    for (const auto& satdata : data)
    {
        string prn = satdata.sat();

        // Skip GLONASS
        if (prn[0] == 'R') continue;

        // 3.1 Ambiguity fix flags (EWL / WL / NL)
        int fix_EWL = 0, fix_WL = 0, fix_NL = 0;
        auto itf = sat_fix.find(prn);
        if (itf != sat_fix.end())
        {
            fix_EWL = get<0>(itf->second);
            fix_WL  = get<1>(itf->second);
            fix_NL  = get<2>(itf->second);
        }

        // 3.2 Ionosphere (SION)
        // Use fixed-solution STEC if EWL/WL/NL is fixed (X_fix = _param_aug),
        // otherwise use float-solution STEC.
        int ision = X_fix.getParam(site, par_type::SION, prn);
        if (ision < 0) continue; // no SION parameter
        
        double stec, stec_std;
        // EWL / WL / NL any fixed -> use X_fix (_param_aug, containing constraints)
        bool is_any_fixed = (fix_NL == 1 || fix_WL == 1 || fix_EWL == 1);
        if (is_any_fixed)
        {
            stec = X_fix.getPar(ision).value();
            stec_std = sqrt(Q_fix(ision + 1, ision + 1));
        }
        else
        {
            int ision_float = X_float.getParam(site, par_type::SION, prn);
            if (ision_float >= 0)
            {
                stec = X_float.getPar(ision_float).value();
                stec_std = sqrt(Q_fix(ision_float + 1, ision_float + 1));
            }
            else
            {
                stec = X_fix.getPar(ision).value();
                stec_std = sqrt(Q_fix(ision + 1, ision + 1));
            }
        }

        // 3.3 Geometry
        double az = satdata.azi() * R2D;
        double el = satdata.ele_deg();

        // Compute IPP
        t_gtriple ipp_ell;
        t_gsatdata tmp_sat = satdata;
        ell2ipp(tmp_sat, blh, ipp_ell);
        double ipp_lat = ipp_ell[0] * R2D;
        double ipp_lon = ipp_ell[1] * R2D;

        // 3.5 Lock epochs and max elevation
        int lck = 0;
        auto itlck = lock_epo.find(prn);
        if (itlck != lock_epo.end()) lck = itlck->second;

        double elmax = 0.0;
        auto itel = el_max.find(prn);
        if (itel != el_max.end()) elmax = itel->second;

        recs.push_back({prn, ipp_lat, ipp_lon, az, el, stec, stec_std,
                        lck, elmax, 0.0, fix_EWL, fix_WL, fix_NL});
    }

    // Sort by system order (GREC) and PRN number
    auto sys_order = [](char sys) -> int {
        switch (sys) {
            case 'G': return 0; // GPS
            case 'R': return 1; // GLO
            case 'E': return 2; // GAL
            case 'C': return 3; // BDS
            case 'J': return 4; // QZSS
            case 'I': return 5; // IRNSS
            case 'S': return 6; // SBAS
            default:  return 7;
        }
    };
    sort(recs.begin(), recs.end(), [&](const AugSatRec& a, const AugSatRec& b) {
        int oa = sys_order(a.prn[0]);
        int ob = sys_order(b.prn[0]);
        if (oa != ob) return oa < ob;
        int na = atoi(a.prn.substr(1).c_str());
        int nb = atoi(b.prn.substr(1).c_str());
        return na < nb;
    });

    if (recs.size() < 1) return false; // no sats to output

    // 4. Format output
    ostringstream os;
    int n_total = 1 + recs.size(); // ZWD line + sat lines
    os << "> " << epoch.year() << " " << epoch.mon() << " " << epoch.day()
       << " " << epoch.hour() << " " << epoch.mins() << " "
       << fixed << setprecision(6) << setw(12) << epoch.dsec() + (epoch.sow() % 60)
       << " " << setw(3) << n_total << "\n";

    os << "ZWD " << fixed << setprecision(4)
       << setw(12) << zwd << setw(12) << grd_n << setw(12) << grd_e
       << setw(12) << zwd_std << "    0.0000    0.0000   9  1\n";

    for (const auto& r : recs)
    {
        os << left << setw(4) << r.prn << " "
           << fixed << setprecision(6)
           << setw(12) << r.ipp_lat << " " << setw(12) << r.ipp_lon << " "
           << fixed << setprecision(4)
           << setw(11) << r.azim << " " << setw(11) << r.elev << " "
           << setw(11) << r.stec << " "
           << fixed << setprecision(4) << setw(10) << r.sigma << " "
           << setw(6) << r.lck << " "
           << fixed << setprecision(4) << setw(11) << r.el_max << " "
           << setw(11) << r.stec_p1p2 << " "
           << r.fix_EWL << " " << r.fix_WL << " " << r.fix_NL << "\n";
    }

    _fout->write(os.str().c_str(), os.str().size());
    _fout->flush();
    return true;
}

void t_gaugwriter::close()
{
    if (_fout)
    {
        _fout->close();
        delete _fout;
        _fout = nullptr;
    }
}
