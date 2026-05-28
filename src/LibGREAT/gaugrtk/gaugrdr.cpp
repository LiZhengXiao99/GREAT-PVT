/**
 * @file         gaugrdr.cpp
 * @brief        AUG file reader for PPP-RTK augmentation
 */
#include "gaugrdr.h"
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>

namespace great {

t_gaugrdr::t_gaugrdr(const string& path, const string& site,
                     const t_gtriple& xyz_rov, t_spdlog spdlog)
    : _path(path), _site(site), _xyz_rov(xyz_rov), _spdlog(spdlog),
      _dist_km(0.0), _is_open(false), _is_lsq(false), _sample(30.0)
{
}

t_gaugrdr::~t_gaugrdr()
{
}

bool t_gaugrdr::init()
{
    if (_is_open) return true;
    return _resolvePath();
}

bool t_gaugrdr::_resolvePath()
{
    struct stat st;
    if (stat(_path.c_str(), &st) != 0) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "AUG path invalid: " + _path);
        return false;
    }

    // File mode
    if (S_ISREG(st.st_mode)) {
        return _parseFile(_path);
    }

    // Directory mode
    if (!S_ISDIR(st.st_mode)) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "AUG path is neither file nor directory: " + _path);
        return false;
    }

    vector<AugHeadInfo> candidates;
    DIR* dir = opendir(_path.c_str());
    if (!dir) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "Cannot open AUG directory: " + _path);
        return false;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string fname = entry->d_name;
        if (fname.length() < 5) continue; // at least "x.aug"

        // Check extension
        string ext = fname.substr(fname.find_last_of('.') + 1);
        transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext != "aug") continue;

        // Full path
        string fpath = _path;
        if (fpath.back() != '/') fpath += '/';
        fpath += fname;

        // Verify it's a regular file
        struct stat fst;
        if (stat(fpath.c_str(), &fst) != 0 || !S_ISREG(fst.st_mode)) continue;

        AugHeadInfo info;
        info.filepath = fpath;

        // Read header only
        ifstream fs(info.filepath);
        if (!fs.is_open()) continue;

        string line;
        bool xyz_ok = false;
        while (getline(fs, line)) {
            if (line.find("END OF HEADER") != string::npos) break;
            if (line.find("ACCURATE POSITION XYZ") != string::npos) {
                istringstream iss(line);
                double x, y, z;
                if (iss >> x >> y >> z) {
                    info.xyz.set(0, x); info.xyz.set(1, y); info.xyz.set(2, z);
                    xyz_ok = true;
                }
            }
        }
        fs.close();

        if (!xyz_ok) continue;

        // Compute distance
        t_gtriple diff = info.xyz - _xyz_rov;
        info.dist_km = diff.norm() / 1000.0;

        // Get filename stem (first 4 chars as sitename)
        string stem = fname.substr(0, fname.find_last_of('.'));
        if (stem.size() >= 4) info.stname = stem.substr(0, 4);
        else info.stname = stem;
        transform(info.stname.begin(), info.stname.end(), info.stname.begin(), ::toupper);

        candidates.push_back(info);
    }
    closedir(dir);

    if (candidates.empty()) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "No .aug files found in: " + _path);
        return false;
    }

    // Priority 1: filename prefix match (first 4 chars = site name)
    string site4 = _site;
    if (site4.size() >= 4) site4 = site4.substr(0, 4);
    transform(site4.begin(), site4.end(), site4.begin(), ::toupper);

    vector<AugHeadInfo> matched;
    for (const auto& c : candidates) {
        if (c.stname == site4) matched.push_back(c);
    }

    const vector<AugHeadInfo>& pool = matched.empty() ? candidates : matched;

    // Priority 2: closest distance
    size_t idx_sel = 0;
    double dmin = 1e9;
    for (size_t i = 0; i < pool.size(); i++) {
        if (pool[i].dist_km < dmin) { dmin = pool[i].dist_km; idx_sel = i; }
    }

    if (_spdlog) {
        SPDLOG_LOGGER_INFO(_spdlog,
            "AUG selected: " + pool[idx_sel].filepath +
            " site=" + pool[idx_sel].stname +
            " dist=" + to_string(int(dmin)) + " km");
    }

    return _parseFile(pool[idx_sel].filepath);
}

bool t_gaugrdr::_parseFile(const string& filepath)
{
    ifstream fs(filepath);
    if (!fs.is_open()) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "Cannot open AUG: " + filepath);
        return false;
    }

    if (!_parseHeader(fs)) {
        fs.close();
        return false;
    }

    if (!_parseBody(fs)) {
        fs.close();
        return false;
    }

    fs.close();
    _is_open = true;
    return true;
}

bool t_gaugrdr::_parseHeader(ifstream& fs)
{
    string line;
    bool xyz_ok = false;

    while (getline(fs, line)) {
        if (line.find("END OF HEADER") != string::npos) break;

        if (line.find("EXECUTABLE NAME") != string::npos) {
            if (line.find("LSQ") != string::npos) _is_lsq = true;
            else _is_lsq = false;
        }
        else if (line.find("ACCURATE POSITION XYZ") != string::npos) {
            istringstream iss(line);
            double x, y, z;
            if (iss >> x >> y >> z) {
                _xyz.set(0, x); _xyz.set(1, y); _xyz.set(2, z);
                xyz_ok = true;
            }
        }
    }

    if (!xyz_ok) return false;

    t_gtriple diff = _xyz - _xyz_rov;
    _dist_km = diff.norm() / 1000.0;

    return true;
}

bool t_gaugrdr::_parseBody(ifstream& fs)
{
    string line;
    vector<AugEpoch> raw_epochs;
    AugEpoch cur_ep;
    bool has_epoch = false;

    // First pass: read all raw epochs and satellites (no filtering)
    while (getline(fs, line)) {
        if (line.empty() || line.size() < 10) continue;

        // Epoch line: > YYYY MM DD HH MM SS.ssssss  n
        if (line[0] == '>') {
            if (has_epoch && !cur_ep.sats.empty()) {
                raw_epochs.push_back(cur_ep);
            }

            istringstream iss(line.substr(1));
            int yr, mn, dy, hr, mi;
            double sec;
            int n_sat;
            iss >> yr >> mn >> dy >> hr >> mi >> sec >> n_sat;

            cur_ep = AugEpoch();
            cur_ep.epoch.from_ymdhms(yr, mn, dy, hr, mi, sec);
            cur_ep.sats.reserve(n_sat);
            has_epoch = true;

            // Read ZWD line
            bool zwd_ok = false;
            if (getline(fs, line)) {
                istringstream zwdiss(line);
                string zwd_tag;
                zwdiss >> zwd_tag >> cur_ep.zwd >> cur_ep.grd_n >> cur_ep.grd_e
                       >> cur_ep.zwd_std0 >> cur_ep.zwd_std1 >> cur_ep.zwd_std2
                       >> cur_ep.sys >> cur_ep.nsys;
                // Validate ZWD line: must have ZWD tag and positive std
                zwd_ok = (zwd_tag == "ZWD" && cur_ep.zwd_std0 > 0.0);
            }
            if (!zwd_ok) {
                if (_spdlog)
                    SPDLOG_LOGGER_WARN(_spdlog,
                        "AUG ZWD line invalid or missing at epoch " +
                        cur_ep.epoch.str_ymdhms() + ", skipping epoch");
                has_epoch = false;
                continue;
            }
            continue;
        }

        // Satellite line
        if (has_epoch) {
            AugSatRec rec;
            istringstream iss(line);
            iss >> rec.prn >> rec.ipp_lat >> rec.ipp_lon
                >> rec.az_deg >> rec.el_deg >> rec.stec_m >> rec.sigma_m
                >> rec.lck >> rec.el_max_deg >> rec.stec_p1p2
                >> rec.fix_ewl >> rec.fix_wl >> rec.fix_nl;

            if (!rec.prn.empty()) {
                cur_ep.sats.push_back(rec);
            }
        }
    }

    if (has_epoch && !cur_ep.sats.empty()) {
        raw_epochs.push_back(cur_ep);
    }

    if (raw_epochs.empty()) return false;

    // Infer sampling interval from first two valid epochs
    _sample = 30.0; // default fallback
    if (raw_epochs.size() >= 2) {
        double dt = fabs(raw_epochs[1].epoch.diff(raw_epochs[0].epoch));
        if (dt > 0.0 && dt < 3600.0) _sample = dt;
    }

    // Second pass: apply quality control per satellite and per epoch
    for (auto& ep : raw_epochs) {
        vector<AugSatRec> filtered;
        for (auto& rec : ep.sats) {
            if (_filterSat(rec)) filtered.push_back(rec);
        }

        // Keep epoch if it has valid satellites OR valid ZWD data.
        // Previously discarded if fewer than 3 sats, but ZWD-only epochs are still useful
        // for PPP-RTK when STEC constraints are scarce (e.g. early epochs with unfixed amb).
        bool has_zwd = (ep.zwd_std0 > 0.0 || ep.nsys > 0);
        if (!filtered.empty() || has_zwd) {
            ep.sats = filtered;
            _epochs.push_back(ep);
        }
    }

    return !_epochs.empty();
}

bool t_gaugrdr::_filterSat(AugSatRec& rec)
{
    // Skip GLONASS
    if (!rec.prn.empty() && rec.prn[0] == 'R') return false;

    // Skip satellites with all-zero fix tags
    if (rec.fix_ewl == 0 && rec.fix_wl == 0 && rec.fix_nl == 0) return false;

    if (_is_lsq) {
        // Strict thresholds for LSQ-generated AUG
        if (rec.lck < 60) return false;
        if (rec.el_max_deg < 20.0) return false;
        // lck * az_deg >= 900 (az_deg in degrees)
        if (rec.lck * rec.az_deg < 900.0) return false;
        // lck * sample >= 900 s (total arc length)
        if (rec.lck * _sample < 900.0) return false;
        // lck * sample * az_deg >= 6000
        if (rec.lck * _sample * rec.az_deg < 6000.0) return false;
    } else {
        // Relaxed thresholds for UDUC/KF-generated AUG
        if (rec.lck < 3) return false;
        if (rec.el_max_deg < 12.0) return false;
        // lck * sample >= 15 s (non-single-epoch mode)
        if (rec.lck * _sample < 15.0) return false;
    }

    // Fix flag check: NL fixed OR (WL fixed AND EWL fixed)
    // Only satellites with reliable ambiguity fix are used for ionospheric constraints.
    if (rec.fix_nl != 1) {
        if (rec.fix_wl != 1) return false;
        else if (rec.fix_ewl != 1) return false;
    }

    return true;
}

bool t_gaugrdr::getEpoch(const t_gtime& t, AugEpoch& out)
{
    if (!_is_open || _epochs.empty()) return false;

    // Find closest epoch
    double dt_min = 1e9;
    int idx = -1;
    for (size_t i = 0; i < _epochs.size(); i++) {
        double dt = fabs(_epochs[i].epoch - t);
        if (dt < dt_min) {
            dt_min = dt;
            idx = i;
        }
    }

    if (idx < 0) return false;

    // Allow max 1.5 * sampling difference
    if (dt_min > _sample * 1.5) {
        if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog,
            "getEpoch: closest dt=" + to_string(dt_min) +
            " aug=" + _epochs[idx].epoch.str_ymdhms() +
            " req=" + t.str_ymdhms());
        return false;
    }

    out = _epochs[idx];
    return true;
}

} // namespace great
