/**
 * @file         t_gppprtk.cpp
 * @brief        PPP-RTK constraint builder for single-station AUG
 * @author       Li Zhengxiao
 */
#include "t_gppprtk.h"
#include "gutils/gtypeconv.h"
#include <algorithm>
#include <cmath>

namespace great {

t_gppprtk::t_gppprtk(const string& path, double ionoCfg, double tropCfg,
                       t_spdlog spdlog)
    : _rdr(nullptr), _spdlog(spdlog), _enabled(false), _ionoCfg(0.0), _tropCfg(0.0), _outlierThres(0.0), _iqrEnabled(false), _stecBoostN(0), _stecBoostM(1.0), _epochCount(0)
{
    if (path.empty()) {
        if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, "AUG path empty. PPP-RTK disabled.");
        return;
    }

    string p = path;
    // Trim whitespace
    auto trim = [](string& s) {
        size_t first = s.find_first_not_of(" \t\n\r");
        if (first == string::npos) { s.clear(); return; }
        size_t last = s.find_last_not_of(" \t\n\r");
        s = s.substr(first, last - first + 1);
    };
    trim(p);

    if (p.empty()) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "AUG path empty after trim. PPP-RTK disabled.");
        return;
    }

    _path = p;
    _ionoCfg = ionoCfg;
    _tropCfg = tropCfg;

    if (_ionoCfg == 0.0 && _tropCfg == 0.0) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "Both iono and trop are 0. PPP-RTK disabled.");
        return;
    }

    _enabled = true;

    if (_spdlog) {
        SPDLOG_LOGGER_INFO(_spdlog,
            "PPP-RTK config: path=" + _path +
            " iono=" + to_string(_ionoCfg) +
            " trop=" + to_string(_tropCfg));
    }
}

t_gppprtk::~t_gppprtk()
{
    if (_rdr) {
        delete _rdr;
        _rdr = nullptr;
    }
}

bool t_gppprtk::initReader(const string& site, const t_gtriple& xyz_rov)
{
    if (!_enabled || _path.empty()) return false;
    if (_rdr) return true; // already initialized

    _rdr = new t_gaugrdr(_path, site, xyz_rov, _spdlog);
    if (!_rdr->init()) {
        if (_spdlog) SPDLOG_LOGGER_WARN(_spdlog, "PPP-RTK AUG reader init failed for path: " + _path);
        delete _rdr;
        _rdr = nullptr;
        return false;
    }
    if (_spdlog) SPDLOG_LOGGER_INFO(_spdlog, "PPP-RTK AUG reader initialized: " + _path);
    return true;
}

bool t_gppprtk::queryZwd(const t_gtime& t, double& zwd, double& zwd_std,
                         double max_dt) const
{
    if (!_enabled || !_rdr) return false;

    AugEpoch aug;
    if (!_rdr->getNearestEpoch(t, aug, max_dt)) return false;

    if (aug.zwd_std0 <= 0.0) return false;

    zwd = aug.zwd;
    zwd_std = aug.zwd_std0;
    return true;
}

bool t_gppprtk::buildConstraints(const t_gtime& epoch,
                                  t_gallpar& param_float,
                                  const vector<t_gsatdata>& data,
                                  const ColumnVector& dx,
                                  const map<string, int>& lock_epo,
                                  const SymmetricMatrix& Qx,
                                  Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                                  int& n_stec, int& n_zwd)
{
    n_stec = 0; n_zwd = 0;
    A_aug.ReSize(0, 0);
    P_aug.ReSize(0);
    l_aug.ReSize(0);

    if (!_enabled || !_rdr) return false;

    AugEpoch aug;
    if (!_rdr->getEpoch(epoch, aug)) return false;

    int nPar = param_float.parNumber();
    int n_aug = 0;

    Matrix A_stec, A_zwd;
    DiagonalMatrix P_stec, P_zwd;
    ColumnVector l_stec, l_zwd;
    int n_rej_stec = 0, n_rej_zwd = 0;

    // Build STEC constraints
    if (_ionoCfg != 0.0) {
        _buildStecConstraints(aug, param_float, data, dx, lock_epo, A_stec, P_stec, l_stec, n_stec, n_rej_stec);
    }

    // Count successful AUG epochs for STEC boost
    if (_stecBoostN > 0 && n_stec > 0) {
        _epochCount++;
    }

    // Build ZWD constraints
    if (_tropCfg != 0.0) {
        _buildZwdConstraints(aug, param_float, dx, Qx, A_zwd, P_zwd, l_zwd, n_zwd, n_rej_zwd);
    }

    n_aug = n_stec + n_zwd;
    if (n_aug == 0) return false;

    // Combine
    A_aug.ReSize(n_aug, nPar); A_aug = 0.0;
    P_aug.ReSize(n_aug);       P_aug = 0.0;
    l_aug.ReSize(n_aug);       l_aug = 0.0;

    for (int i = 1; i <= n_stec; i++) {
        for (int j = 1; j <= nPar; j++) A_aug(i, j) = A_stec(i, j);
        P_aug(i) = P_stec(i);
        l_aug(i) = l_stec(i);
    }
    for (int i = 1; i <= n_zwd; i++) {
        int row = n_stec + i;
        for (int j = 1; j <= nPar; j++) A_aug(row, j) = A_zwd(i, j);
        P_aug(row) = P_zwd(i);
        l_aug(row) = l_zwd(i);
    }

    // Statistics output
    if (_spdlog) {
        auto stats = [](const ColumnVector& v) -> string {
            int n = v.Nrows();
            if (n == 0) return "n=0";
            vector<double> vals;
            vals.reserve(n);
            double sum = 0.0, sum2 = 0.0;
            double vmin = v(1), vmax = v(1);
            for (int i = 1; i <= n; i++) {
                double x = v(i);
                vals.push_back(x);
                sum += x;
                sum2 += x * x;
                if (x < vmin) vmin = x;
                if (x > vmax) vmax = x;
            }
            double mean = sum / n;
            double variance = sum2 / n - mean * mean;
            if (variance < 0.0 && variance > -1e-12) variance = 0.0;
            double stddev = sqrt(variance);
            sort(vals.begin(), vals.end());
            double p50 = vals[n / 2];
            double p16 = vals[n * 16 / 100];
            double p84 = vals[n * 84 / 100];
            return "n=" + to_string(n) + " mean=" + to_string(mean) + "m std=" + to_string(stddev)
                 + "m min=" + to_string(vmin) + "m max=" + to_string(vmax) + "m"
                 + " p16=" + to_string(p16) + "m p50=" + to_string(p50) + "m p84=" + to_string(p84) + "m";
        };

        string msg = "PPP-RTK constraints: STEC=" + to_string(n_stec) + " ZWD=" + to_string(n_zwd);
        if (n_rej_stec > 0 || n_rej_zwd > 0)
            msg += " rejected(STEC=" + to_string(n_rej_stec) + " ZWD=" + to_string(n_rej_zwd) + ")";
        SPDLOG_LOGGER_INFO(_spdlog, msg);

        if (n_stec > 0) {
            SPDLOG_LOGGER_INFO(_spdlog, "  STEC residual stats: " + stats(l_stec));
        }
        if (n_zwd > 0) {
            SPDLOG_LOGGER_INFO(_spdlog, "  ZWD  residual stats: " + stats(l_zwd));
        }
    }

    return true;
}

string t_gppprtk::_selectRefSat(const AugEpoch& aug, const vector<t_gsatdata>& data,
                                     char sys, const map<string, int>& lock_epo)
{
    // Select reference sat by max(lck * el) every epoch (no persistence)
    double max_lck_el = 0.0;
    string ref;

    for (const auto& sat : data) {
        string prn = sat.sat();
        if (prn.empty() || prn[0] != sys) continue; // only specified system
        if (prn[0] == 'R') continue; // skip GLONASS

        // Check if this sat exists in AUG
        const AugSatRec* rec = nullptr;
        for (const auto& a : aug.sats) {
            if (a.prn == prn) { rec = &a; break; }
        }
        if (!rec) continue;

        // Use rover current elevation and lock epochs
        int lck = 0;
        auto itlck = lock_epo.find(prn);
        if (itlck != lock_epo.end()) lck = itlck->second;

        double el_deg = sat.ele_deg();
        double lck_el = (lck > 0) ? (el_deg * lck) : el_deg;
        if (lck_el > max_lck_el) {
            max_lck_el = lck_el;
            ref = prn;
        }
    }

    // Cache selected ref for next epoch
    if (!ref.empty()) {
        _last_ref[sys] = ref;
    } else {
        _last_ref.erase(sys);
    }
    return ref;
}

bool t_gppprtk::_buildStecConstraints(const AugEpoch& aug,
                                       t_gallpar& param,
                                       const vector<t_gsatdata>& data,
                                       const ColumnVector& dx,
                                       const map<string, int>& lock_epo,
                                       Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                                       int& n_stec, int& n_rejected)
{
    n_stec = 0; n_rejected = 0;
    int nPar = param.parNumber();
    vector<int> idx_list;
    vector<double> v_list, p_list;

    // Process each system independently (GPS, GAL, BDS, QZS, IRN)
    const char syss[] = {'G', 'E', 'C', 'J', 'I'};
    const int nsys = sizeof(syss) / sizeof(syss[0]);

    for (int isys = 0; isys < nsys; isys++) {
        char sys = syss[isys];
        string ref = _selectRefSat(aug, data, sys, lock_epo);
        if (ref.empty()) continue;

        // Find reference sat in AUG
        const AugSatRec* ref_rec = nullptr;
        for (const auto& a : aug.sats) { if (a.prn == ref) { ref_rec = &a; break; } }
        if (!ref_rec) continue;

        int idx_rs = param.getParam(_site, par_type::SION, ref);
        if (idx_rs < 0) continue;

        for (const auto& sat : data) {
            string prn = sat.sat();
            if (prn == ref) continue;
            if (prn.empty() || prn[0] != sys) continue; // only same system
            if (prn[0] == 'R') continue;

            const AugSatRec* rec = nullptr;
            for (const auto& a : aug.sats) { if (a.prn == prn) { rec = &a; break; } }
            if (!rec) continue;

            int idx_s0 = param.getParam(_site, par_type::SION, prn);
            if (idx_s0 < 0) continue;

            // SD STEC residual: v = (STEC_s0 - STEC_rs)_aug - (x_post_s0 - x_post_rs)
            // param_float is the apriori state and dx is the PPP update vector.
            // Posterior state = param + dx (post-PPP float solution).
            double ion_aug_sd = rec->stec_m - ref_rec->stec_m;
            double ion_est_sd = (param[idx_s0].value() + dx(param[idx_s0].index))
                              - (param[idx_rs].value() + dx(param[idx_rs].index));
            double v = ion_aug_sd - ion_est_sd;

            // Outlier rejection
            if (_outlierThres > 0.0 && fabs(v) > _outlierThres) {
                n_rejected++;
                if (_spdlog)
                    SPDLOG_LOGGER_WARN(_spdlog,
                        "PPP-RTK STEC outlier rejected: " + ref + "->" + prn +
                        " aug=" + to_string(ion_aug_sd) + "m est=" + to_string(ion_est_sd) +
                        "m v=" + to_string(v) + "m (thres=" + to_string(_outlierThres) + "m)");
                continue;
            }

            // Variance: SD STEC variance = var(s0) + var(rs)
            double var_s0_base = _resolveVar(_ionoCfg, rec->sigma_m);
            double var_rs_base = _resolveVar(_ionoCfg, ref_rec->sigma_m);
            double var;
            if (_ionoCfg < 0.0) {
                // Sigma mode (auto sigma): direct variance sum, no elevation weighting.
                // var = var_s0_base + var_rs_base;
                var = var_s0_base;
            } else {
                // User-config mode: apply elevation-dependent weighting per satellite, then sum.
                double el_s0 = rec->el_deg * G_PI / 180.0;
                double el_rs = ref_rec->el_deg * G_PI / 180.0;
                var = _calIonoVar(el_s0, var_s0_base) + _calIonoVar(el_rs, var_rs_base);
            }

            // Boost STEC weight for first N epochs
            if (_stecBoostN > 0 && _epochCount < _stecBoostN) {
                var /= (_stecBoostM * _stecBoostM);
            }

            // NL un-fixed penalty removed per user request.
            // AUG STEC now comes from fixed-solution for fixed satellites
            // and float-solution for unfixed ones, so NL status no longer
            // indicates STEC quality.
            // if (!(rec->fix_nl == 1 && ref_rec->fix_nl == 1))
            //     var += 0.02 * 0.02;

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog,
                    "PPP-RTK STEC: " + ref + "->" + prn +
                    " aug=" + to_string(ion_aug_sd) + "m est=" + to_string(ion_est_sd) +
                    "m v=" + to_string(v) + "m var=" + to_string(var));

            idx_list.push_back(idx_s0);
            idx_list.push_back(idx_rs);
            v_list.push_back(v);
            p_list.push_back(1.0 / var);
            n_stec++;
        }
    }

    // IQR-based outlier rejection for STEC constraints
    // Minimum sample size raised to 10 to avoid degeneracy at n=4 (IQR=0 rejects all non-equal values).
    if (_iqrEnabled && n_stec >= 10) {
        vector<double> sorted_v = v_list;
        sort(sorted_v.begin(), sorted_v.end());
        int n = sorted_v.size();
        double Q1 = sorted_v[n / 4];
        double Q3 = sorted_v[3 * n / 4];
        double IQR = Q3 - Q1;
        // Protect against IQR collapse (repeated values) which would reject everything.
        if (IQR < 1e-6) {
            if (_spdlog)
                SPDLOG_LOGGER_INFO(_spdlog,
                    "PPP-RTK STEC IQR skipped: IQR=" + to_string(IQR) +
                    "m too small (n=" + to_string(n) + ")");
        } else {
            double lower = Q1 - 1.5 * IQR;
            double upper = Q3 + 1.5 * IQR;

        vector<int> new_idx_list;
        vector<double> new_v_list, new_p_list;
        int new_n_stec = 0;
            for (int i = 0; i < n_stec; i++) {
                double v = v_list[i];
                if (v < lower || v > upper) {
                    n_rejected++;
                    if (_spdlog)
                        SPDLOG_LOGGER_INFO(_spdlog,
                            "PPP-RTK STEC IQR outlier rejected: v=" + to_string(v) +
                            "m Q1=" + to_string(Q1) + "m Q3=" + to_string(Q3) +
                            "m IQR=" + to_string(IQR) + "m bounds=[" +
                            to_string(lower) + "," + to_string(upper) + "]");
                    continue;
                }
                new_idx_list.push_back(idx_list[i * 2]);
                new_idx_list.push_back(idx_list[i * 2 + 1]);
                new_v_list.push_back(v);
                new_p_list.push_back(p_list[i]);
                new_n_stec++;
            }
            idx_list = new_idx_list;
            v_list = new_v_list;
            p_list = new_p_list;
            n_stec = new_n_stec;
        }
    }

    if (n_stec == 0) return false;

    A_aug.ReSize(n_stec, nPar); A_aug = 0.0;
    P_aug.ReSize(n_stec);       P_aug = 0.0;
    l_aug.ReSize(n_stec);       l_aug = 0.0;

    for (int i = 0; i < n_stec; i++) {
        int idx_s0 = idx_list[i * 2];
        int idx_rs = idx_list[i * 2 + 1];
        A_aug(i + 1, param[idx_s0].index) = 1.0;      // s0
        A_aug(i + 1, param[idx_rs].index) = -1.0;     // rs
        l_aug(i + 1) = v_list[i];
        P_aug(i + 1) = p_list[i];
    }
    return true;
}

bool t_gppprtk::_buildZwdConstraints(const AugEpoch& aug,
                                      t_gallpar& param,
                                      const ColumnVector& dx,
                                      const SymmetricMatrix& Qx,
                                      Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                                      int& n_zwd, int& n_rejected)
{
    n_zwd = 0; n_rejected = 0;
    int idx_zwd = param.getParam(_site, par_type::TRP, "");
    if (idx_zwd < 0) return false;

    // Posterior state = param + dx (post-PPP float solution).
    double est_zwd = param[idx_zwd].value() + dx(param[idx_zwd].index);
    double v = aug.zwd - est_zwd;

    // Outlier rejection
    if (_outlierThres > 0.0 && fabs(v) > _outlierThres) {
        n_rejected++;
        if (_spdlog)
            SPDLOG_LOGGER_WARN(_spdlog,
                "PPP-RTK ZWD outlier rejected: aug=" + to_string(aug.zwd) +
                "m est=" + to_string(est_zwd) + "m v=" + to_string(v) +
                "m (thres=" + to_string(_outlierThres) + "m)");
        return false;
    }

    double var;
    if (_tropCfg > 0.0) {
        var = _tropCfg * _tropCfg;
    } else if (_tropCfg < 0.0) {
        double s = fabs(_tropCfg) * aug.zwd_std0;
        if (s < 1e-6) s = 1e-6; // fallback to avoid extreme small variance
        var = s * s;
    } else {
        return false; // ZWD disabled when _tropCfg == 0
    }

    if (_spdlog)
        SPDLOG_LOGGER_DEBUG(_spdlog,
            "PPP-RTK ZWD: aug=" + to_string(aug.zwd) + "m est=" + to_string(est_zwd) +
            "m v=" + to_string(v) + "m var=" + to_string(var));

    int nPar = param.parNumber();
    A_aug.ReSize(1, nPar); A_aug = 0.0;
    P_aug.ReSize(1);       P_aug = 0.0;
    l_aug.ReSize(1);       l_aug = 0.0;

    A_aug(1, param[idx_zwd].index) = 1.0;
    l_aug(1) = v;
    P_aug(1) = 1.0 / var;
    n_zwd = 1;
    return true;
}

double t_gppprtk::_calIonoVar(double el_rad, double var0)
{
    double s = sin(el_rad);
    if (s < 1e-6) s = 1e-6;

    return var0 / (s * s);
}

double t_gppprtk::_resolveVar(double cfgVal, double augSigma)
{
    if (cfgVal > 0.0) {
        return cfgVal * cfgVal;
    } else if (cfgVal < 0.0) {
        double s = fabs(cfgVal) * augSigma;
        return s * s;
    }
    return 999999.0;
}

} // namespace great
