/**
 * @file         t_gppprtk.h
 * @brief        PPP-RTK constraint builder for single-station AUG
 * @author       Li Zhengxiao
 * @version      1.0
 * @date         2025-05-23
 */
#ifndef T_GPPPRTK_H
#define T_GPPPRTK_H

#include <string>
#include <vector>
#include <map>

#include "gexport/ExportLibGREAT.h"
#include "gaugrdr.h"
#include "gall/gallpar.h"
#include "gdata/gsatdata.h"
#include "gmodels/gpar.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gconst.h"
#include "spdlog/spdlog.h"
#include "newmat/newmat.h"
#include "newmat/newmatap.h"

using namespace std;
using namespace gnut;
// NEWMAT namespace is only available when use_namespace is defined;
// types like Matrix/SymmetricMatrix are in global scope by default.

namespace great {

class LibGREAT_LIBRARY_EXPORT t_gppprtk {
public:
    // Constructor with parsed config values
    t_gppprtk(const string& path, double ionoCfg, double tropCfg,
              t_spdlog spdlog = nullptr);
    ~t_gppprtk();

    void setReader(t_gaugrdr* rdr) { _rdr = rdr; }
    void setSite(const string& site) { _site = site; }

    // Initialize AUG reader with rover coordinates (call after construction)
    bool initReader(const string& site, const t_gtriple& xyz_rov);

    // Main entry: build AUG constraints for this epoch
    // Returns true if constraints are built, false otherwise
    // Output: A_aug, P_aug, l_aug (empty if no constraints)
    bool buildConstraints(const t_gtime& epoch,
                          t_gallpar& param_float,
                          const vector<t_gsatdata>& data,
                          const ColumnVector& dx,
                          const map<string, int>& lock_epo,
                          const SymmetricMatrix& Qx,
                          Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                          int& n_stec, int& n_zwd);

    bool enabled() const { return _enabled; }

    // Query AUG ZWD and its std for initialization (nearest epoch within max_dt).
    // Returns false if no valid AUG data is available or zwd_std <= 0.
    bool queryZwd(const t_gtime& t, double& zwd, double& zwd_std,
                  double max_dt = 30.0) const;

    // Outlier rejection threshold for AUG constraints (m). 0 = disabled.
    void setOutlierThres(double thres) { _outlierThres = thres; }

    // IQR-based outlier rejection for STEC constraints. 0 = disabled, 1 = enabled.
    void setIqrEnabled(bool enabled) { _iqrEnabled = enabled; }

    // Boost STEC weight for first N epochs by factor M (variance divided by M^2).
    // n <= 0 disables boost.
    void setStecBoost(int n, double m) { _stecBoostN = n; _stecBoostM = m; }

private:
    bool _buildStecConstraints(const AugEpoch& aug,
                               t_gallpar& param,
                               const vector<t_gsatdata>& data,
                               const ColumnVector& dx,
                               const map<string, int>& lock_epo,
                               Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                               int& n_stec, int& n_rejected);
    bool _buildZwdConstraints(const AugEpoch& aug,
                              t_gallpar& param,
                              const ColumnVector& dx,
                              const SymmetricMatrix& Qx,
                              Matrix& A_aug, DiagonalMatrix& P_aug, ColumnVector& l_aug,
                              int& n_zwd, int& n_rejected);
    string _selectRefSat(const AugEpoch& aug, const vector<t_gsatdata>& data, char sys,
                         const map<string, int>& lock_epo);
    double _calIonoVar(double el_rad, double var0);
    double _resolveVar(double cfgVal, double augSigma);

    t_gaugrdr* _rdr;
    t_spdlog _spdlog;
    string _site;
    string _path;
    bool _enabled;
    double _ionoCfg;
    double _tropCfg;
    map<char, string> _last_ref;  ///< last epoch reference sat per system (G/E/C/J/I)
    double _outlierThres;         ///< outlier rejection threshold for AUG constraints (m), 0=disabled
    bool _iqrEnabled;             ///< IQR-based outlier rejection for STEC constraints
    int _stecBoostN;              ///< boost STEC weight for first N epochs, 0=disabled
    double _stecBoostM;           ///< noise shrink factor (variance divided by M^2)
    int _epochCount;              ///< number of epochs with AUG constraints processed
};

} // namespace great

#endif
