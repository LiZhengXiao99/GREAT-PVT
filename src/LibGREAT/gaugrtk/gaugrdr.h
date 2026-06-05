/**
 * @file         gaugrdr.h
 * @brief        AUG file reader for PPP-RTK augmentation
 * @author       Li Zhengxiao
 * @version      1.0
 * @date         2025-05-23
 */
#ifndef GAUGRDR_H
#define GAUGRDR_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "gexport/ExportLibGREAT.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gsys.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "spdlog/spdlog.h"

using namespace std;
using namespace gnut;

namespace great {

struct AugSatRec {
    string prn;
    double ipp_lat;
    double ipp_lon;
    double az_deg;
    double el_deg;
    double stec_m;
    double sigma_m;
    int lck;
    double el_max_deg;
    double stec_p1p2;
    int fix_ewl;
    int fix_wl;
    int fix_nl;
};

struct AugEpoch {
    t_gtime epoch;
    double zwd;
    double grd_n;
    double grd_e;
    double zwd_std0;
    double zwd_std1;
    double zwd_std2;
    int sys;
    int nsys;
    vector<AugSatRec> sats;
};

struct AugHeadInfo {
    string filepath;
    string stname;
    t_gtriple xyz;
    double dist_km;
};

class LibGREAT_LIBRARY_EXPORT t_gaugrdr {
public:
    t_gaugrdr(const string& path, const string& site,
              const t_gtriple& xyz_rov, t_spdlog spdlog = nullptr);
    ~t_gaugrdr();
    bool init();
    bool getEpoch(const t_gtime& t, AugEpoch& out);
    bool getNearestEpoch(const t_gtime& t, AugEpoch& out, double max_dt = 30.0) const;
    t_gtriple xyz_svr() const { return _xyz; }
    double distance_km() const { return _dist_km; }
    bool is_open() const { return _is_open; }

private:
    bool _resolvePath();
    bool _parseFile(const string& filepath);
    bool _parseHeader(ifstream& fs);
    bool _parseBody(ifstream& fs);
    bool _filterSat(AugSatRec& rec);

    string _path;
    string _site;
    t_gtriple _xyz_rov;
    t_spdlog _spdlog;
    t_gtriple _xyz;
    double _dist_km;
    bool _is_open;
    bool _is_lsq;
    double _sample;      ///< sampling interval (s), inferred from AUG epochs
    vector<AugEpoch> _epochs;
};

} // namespace great

#endif
