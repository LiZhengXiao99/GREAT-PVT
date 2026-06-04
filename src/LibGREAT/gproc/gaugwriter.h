/**
 * @file         gaugwriter.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        AUG file writer for PPP-AR augmentation
 * @version      1.0
 * @date         2024-08-29
 * 
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 * 
 */
#ifndef GAUGWRITER_H
#define GAUGWRITER_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "gexport/ExportLibGREAT.h"
#include "gall/gallpar.h"
#include "gmodels/gpar.h"
#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gutils/gtime.h"
#include "gutils/gsys.h"
#include "gutils/gconst.h"
#include "gambfix/gambiguity.h"
#include "gambfix/gambcommon.h"
#include "gio/giof.h"
#include "gset/gsetout.h"

using namespace std;
using namespace gnut;

namespace great {

class LibGREAT_LIBRARY_EXPORT t_gaugwriter {
public:
    t_gaugwriter(t_gsetbase *gset, t_spdlog spdlog);
    ~t_gaugwriter();

    void setSite(const string& site);
    void setPath(const string& path);
    bool setHeader(const t_gtriple& xyz, const t_gtriple& blh,
                   const string& rcv_type, const string& ant_type, const string& ant_radome,
                   const t_gtime& beg, const t_gtime& end,
                   const set<string>& sys);

    bool writeEpoch(const t_gtime& epoch,
                    const t_gallpar& X_fix,
                    const t_gallpar& X_float,
                    const SymmetricMatrix& Q_fix,
                    const vector<t_gsatdata>& data,
                    const map<string, double>& ele,
                    t_gambiguity* ambfix,
                    const map<string, int>& lock_epo,
                    const map<string, double>& el_max,
                    const string& site);

    void close();

private:
    bool _writeHeader();
    bool _writeObsCodePairs(const vector<t_gsatdata>& data);
    bool _writeEpochBody(const t_gtime& epoch,
                         const t_gallpar& X_fix,
                         const t_gallpar& X_float,
                         const SymmetricMatrix& Q_fix,
                         const vector<t_gsatdata>& data,
                         const map<string, double>& ele,
                         t_gambiguity* ambfix,
                         const map<string, int>& lock_epo,
                         const map<string, double>& el_max,
                         const string& site);

    t_giof* _fout;
    string _site;
    string _path;
    bool _header_written;

    // header cache
    t_gtriple _xyz;
    t_gtriple _blh;
    string _rcv;
    string _ant;
    string _radome;
    t_gtime _beg;
    t_gtime _end;
    string _sys_str;
    t_spdlog _spdlog;
};

} // namespace great

#endif
