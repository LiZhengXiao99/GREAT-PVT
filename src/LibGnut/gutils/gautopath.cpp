/**
 *
 * @verbatim
    History
    2025-05-27  Auto-path file discovery for GREAT-PVT
  @endverbatim
 *
 * Copyright (c) 2025
 *
 * @file        gautopath.cpp
 * @brief       Auto-discover GNSS input files from a base directory
 * @author      Auto-generated
 * @version     1.0.0
 * @date        2025-05-27
 *
 */

#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "gutils/gautopath.h"
#include "gutils/gfileconv.h"

using namespace std;

namespace gnut
{

    // ------------------------------------------------------------------
    // Helper: list directory contents
    // ------------------------------------------------------------------
    vector<string> listDirFiles(const string &path)
    {
        vector<string> result;
        DIR *dir = opendir(path.c_str());
        if (!dir)
            return result;

        struct dirent *entry;
        while ((entry = readdir(dir)) != nullptr)
        {
            string name(entry->d_name);
            if (name == "." || name == "..")
                continue;
            result.push_back(name);
        }
        closedir(dir);
        return result;
    }

    // ------------------------------------------------------------------
    // Helper: simple wildcard match (* = any sequence, ? = single char)
    // ------------------------------------------------------------------
    bool matchWildcard(const string &text, const string &pattern)
    {
        size_t i = 0, j = 0;
        size_t starIdx = string::npos, matchIdx = 0;

        while (i < text.size())
        {
            if (j < pattern.size() && (pattern[j] == '?' || pattern[j] == text[i]))
            {
                ++i;
                ++j;
            }
            else if (j < pattern.size() && pattern[j] == '*')
            {
                starIdx = j;
                matchIdx = i;
                ++j;
            }
            else if (starIdx != string::npos)
            {
                j = starIdx + 1;
                i = ++matchIdx;
            }
            else
            {
                return false;
            }
        }

        while (j < pattern.size() && pattern[j] == '*')
            ++j;

        return j == pattern.size();
    }

    // ------------------------------------------------------------------
    // Helper: GPS week / DOW boundary crossing
    // ------------------------------------------------------------------
    static void _gpsWeekDowBoundary(int &gpsWeek, int &dow, int deltaDays)
    {
        dow += deltaDays;
        while (dow < 0)
        {
            dow += 7;
            gpsWeek -= 1;
        }
        while (dow > 6)
        {
            dow -= 7;
            gpsWeek += 1;
        }
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------
    static string _joinPath(const string &dir, const string &file)
    {
        if (dir.empty())
            return file;
        string d = dir;
        if (d.back() != '/' && d.back() != '\\')
            d += PATH_SEPARATOR;
        return d + file;
    }

    static string _int2doy(int doy)
    {
        ostringstream oss;
        oss << setw(3) << setfill('0') << doy;
        return oss.str();
    }

    static string _int2yy(int year)
    {
        int yy = year % 100;
        ostringstream oss;
        oss << setw(2) << setfill('0') << yy;
        return oss.str();
    }

    static string _int2yyyy(int year)
    {
        ostringstream oss;
        oss << setw(4) << setfill('0') << year;
        return oss.str();
    }

    static string _int2gwk(int week)
    {
        ostringstream oss;
        oss << setw(4) << setfill('0') << week;
        return oss.str();
    }

    static vector<string> _searchPaths(const string &basepath)
    {
        vector<string> paths;
        if (basepath.empty()) return paths;
        paths.push_back(basepath);
        paths.push_back(_joinPath(basepath, "obs"));
        paths.push_back(_joinPath(basepath, "usr_20_obs"));
        paths.push_back(_joinPath(basepath, "gnss"));
        paths.push_back(_joinPath(basepath, "model"));
        paths.push_back(_joinPath(basepath, "upd"));
        paths.push_back(_joinPath(basepath, "poleut1"));
        return paths;
    }

    static vector<string> _findAllMatches(const string &basepath, const vector<string> &patterns)
    {
        vector<string> result;
        if (basepath.empty())
            return result;

        for (const string &sp : _searchPaths(basepath))
        {
            vector<string> files = listDirFiles(sp);
            for (const string &pat : patterns)
            {
                string patLower = pat;
                transform(patLower.begin(), patLower.end(), patLower.begin(), ::tolower);
                for (const string &f : files)
                {
                    string fLower = f;
                    transform(fLower.begin(), fLower.end(), fLower.begin(), ::tolower);
                    if (matchWildcard(fLower, patLower))
                    {
                        string full = _joinPath(sp, f);
                        if (find(result.begin(), result.end(), full) == result.end())
                            result.push_back(full);
                    }
                }
            }
        }
        return result;
    }

    static string _findFirstMatch(const string &basepath, const vector<string> &patterns)
    {
        if (basepath.empty())
            return "";

        for (const string &sp : _searchPaths(basepath))
        {
            vector<string> files = listDirFiles(sp);
            for (const string &pat : patterns)
            {
                string patLower = pat;
                transform(patLower.begin(), patLower.end(), patLower.begin(), ::tolower);
                for (const string &f : files)
                {
                    string fLower = f;
                    transform(fLower.begin(), fLower.end(), fLower.begin(), ::tolower);
                    if (matchWildcard(fLower, patLower))
                        return _joinPath(sp, f);
                }
            }
        }
        return "";
    }

    // Search across multiple days (for SP3/CLK/ERP/IFCB)
    static vector<string> _findAcrossDays(
        const string &basepath,
        int year, int doy, int gpsWeek, int dow,
        const vector<string> &patternsToday,
        const vector<string> &patternsYesterday,
        const vector<string> &patternsTomorrow)
    {
        vector<string> result;

        // Today
        vector<string> today = _findAllMatches(basepath, patternsToday);
        result.insert(result.end(), today.begin(), today.end());

        // Yesterday
        if (!patternsYesterday.empty())
        {
            int w = gpsWeek, d = dow;
            _gpsWeekDowBoundary(w, d, -1);
            int y = year, dy = doy - 1;
            if (dy < 1)
            {
                dy = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
                    dy = 366;
                y = year - 1;
            }
            vector<string> pats;
            for (const string &pat : patternsYesterday)
            {
                string tmp = pat;
                size_t pos = 0;
                while ((pos = tmp.find("$(YESTERDAY_YEAR)", pos)) != string::npos)
                {
                    tmp.replace(pos, 17, _int2yyyy(y));
                    pos += 4;
                }
                pos = 0;
                while ((pos = tmp.find("$(YESTERDAY_DOY)", pos)) != string::npos)
                {
                    tmp.replace(pos, 16, _int2doy(dy));
                    pos += 3;
                }
                pos = 0;
                while ((pos = tmp.find("$(YESTERDAY_GWK)", pos)) != string::npos)
                {
                    tmp.replace(pos, 16, _int2gwk(w));
                    pos += 4;
                }
                pos = 0;
                while ((pos = tmp.find("$(YESTERDAY_DOW)", pos)) != string::npos)
                {
                    ostringstream oss;
                    oss << d;
                    tmp.replace(pos, 16, oss.str());
                    pos += 1;
                }
                pats.push_back(tmp);
            }
            vector<string> yest = _findAllMatches(basepath, pats);
            result.insert(result.end(), yest.begin(), yest.end());
        }

        // Tomorrow
        if (!patternsTomorrow.empty())
        {
            int w = gpsWeek, d = dow;
            _gpsWeekDowBoundary(w, d, +1);
            int y = year, dy = doy + 1;
            int maxDoy = 365;
            if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
                maxDoy = 366;
            if (dy > maxDoy)
            {
                dy = 1;
                y = year + 1;
            }
            vector<string> pats;
            for (const string &pat : patternsTomorrow)
            {
                string tmp = pat;
                size_t pos = 0;
                while ((pos = tmp.find("$(TOMORROW_YEAR)", pos)) != string::npos)
                {
                    tmp.replace(pos, 16, _int2yyyy(y));
                    pos += 4;
                }
                pos = 0;
                while ((pos = tmp.find("$(TOMORROW_DOY)", pos)) != string::npos)
                {
                    tmp.replace(pos, 15, _int2doy(dy));
                    pos += 3;
                }
                pos = 0;
                while ((pos = tmp.find("$(TOMORROW_GWK)", pos)) != string::npos)
                {
                    tmp.replace(pos, 15, _int2gwk(w));
                    pos += 4;
                }
                pos = 0;
                while ((pos = tmp.find("$(TOMORROW_DOW)", pos)) != string::npos)
                {
                    ostringstream oss;
                    oss << d;
                    tmp.replace(pos, 15, oss.str());
                    pos += 1;
                }
                pats.push_back(tmp);
            }
            vector<string> tom = _findAllMatches(basepath, pats);
            result.insert(result.end(), tom.begin(), tom.end());
        }

        return result;
    }

    // ------------------------------------------------------------------
    // Find RINEX observation files
    // ------------------------------------------------------------------
    vector<string> findObsFiles(const string &basepath, const string &site, int year, int doy)
    {
        if (basepath.empty())
            return {};

        string s = site;
        transform(s.begin(), s.end(), s.begin(), ::toupper);

        vector<string> result;

        // Broad search when year/doy is not specified (e.g., beg=0)
        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                s + "*.??o", s + "*.??O",
                "*" + s + "*.??o", "*" + s + "*.??O",
                "*" + s + "*_MO*.rnx", "*" + s + "*_GO*.rnx", "*" + s + "*_RO*.rnx",
                "*" + s + "*_MO*.RNX", "*" + s + "*_GO*.RNX", "*" + s + "*_RO*.RNX",
                "*" + s + "*.rnx", "*" + s + "*.RNX",
            };
            result = _findAllMatches(basepath, patterns);
        }
        else
        {
            string yy = _int2yy(year);
            string ddd = _int2doy(doy);

            vector<string> patterns = {
                s + "*" + ddd + "*" + yy + "o",
                s + "*" + ddd + "*" + yy + "O",
                "*" + s + "*" + ddd + "*" + yy + "o",
                "*" + s + "*" + ddd + "*" + yy + "O",
                "*" + s + "*_MO*.rnx",
                "*" + s + "*_GO*.rnx",
                "*" + s + "*_RO*.rnx",
                "*" + s + "*_MO*.RNX",
                "*" + s + "*_GO*.RNX",
                "*" + s + "*_RO*.RNX",
                "*" + s + "*.rnx",
                "*" + s + "*.RNX",
            };

            result = _findAllMatches(basepath, patterns);
        }

        // Filter out navigation files to avoid conflict with findEphFiles
        vector<string> filtered;
        for (const string &f : result)
        {
            string fname = f;
            transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
            if (fname.find("_mn.") != string::npos) continue; // Mixed Navigation
            if (fname.find("_gn.") != string::npos) continue; // GLONASS Navigation
            filtered.push_back(f);
        }
        return filtered;
    }

    // ------------------------------------------------------------------
    // Find RINEX navigation files
    // ------------------------------------------------------------------
    vector<string> findEphFiles(const string &basepath, int year, int doy)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                "*.??p", "*.??P",
                "*_MN.rnx", "*_MN.RNX",
            };
            return _findAllMatches(basepath, patterns);
        }

        string yy = _int2yy(year);
        string ddd = _int2doy(doy);
        string yyyy = _int2yyyy(year);

        vector<string> patterns = {
            "brdm" + ddd + "0." + yy + "p",
            "brdm" + ddd + "0." + yy + "P",
            "brdc" + ddd + "0." + yy + "p",
            "brdc" + ddd + "0." + yy + "P",
            "BRDM" + ddd + "0." + yy + "p",
            "BRDC" + ddd + "0." + yy + "p",
            "*" + yyyy + ddd + "0000_01D_MN.rnx",
            "*" + yyyy + ddd + "0000_01D_MN.RNX",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find SP3 orbit files
    // ------------------------------------------------------------------
    vector<string> findSp3Files(const string &basepath, int year, int doy, int gpsWeek, int dow)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                "*.sp3", "*.SP3",
            };
            return _findAllMatches(basepath, patterns);
        }

        string gwk = _int2gwk(gpsWeek);
        string yyyy = _int2yyyy(year);
        string ddd = _int2doy(doy);
        ostringstream oss;
        oss << dow;
        string sdow = oss.str();

        vector<string> today = {
            "*" + gwk + sdow + ".sp3",
            "*" + gwk + sdow + ".SP3",
            "*" + yyyy + ddd + "0000_*_ORB.SP3",
            "*" + yyyy + ddd + "0000_*_ORB.sp3",
        };

        vector<string> yesterday = {
            "*$(YESTERDAY_GWK)$(YESTERDAY_DOW).sp3",
            "*$(YESTERDAY_GWK)$(YESTERDAY_DOW).SP3",
            "*$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_*_ORB.SP3",
            "*$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_*_ORB.sp3",
        };

        vector<string> tomorrow = {
            "*$(TOMORROW_GWK)$(TOMORROW_DOW).sp3",
            "*$(TOMORROW_GWK)$(TOMORROW_DOW).SP3",
            "*$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_*_ORB.SP3",
            "*$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_*_ORB.sp3",
        };

        return _findAcrossDays(basepath, year, doy, gpsWeek, dow, today, yesterday, tomorrow);
    }

    // ------------------------------------------------------------------
    // Find RINEX clock files
    // ------------------------------------------------------------------
    vector<string> findClkFiles(const string &basepath, int year, int doy, int gpsWeek, int dow)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                "*.clk", "*.CLK",
            };
            return _findAllMatches(basepath, patterns);
        }

        string gwk = _int2gwk(gpsWeek);
        string yyyy = _int2yyyy(year);
        string ddd = _int2doy(doy);
        ostringstream oss;
        oss << dow;
        string sdow = oss.str();

        vector<string> today = {
            "*" + gwk + sdow + ".clk",
            "*" + gwk + sdow + ".CLK",
            "*" + yyyy + ddd + "0000_*_CLK.CLK",
            "*" + yyyy + ddd + "0000_*_CLK.clk",
        };

        vector<string> yesterday = {
            "*$(YESTERDAY_GWK)$(YESTERDAY_DOW).clk",
            "*$(YESTERDAY_GWK)$(YESTERDAY_DOW).CLK",
            "*$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_*_CLK.CLK",
            "*$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_*_CLK.clk",
        };

        vector<string> tomorrow = {
            "*$(TOMORROW_GWK)$(TOMORROW_DOW).clk",
            "*$(TOMORROW_GWK)$(TOMORROW_DOW).CLK",
            "*$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_*_CLK.CLK",
            "*$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_*_CLK.clk",
        };

        return _findAcrossDays(basepath, year, doy, gpsWeek, dow, today, yesterday, tomorrow);
    }

    // ------------------------------------------------------------------
    // Find SINEX files
    // ------------------------------------------------------------------
    vector<string> findSnxFiles(const string &basepath, int gpsWeek)
    {
        if (basepath.empty())
            return {};

        string gwk = _int2gwk(gpsWeek);

        vector<string> patterns = {
            "*" + gwk + ".snx",
            "*" + gwk + ".SNX",
            "igs*" + gwk + ".snx",
            "igs*" + gwk + ".SNX",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find ERP/EOP files
    // ------------------------------------------------------------------
    vector<string> findErpFiles(const string &basepath, int year, int doy)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                "poleut1", "poleut1_new",
                "*.erp", "*.ERP",
            };
            return _findAllMatches(basepath, patterns);
        }

        string yyyy = _int2yyyy(year);
        string ddd = _int2doy(doy);

        vector<string> today = {
            "poleut1",
            "*0MGX_" + yyyy + ddd + "0000_01D_01D_ERP.ERP",
            "*0MGX_" + yyyy + ddd + "0000_01D_01D_ERP.erp",
            "*.erp",
            "*.ERP",
        };

        vector<string> yesterday = {
            "*0MGX_$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_01D_01D_ERP.ERP",
            "*0MGX_$(YESTERDAY_YEAR)$(YESTERDAY_DOY)0000_01D_01D_ERP.erp",
        };

        vector<string> tomorrow = {
            "*0MGX_$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_01D_01D_ERP.ERP",
            "*0MGX_$(TOMORROW_YEAR)$(TOMORROW_DOY)0000_01D_01D_ERP.erp",
        };

        // For ERP, try today first, then yesterday/tomorrow
        vector<string> result = _findAllMatches(basepath, today);
        if (!result.empty())
            return result;

        return _findAcrossDays(basepath, year, doy, 0, 0, {}, yesterday, tomorrow);
    }

    // ------------------------------------------------------------------
    // Find BIAS/OSB files
    // ------------------------------------------------------------------
    vector<string> findBiasFiles(const string &basepath, int year, int doy)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            vector<string> patterns = {
                "*.BSX", "*.bsx", "*.BIA", "*.bia",
            };
            return _findAllMatches(basepath, patterns);
        }

        string yyyy = _int2yyyy(year);
        string ddd = _int2doy(doy);

        vector<string> patterns = {
            "CAS0MGXRAP_" + yyyy + ddd + "0000_01D_01D_DCB.BSX",
            "CAS0MGXRAP_" + yyyy + ddd + "0000_01D_01D_DCB.bsx",
            "*0MGX_" + yyyy + ddd + "0000_01D_01D_OSB.BIA",
            "*0MGX_" + yyyy + ddd + "0000_01D_01D_OSB.bia",
            "WUM0MGXRAP_" + yyyy + ddd + "0000_01D_01D_ABS.BIA",
            "WUM0MGXRAP_" + yyyy + ddd + "0000_01D_01D_ABS.bia",
            "*.BSX",
            "*.bsx",
            "*.BIA",
            "*.bia",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find SGG IFCB files
    // ------------------------------------------------------------------
    vector<string> findSGGIfcbs(const string &basepath, int year, int doy)
    {
        if (basepath.empty())
            return {};

        ostringstream oss;
        oss << year << doy;
        string prefix = "ifcb_" + oss.str();

        vector<string> today = {
            prefix + "*",
        };

        // Yesterday
        int y1 = year, d1 = doy - 1;
        if (d1 < 1)
        {
            d1 = 365;
            if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
                d1 = 366;
            y1 = year - 1;
        }
        ostringstream oss1;
        oss1 << y1 << d1;
        vector<string> yesterday = {
            "ifcb_" + oss1.str() + "*",
        };

        // Tomorrow
        int y2 = year, d2 = doy + 1;
        int maxDoy = 365;
        if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
            maxDoy = 366;
        if (d2 > maxDoy)
        {
            d2 = 1;
            y2 = year + 1;
        }
        ostringstream oss2;
        oss2 << y2 << d2;
        vector<string> tomorrow = {
            "ifcb_" + oss2.str() + "*",
        };

        return _findAcrossDays(basepath, year, doy, 0, 0, today, yesterday, tomorrow);
    }

    // ------------------------------------------------------------------
    // Find UPD files
    // ------------------------------------------------------------------
    vector<string> findUPDFiles(const string &basepath, int year, int doy, const string &type)
    {
        if (basepath.empty())
            return {};

        if (year == 0 || doy == 0)
        {
            string t = type;
            if (!t.empty()) t = t + "*";
            vector<string> patterns = {
                "upd_" + t + "*",
            };
            return _findAllMatches(basepath, patterns);
        }

        ostringstream oss;
        oss << "upd_" << type << "_" << year << doy << "*";

        vector<string> patterns = {
            oss.str(),
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find ATX files
    // ------------------------------------------------------------------
    vector<string> findATXFiles(const string &basepath)
    {
        if (basepath.empty())
            return {};

        vector<string> patterns = {
            "igs20_*.atx",
            "igs20_*.ATX",
            "igs14_*.atx",
            "igs14_*.ATX",
            "*.atx",
            "*.ATX",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find BLQ files
    // ------------------------------------------------------------------
    vector<string> findBLQFiles(const string &basepath)
    {
        if (basepath.empty())
            return {};

        vector<string> patterns = {
            "oceanload",
            "OCEANLOAD",
            "*.blq",
            "*.BLQ",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find DE405 files
    // ------------------------------------------------------------------
    vector<string> findDEFiles(const string &basepath)
    {
        if (basepath.empty())
            return {};

        vector<string> patterns = {
            "jpleph_de405*",
            "DE405*",
            "de405*",
        };

        return _findAllMatches(basepath, patterns);
    }

    // ------------------------------------------------------------------
    // Find leap second files
    // ------------------------------------------------------------------
    vector<string> findLeapFiles(const string &basepath)
    {
        if (basepath.empty())
            return {};

        vector<string> patterns = {
            "Leap_Second*",
            "leapsec*",
            "*.dat",
            "*.DAT",
        };

        return _findAllMatches(basepath, patterns);
    }

} // namespace
