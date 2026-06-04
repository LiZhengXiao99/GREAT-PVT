/**
*
* @verbatim
    History
    2025-05-27: created for auto file discovery
  @endverbatim
*
* @file        gautopath.h
* @brief       auto file discovery based on basepath + time info
* @author      AI Assistant
* @version     1.0.0
* @date        2025-05-27
*
*/

#ifndef GAUTOPATH_H
#define GAUTOPATH_H

#include <string>
#include <vector>
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut
{
    // ------------------------------------------------------------------
    // Bulk finders (return all matches) -- declared FIRST so inline
    // wrappers below can see them.
    // ------------------------------------------------------------------
    LibGnut_LIBRARY_EXPORT vector<string> findObsFiles(const string& basepath, const string& site, int year, int doy);
    LibGnut_LIBRARY_EXPORT vector<string> findEphFiles(const string& basepath, int year, int doy);
    LibGnut_LIBRARY_EXPORT vector<string> findSp3Files(const string& basepath, int year, int doy, int gpsWeek, int dow);
    LibGnut_LIBRARY_EXPORT vector<string> findClkFiles(const string& basepath, int year, int doy, int gpsWeek, int dow);
    LibGnut_LIBRARY_EXPORT vector<string> findSnxFiles(const string& basepath, int gpsWeek);
    LibGnut_LIBRARY_EXPORT vector<string> findErpFiles(const string& basepath, int year, int doy);
    LibGnut_LIBRARY_EXPORT vector<string> findBiasFiles(const string& basepath, int year, int doy);
    LibGnut_LIBRARY_EXPORT vector<string> findSGGIfcbs(const string& basepath, int year, int doy);
    LibGnut_LIBRARY_EXPORT vector<string> findUPDFiles(const string& basepath, int year, int doy, const string& type = "");
    LibGnut_LIBRARY_EXPORT vector<string> findATXFiles(const string& basepath);
    LibGnut_LIBRARY_EXPORT vector<string> findBLQFiles(const string& basepath);
    LibGnut_LIBRARY_EXPORT vector<string> findDEFiles(const string& basepath);
    LibGnut_LIBRARY_EXPORT vector<string> findLeapFiles(const string& basepath);
    LibGnut_LIBRARY_EXPORT vector<string> findTrueCrdFiles(const string& basepath);

    // ------------------------------------------------------------------
    // Convenience wrappers returning the first match (empty if none)
    // These mirror the interface used by gsetinp.cpp
    // ------------------------------------------------------------------
    inline string findObsFile(const string& basepath, const string& site, int year, int doy)
    {
        vector<string> v = findObsFiles(basepath, site, year, doy);
        return v.empty() ? "" : v[0];
    }

    inline string findEphFile(const string& basepath, int year, int doy)
    {
        vector<string> v = findEphFiles(basepath, year, doy);
        return v.empty() ? "" : v[0];
    }

    inline string findSp3File(const string& basepath, int year, int doy, int gpsWeek, int dow)
    {
        vector<string> v = findSp3Files(basepath, year, doy, gpsWeek, dow);
        return v.empty() ? "" : v[0];
    }

    inline string findClkFile(const string& basepath, int year, int doy, int gpsWeek, int dow)
    {
        vector<string> v = findClkFiles(basepath, year, doy, gpsWeek, dow);
        return v.empty() ? "" : v[0];
    }

    inline string findSnxFile(const string& basepath, int gpsWeek)
    {
        vector<string> v = findSnxFiles(basepath, gpsWeek);
        return v.empty() ? "" : v[0];
    }

    inline string findErpFile(const string& basepath, int year, int doy)
    {
        vector<string> v = findErpFiles(basepath, year, doy);
        return v.empty() ? "" : v[0];
    }

    inline string findBiasFile(const string& basepath, int year, int doy)
    {
        vector<string> v = findBiasFiles(basepath, year, doy);
        return v.empty() ? "" : v[0];
    }

    inline string findDCBFile(const string& basepath, int year, int doy)
    {
        return findBiasFile(basepath, year, doy);
    }

    inline string findSGGIfcb(const string& basepath, int year, int doy)
    {
        vector<string> v = findSGGIfcbs(basepath, year, doy);
        return v.empty() ? "" : v[0];
    }

    inline string findUPDFile(const string& basepath, int year, int doy, const string& type)
    {
        vector<string> v = findUPDFiles(basepath, year, doy, type);
        return v.empty() ? "" : v[0];
    }

    inline string findATXFile(const string& basepath)
    {
        vector<string> v = findATXFiles(basepath);
        return v.empty() ? "" : v[0];
    }

    inline string findBLQFile(const string& basepath)
    {
        vector<string> v = findBLQFiles(basepath);
        return v.empty() ? "" : v[0];
    }

    inline string findDEFile(const string& basepath)
    {
        vector<string> v = findDEFiles(basepath);
        return v.empty() ? "" : v[0];
    }

    inline string findLeapFile(const string& basepath)
    {
        vector<string> v = findLeapFiles(basepath);
        return v.empty() ? "" : v[0];
    }

    inline string findTrueCrdFile(const string& basepath)
    {
        vector<string> v = findTrueCrdFiles(basepath);
        return v.empty() ? "" : v[0];
    }

    // ------------------------------------------------------------------
    // Low-level helpers
    // ------------------------------------------------------------------
    LibGnut_LIBRARY_EXPORT vector<string> listDirFiles(const string& path);
    LibGnut_LIBRARY_EXPORT bool matchWildcard(const string& text, const string& pattern);

} // namespace

#endif
