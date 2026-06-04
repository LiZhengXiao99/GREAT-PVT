/*
*
* @file        gtruecrd.h
* @brief       true coordinate reader for .true_crd files
* @author      GREAT-WHU
* @version     1.0.0
* @date        2024-08-29
*
* @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
*
*/

#ifndef GTRUECRD_H
#define GTRUECRD_H

#include <map>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "gutils/gtriple.h"
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut
{

    class LibGnut_LIBRARY_EXPORT t_gtruecrd
    {
    public:
        /** @brief default constructor */
        t_gtruecrd();

        /** @brief constructor with file path */
        explicit t_gtruecrd(const string &path);

        /** @brief destructor */
        ~t_gtruecrd();

        /** @brief load true_crd file */
        bool load(const string &path);

        /** @brief check if file loaded successfully */
        bool is_loaded() const { return !_crd_table.empty(); }

        /** @brief check if site exists */
        bool has_site(const string &site) const;

        /** @brief get XYZ coordinate for a site (returns zero triple if not found) */
        t_gtriple get_crd_xyz(const string &site) const;

        /** @brief get all available sites */
        set<string> sites() const;

    private:
        map<string, t_gtriple> _crd_table;  ///< site -> XYZ
        bool _loaded;
    };

} // namespace gnut

#endif // GTRUECRD_H
