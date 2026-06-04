/*
*
* @file        gtruecrd.cpp
* @brief       true coordinate reader for .true_crd files
* @author      GREAT-WHU
* @version     1.0.0
* @date        2024-08-29
*
* @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
*
*/

#include "gset/gtruecrd.h"

using namespace gnut;

t_gtruecrd::t_gtruecrd()
    : _loaded(false)
{
}

t_gtruecrd::t_gtruecrd(const string &path)
    : _loaded(false)
{
    load(path);
}

t_gtruecrd::~t_gtruecrd()
{
}

bool t_gtruecrd::load(const string &path)
{
    _crd_table.clear();
    _loaded = false;

    ifstream ifs(path);
    if (!ifs.is_open())
    {
        return false;
    }

    string line;
    while (getline(ifs, line))
    {
        // Skip empty lines and comment lines
        if (line.empty())
            continue;

        string trimmed = line;
        // Trim leading whitespace
        size_t start = trimmed.find_first_not_of(" \t\r\n");
        if (start == string::npos)
            continue;
        trimmed = trimmed.substr(start);

        // Skip comment lines starting with '%'
        if (trimmed[0] == '%')
            continue;

        istringstream iss(trimmed);
        string site;
        double x, y, z;

        if (!(iss >> site >> x >> y >> z))
        {
            continue;  // Parse failed, skip this line
        }

        // Normalize site name to uppercase
        transform(site.begin(), site.end(), site.begin(), ::toupper);

        _crd_table[site] = t_gtriple(x, y, z);
    }

    ifs.close();
    _loaded = !_crd_table.empty();
    return _loaded;
}

bool t_gtruecrd::has_site(const string &site) const
{
    string key = site;
    transform(key.begin(), key.end(), key.begin(), ::toupper);
    return _crd_table.find(key) != _crd_table.end();
}

t_gtriple t_gtruecrd::get_crd_xyz(const string &site) const
{
    string key = site;
    transform(key.begin(), key.end(), key.begin(), ::toupper);

    auto it = _crd_table.find(key);
    if (it != _crd_table.end())
    {
        return it->second;
    }

    return t_gtriple(0.0, 0.0, 0.0);
}

set<string> t_gtruecrd::sites() const
{
    set<string> result;
    for (const auto &it : _crd_table)
    {
        result.insert(it.first);
    }
    return result;
}
