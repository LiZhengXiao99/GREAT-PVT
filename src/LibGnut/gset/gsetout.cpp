/**
*
* @verbatim
     (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

     (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
        Research Institute of Geodesy, Topography and Cartography
        Ondrejov 244, 251 65, Czech Republic
  @endverbatim
*
* @file        gsetout.cpp
* @brief       implements output setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetout.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut
{
    OFMT t_gsetout::str2ofmt(const string &s)
    {
        string tmp = s;
        transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (tmp == "OUT")
            return XXX_OUT;
        if (tmp == "LOG")
            return LOG_OUT;
        if (tmp == "PPP")
            return PPP_OUT;
        if (tmp == "FLT")
            return FLT_OUT;
        if (tmp == "FLT_FLOAT")
            return FLT_FLOAT_OUT;
        if (tmp == "AUG")
            return AUG_OUT;
        if (tmp == "FLT_PPPRTK")
            return FLT_PPPRTK_OUT;
        return OFMT(-1);
    }

    string t_gsetout::ofmt2str(const OFMT &f)
    {
        switch (f)
        {
        case XXX_OUT:
            return "OUT";
        case LOG_OUT:
            return "LOG";
        case PPP_OUT:
            return "PPP";
        case FLT_OUT:
            return "FLT";
        case FLT_FLOAT_OUT:
            return "FLT_FLOAT";
        case AUG_OUT:
            return "AUG";
        case FLT_PPPRTK_OUT:
            return "FLT_PPPRTK";
        default:
            return "UNDEF";
        }
        return "UNDEF";
    }

    t_gsetout::t_gsetout()
        : t_gsetbase(),
          _append(false),
          _verb(0),
          _ctx_year(0),
          _ctx_doy(0),
          _ctx_set(false)
    {
        _set.insert(XMLKEY_OUT);
    }

    t_gsetout::~t_gsetout()
    {
    }

    int t_gsetout::output_size(const string &fmt)
    {
        _gmutex.lock();

        int tmp = _outputs(fmt).size();

        _gmutex.unlock();
        return tmp;
    }

    int t_gsetout::verb()
    {
        _gmutex.lock();

        int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("verb").as_int();

        _gmutex.unlock();
        return tmp;
    }

    bool t_gsetout::append()
    {
        _gmutex.lock();

        bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("append").as_bool();

        _gmutex.unlock();
        return tmp;
    }

    string t_gsetout::outputs(const string &fmt)
    {
        _gmutex.lock();

        string tmp = _outputs(fmt);
        if (tmp == "AUTO")
        {
            if (_ctx_set)
            {
                string mode_str = _ctx_mode.empty() ? ofmt2str(str2ofmt(fmt)) : _ctx_mode;
                tmp = _default_output_path(fmt, _ctx_site, _ctx_year, _ctx_doy,
                                           mode_str, _ctx_iono, _ctx_sys);
            }
            else
            {
                tmp = "";
            }
        }

        _gmutex.unlock();
        return tmp;
    }

    string t_gsetout::outputs(const string &fmt, const string &site,
                              int year, int doy,
                              const string &mode_str,
                              const string &iono_str,
                              const string &sys_str)
    {
        _gmutex.lock();

        string tmp = _outputs(fmt);
        if (tmp == "AUTO")
        {
            tmp = _default_output_path(fmt, site, year, doy, mode_str, iono_str, sys_str);
        }
        else if (!tmp.empty())
        {
            tmp = substitute_placeholders(tmp, site, year, doy, mode_str, iono_str, sys_str);
        }

        _gmutex.unlock();
        return tmp;
    }

    string t_gsetout::substitute_placeholders(const string &path,
                                              const string &site,
                                              int year, int doy,
                                              const string &mode_str,
                                              const string &iono_str,
                                              const string &sys_str)
    {
        string result = path;
        if (result.empty()) return result;

        // strip file:// prefix if present, substitute, then re-add
        bool hasPrefix = false;
        if (result.find(GFILE_PREFIX) == 0)
        {
            hasPrefix = true;
            result = result.substr(strlen(GFILE_PREFIX));
        }

        // Helper to replace all occurrences
        auto replaceAll = [](string &s, const string &from, const string &to)
        {
            size_t pos = 0;
            while ((pos = s.find(from, pos)) != string::npos)
            {
                s.replace(pos, from.length(), to);
                pos += to.length();
            }
        };

        replaceAll(result, "$(MODE)",    mode_str);
        replaceAll(result, "$(mode)",    mode_str);
        replaceAll(result, "$(IONO)",    iono_str);
        replaceAll(result, "$(iono)",    iono_str);
        replaceAll(result, "$(SYSTEM)",  sys_str);
        replaceAll(result, "$(system)",  sys_str);
        replaceAll(result, "$(REC)",     site);
        replaceAll(result, "$(rec)",     site);
        replaceAll(result, "$(STATION)", site);
        replaceAll(result, "$(station)", site);
        replaceAll(result, "$(YEAR)",    int2str(year, 4));
        replaceAll(result, "$(year)",    int2str(year, 4));
        replaceAll(result, "$(DOY)",     int2str(doy, 3));
        replaceAll(result, "$(doy)",     int2str(doy, 3));

        // Extension placeholder: derive from path or from known format
        string ext;
        size_t dot = result.find_last_of('.');
        if (dot != string::npos) ext = result.substr(dot + 1);
        replaceAll(result, "$(EXT)", ext);

        if (hasPrefix)
            result = string(GFILE_PREFIX) + result;

        return result;
    }

    string t_gsetout::log_type()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child("log").attribute("type").as_string();
        str_erase(tmp);
        transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        _gmutex.unlock();
        if (tmp.empty())
        {
            return "CONSOLE";
        }
        else
        {
            return tmp;
        }
    }

    string t_gsetout::log_name()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child("log").attribute("name").as_string();
        str_erase(tmp);
        _gmutex.unlock();
        if (tmp.empty())
        {
            return "my_logger";
        }
        else
        {
            return tmp;
        }
    }

    string t_gsetout::log_pattern()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child("log").attribute("pattern").as_string();
        str_erase(tmp);
        transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        _gmutex.unlock();
        if (tmp.empty())
        {
            return string("[%Y-%m-%d %H:%M:%S] <thread %t> [%l] [%@] %v");
        }
        else
        {
            return tmp;
        }
    }

    level::level_enum t_gsetout::log_level()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child("log").attribute("level").as_string();
        str_erase(tmp);
        transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        _gmutex.unlock();
        if (tmp.empty())
        {
            return level::level_enum::info;
        }
        else
        {
            if (tmp.find("ERROR") != string::npos)
            {
                return level::level_enum::err;
            }
            else if (tmp.find("DEBUG") != string::npos)
            {
                return level::level_enum::debug;
            }
            else if (tmp.find("WARN") != string::npos)
            {
                return level::level_enum::warn;
            }
            else if (tmp.find("CRITICAL") != string::npos)
            {
                return level::level_enum::critical;
            }
            else if (tmp.find("TRACE") != string::npos)
            {
                return level::level_enum::trace;
            }
            else if (tmp.find("INFO") != string::npos)
            {
                return level::level_enum::info;
            }
            else
            {
                return level::level_enum::off;
            }
        }
    }

    string t_gsetout::version(const string &fmt)
    {
        _gmutex.lock();

        string ver = DEFAULT_FILE_VER;
        xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
        if (!fmt.empty() &&
            !node.attribute("ver").empty())
        {
            ver = node.attribute("ver").as_string();
        }

        _gmutex.unlock();
        return ver;
    }

    set<string> t_gsetout::oformats()
    {
        return _oformats();
    }

    set<string> t_gsetout::_oformats()
    {
        set<string> tmp;
        for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling())
        {
            tmp.insert(node.name());
        }
        return tmp;
    }

    void t_gsetout::set_context(const string &site, int year, int doy,
                                const string &mode, const string &iono, const string &sys)
    {
        _gmutex.lock();
        _ctx_site = site;
        _ctx_year = year;
        _ctx_doy = doy;
        _ctx_mode = mode;
        _ctx_iono = iono;
        _ctx_sys = sys;
        _ctx_set = true;
        _gmutex.unlock();
    }

    bool t_gsetout::ctx_set() const
    {
        return _ctx_set;
    }

    int t_gsetout::ctx_year() const
    {
        return _ctx_year;
    }

    int t_gsetout::ctx_doy() const
    {
        return _ctx_doy;
    }

    string t_gsetout::_outputs(const string &fmt)
    {
        string str;
        for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling())
        {
            if (node.name() == fmt)
            {
                string val = node.child_value();
                // Trim leading/trailing whitespace only
                size_t start = val.find_first_not_of(" \t\n\r");
                if (start == string::npos) return "";
                size_t end = val.find_last_not_of(" \t\n\r");
                string trimmed = val.substr(start, end - start + 1);

                istringstream is(trimmed);
                string first;
                is >> first;
                if (first == "0")
                {
                    return "AUTO";
                }
                else if (first == "1")
                {
                    if (is >> str && !is.fail())
                    {
                        if (str.find("://") == string::npos)
                            str = GFILE_PREFIX + str;
                        if (_ctx_set)
                        {
                            string mode_str = ofmt2str(str2ofmt(fmt));
                            str = substitute_placeholders(str, _ctx_site, _ctx_year, _ctx_doy,
                                                          mode_str, _ctx_iono, _ctx_sys);
                        }
                        return str;
                    }
                    return "";
                }
                else
                {
                    // Old format: remove all whitespace
                    str_erase(val);
                    if (val.empty()) return "";
                    istringstream is_old(val);
                    while (is_old >> str && !is_old.fail())
                    {
                        if (str.find("://") == string::npos)
                            str = GFILE_PREFIX + str;
                        if (_ctx_set)
                        {
                            string mode_str = ofmt2str(str2ofmt(fmt));
                            str = substitute_placeholders(str, _ctx_site, _ctx_year, _ctx_doy,
                                                          mode_str, _ctx_iono, _ctx_sys);
                        }
                        return str;
                    }
                }
            }
        }
        return "";
    }

    string t_gsetout::_default_output_path(const string &fmt, const string &site,
                                           int year, int doy,
                                           const string &mode_str,
                                           const string &iono_str,
                                           const string &sys_str)
    {
        string basepath;
        {
            string tmp = _doc.child(XMLKEY_ROOT).child("inputs").child_value("basepath");
            str_erase(tmp);
            basepath = tmp;
        }
        if (basepath.empty()) basepath = ".";

        string ext;
        if (fmt == "flt" || fmt == "flt_float" || fmt == "flt_ppprtk") ext = "flt";
        else if (fmt == "kml") ext = "kml";
        else if (fmt == "aug") ext = "aug";
        else if (fmt == "ppp") ext = "log";
        else ext = "txt";

        string suffix;
        if (fmt == "flt") suffix = "_pppar";
        else if (fmt == "flt_float") suffix = "_float";
        else if (fmt == "flt_ppprtk") suffix = "_ppprtk";

        ostringstream oss;
        oss << basepath << "/result_" << mode_str << "_" << iono_str << "_" << sys_str << "/";
        if (year > 0 && doy > 0)
            oss << site << setfill('0') << setw(4) << year << setw(3) << doy;
        else
            oss << site << "0000000";
        oss << suffix << "." << ext;

        string result = oss.str();
        if (result.find("://") == string::npos)
            result = GFILE_PREFIX + result;
        return result;
    }

    void t_gsetout::check()
    {
        _gmutex.lock();

        // check existence of nodes/attributes
        xml_node parent = _doc.child(XMLKEY_ROOT);
        xml_node node = _default_node(parent, XMLKEY_OUT);

        // check existence of attributes
        _default_attr(node, "append", _append);

        // check supported input formats (see OFMT enum !)
        set<string> ofmt = _oformats();
        set<string>::const_iterator itFMT = ofmt.begin();
        while (itFMT != ofmt.end())
        {
            string fmt = *itFMT;
            OFMT ofmt = str2ofmt(fmt);
            if (ofmt < 0)
            {
                _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child(node.child(fmt.c_str()));
                cout << "Warning: " + fmt + " out format not implemented [gsetout::check()]!" << endl;
                itFMT++;
                continue;
            }

            // check application-specific output format
            if (_OFMT_supported.find(ofmt) == _OFMT_supported.end())
            {
                _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child(node.child(fmt.c_str()));

                cout << "Warning: " + fmt + " out format not supported by this application!" << endl;
            }
            itFMT++;
        }

        _gmutex.unlock();
        return;
    }

    void t_gsetout::help()
    {
        _gmutex.lock();

        cerr << " <outputs append=\"" << _append << "\" verb=\"" << _verb << "\" >\n"
             << "   <flt> file://dir/name </flt>    \t\t <!-- filter output encoder -->\n"
             << " </outputs>\n";

        _gmutex.unlock();
        return;
    }

} // namespace
