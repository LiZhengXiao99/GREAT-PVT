/**
*
* @verbatim
(c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
Ondrejov 244, 251 65, Czech Republic
@endverbatim
*
* @file        gsetinp.cpp
* @brief       implements input setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"
#include "gutils/gmutex.h"
#include "gutils/gautopath.h"
using namespace std;
using namespace pugi;

namespace gnut
{
    IFMT t_gsetinp::str2ifmt(const string &s)
    {
        string tmp = s;
        transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (tmp == "RINEXC" || tmp == "RNC")
            return IFMT::RINEXC_INP;
        if (tmp == "RINEXO" || tmp == "RNO")
            return IFMT::RINEXO_INP;
        if (tmp == "RINEXN" || tmp == "RNN")
            return IFMT::RINEXN_INP;
        if (tmp == "SP3")
            return IFMT::SP3_INP;
        if (tmp == "ATX")
            return IFMT::ATX_INP;
        if (tmp == "BLQ")
            return IFMT::BLQ_INP;
        if (tmp == "BIASINEX")
            return IFMT::BIASINEX_INP;
        if (tmp == "BIAS" || tmp == "BIABERN")
            return IFMT::BIAS_INP;
        if (tmp == "DE")
            return IFMT::DE_INP;
        if (tmp == "EOP" || tmp == "POLEUT1") 
            return IFMT::EOP_INP; // optinal for xml node
        if (tmp == "LEAPSECOND")
            return IFMT::LEAPSECOND_INP;
        if (tmp == "UPD")
            return IFMT::UPD_INP;
        if (tmp == "IFCB")
            return IFMT::IFCB_INP;

        string message = "The Type : " + tmp + " is not support, check your xml";
        spdlog::warn(message);
        throw logic_error(message);
    }

    string t_gsetinp::ifmt2str(const IFMT &f)
    {
        switch (f)
        {
        case IFMT::RINEXO_INP:
            return "RINEXO";
        case IFMT::RINEXC_INP:
            return "RINEXC";
        case IFMT::RINEXN_INP:
            return "RINEXN";
        case IFMT::SP3_INP:
            return "SP3";
        case IFMT::ATX_INP:
            return "ATX";
        case IFMT::BLQ_INP:
            return "BLQ";
        case IFMT::BIASINEX_INP:
            return "BIASINEX";
        case IFMT::BIAS_INP:
            return "BIAS";
        case IFMT::UPD_INP:
            return "UPD";
        case IFMT::IFCB_INP:
            return "IFCB";
        case IFMT::LEAPSECOND_INP:
            return "LEAPSECOND";
        case IFMT::DE_INP:
            return "DE";
		case IFMT::EOP_INP:
			return "EOP";
        default:
            spdlog::critical("No fmt for {}, check your inp.", f);
            throw logic_error("check your inp");
        }
    }

    t_gsetinp::t_gsetinp()
        : t_gsetbase()
    {
        _set.insert(XMLKEY_INP);
        _chkNavig = true;
        _chkHealth = true;
        _corrStream = "";
    }

    t_gsetinp::~t_gsetinp()
    {
    }

    int t_gsetinp::input_size(const string &fmt)
    {
        _gmutex.lock();

        int tmp = _inputs(fmt).size();

        _gmutex.unlock();
        return tmp;
    }

    bool t_gsetinp::check_input(const string &fmt)
    {
        _gmutex.lock();
        int tmp = _inputs(fmt).size();
        _gmutex.unlock();
        return (tmp > 0);
    }

    void t_gsetinp::check_input(const string &fmt, const string &message)
    {
        _gmutex.lock();
        int tmp = _inputs(fmt).size();
        _gmutex.unlock();
        if (tmp < 0)
        {
            spdlog::critical(message);
            throw logic_error(message);
        }
    }

    multimap<IFMT, string> t_gsetinp::inputs_all()
    {
        _gmutex.lock();

        multimap<IFMT, string> map;

        set<string> ifmt = _iformats();
        set<string>::const_iterator itFMT = ifmt.begin();

        while (itFMT != ifmt.end())
        {
            string fmt = *itFMT;
            if (fmt.empty() || fmt == "basepath" || fmt == "tbl")
            {
                itFMT++;
                continue;
            }

            IFMT ifmt = str2ifmt(fmt);
            vector<string> inputs = _inputs(fmt); //get file name in input node
            vector<string>::const_iterator itINP = inputs.begin();
            while (itINP != inputs.end())
            {
                map.insert(map.end(), pair<IFMT, string>(ifmt, *itINP));
                itINP++;
            }
            itFMT++;
        }
        _gmutex.unlock();
        return map;
    }

    string t_gsetinp::basepath()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child_value("basepath");
        str_erase(tmp);
        _gmutex.unlock();
        return tmp;
    }

    string t_gsetinp::tbl()
    {
        _gmutex.lock();
        string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child_value("tbl");
        str_erase(tmp);
        _gmutex.unlock();
        return tmp;
    }

    void t_gsetinp::auto_discover(int year, int doy, int gpsWeek, int dow,
                                  const set<string> &sites)
    {
        _gmutex.lock();

        string bp = basepath();
        string tblPath = tbl();
        if (bp.empty() && tblPath.empty())
        {
            _gmutex.unlock();
            return;
        }

        xml_node inpNode = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP);

        // Helper: check if a child node exists and has non-empty text content
        auto has_content = [&](const string &name) -> bool
        {
            for (xml_node n = inpNode.first_child(); n; n = n.next_sibling())
            {
                if (string(n.name()) == name)
                {
                    string val = n.child_value();
                    str_erase(val);
                    if (val == "0") return false; // auto-discovery marker, treat as empty
                    return !val.empty();
                }
            }
            return false;
        };

        // Helper: add discovered path as a new child node
        auto add_to_xml = [&](const string &name, const string &path)
        {
            if (path.empty()) return;
            string p = path;
            if (p.find("://") == string::npos)
                p = GFILE_PREFIX + p;
            xml_node n = inpNode.append_child(name.c_str());
            n.append_child(node_pcdata).set_value(p.c_str());
        };

        // RINEXO: find obs for each site (or all obs if no site restriction)
        if (!has_content("rinexo"))
        {
            if (!sites.empty())
            {
                for (const string &site : sites)
                {
                    string path = gnut::findObsFile(bp, site, year, doy);
                    if (path.empty() && !tblPath.empty()) path = gnut::findObsFile(tblPath, site, year, doy);
                    add_to_xml("rinexo", path);
                }
            }
            else
            {
                // No site restriction (e.g., rec=0): discover all obs files
                vector<string> paths = gnut::findObsFiles(bp, "", year, doy);
                if (paths.empty() && !tblPath.empty()) paths = gnut::findObsFiles(tblPath, "", year, doy);
                for (const string &p : paths) add_to_xml("rinexo", p);
            }
        }

        // RINEXN
        if (!has_content("rinexn"))
        {
            string path = gnut::findEphFile(bp, year, doy);
            if (path.empty() && !tblPath.empty()) path = gnut::findEphFile(tblPath, year, doy);
            add_to_xml("rinexn", path);
        }

        // SP3
        if (!has_content("sp3"))
        {
            string path = gnut::findSp3File(bp, year, doy, gpsWeek, dow);
            if (path.empty() && !tblPath.empty()) path = gnut::findSp3File(tblPath, year, doy, gpsWeek, dow);
            add_to_xml("sp3", path);
        }

        // RINEXC (CLK)
        if (!has_content("rinexc"))
        {
            string path = gnut::findClkFile(bp, year, doy, gpsWeek, dow);
            if (path.empty() && !tblPath.empty()) path = gnut::findClkFile(tblPath, year, doy, gpsWeek, dow);
            add_to_xml("rinexc", path);
        }

        // EOP/ERP
        if (!has_content("eop"))
        {
            string path = gnut::findErpFile(bp, year, doy);
            if (path.empty() && !tblPath.empty()) path = gnut::findErpFile(tblPath, year, doy);
            add_to_xml("eop", path);
        }

        // BIAS / BIASINEX
        if (!has_content("bias") && !has_content("biasinex"))
        {
            vector<string> paths = gnut::findBiasFiles(bp, year, doy);
            if (paths.empty() && !tblPath.empty()) paths = gnut::findBiasFiles(tblPath, year, doy);
            for (const string &p : paths)
            {
                string fname = gnut::base_name(p);
                size_t dot = fname.rfind('.');
                if (dot != string::npos)
                {
                    string ext = fname.substr(dot);
                    transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                    // .bsx files are Bias-SINEX (Bernese) format and should be parsed
                    // by t_biabernese (BIAS_INP), not t_biasinex (BIASINEX_INP).
                    // t_biabernese correctly sets the _A/_R suffix for _ac, which is
                    // required by gallbias AC priority lookup.
                    add_to_xml("bias", p);
                }
                else
                {
                    add_to_xml("bias", p);
                }
            }
        }

        // IFCB (optional, only auto-discover if frequency >= 3 and GPS is used)
        // Skip auto-discovery for IFCB to avoid compatibility issues with 2-frequency setups
        // if (!has_content("ifcb"))
        // {
        //     string path = gnut::findSGGIfcb(bp, year, doy);
        //     if (path.empty() && !tblPath.empty()) path = gnut::findSGGIfcb(tblPath, year, doy);
        //     add_to_xml("ifcb", path);
        // }

        // UPD: search common types and all systems
        if (!has_content("upd"))
        {
            vector<string> upd_types = {"ewl", "wl", "nl"};
            for (const string &t : upd_types)
            {
                vector<string> paths = gnut::findUPDFiles(bp, year, doy, t);
                if (paths.empty() && !tblPath.empty()) paths = gnut::findUPDFiles(tblPath, year, doy, t);
                for (const string &p : paths)
                {
                    add_to_xml("upd", p);
                }
            }
        }

        // ATX
        if (!has_content("atx"))
        {
            string path = gnut::findATXFile(bp);
            if (path.empty() && !tblPath.empty()) path = gnut::findATXFile(tblPath);
            add_to_xml("atx", path);
        }

        // BLQ
        if (!has_content("blq"))
        {
            string path = gnut::findBLQFile(bp);
            if (path.empty() && !tblPath.empty()) path = gnut::findBLQFile(tblPath);
            add_to_xml("blq", path);
        }

        // DE
        if (!has_content("de"))
        {
            string path = gnut::findDEFile(bp);
            if (path.empty() && !tblPath.empty()) path = gnut::findDEFile(tblPath);
            add_to_xml("de", path);
        }

        // LEAPSECOND
        if (!has_content("leapsecond"))
        {
            string path = gnut::findLeapFile(bp);
            if (path.empty() && !tblPath.empty()) path = gnut::findLeapFile(tblPath);
            add_to_xml("leapsecond", path);
        }

        _gmutex.unlock();
    }

    vector<string> t_gsetinp::inputs(const string &fmt)
    {
        IFMT ifmt = str2ifmt(fmt);
        return _inputs(ifmt);
    }

    vector<string> t_gsetinp::inputs(const IFMT &ifmt)
    {
        return _inputs(ifmt);
    }

    vector<string> t_gsetinp::_inputs(const string &fmt)
    {
        vector<string> tmp;
        set<string> list;
        string str;

        for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
        {
            if (node.name() == fmt)
            {
                string val = node.child_value();
                size_t s = val.find_first_not_of(" \t\n\r");
                if (s == string::npos) continue;
                size_t e = val.find_last_not_of(" \t\n\r");
                string trimmed = val.substr(s, e - s + 1);

                istringstream is(trimmed);
                if (is >> str && !is.fail())
                {
                    if (str == "0") continue; // auto-discovery marker, skip
                    if (str == "1")
                    {
                        // Use remaining tokens as explicit paths
                        while (is >> str && !is.fail())
                        {
                            if (str.find("://") == string::npos)
                                str = GFILE_PREFIX + str;
                            if (list.find(str) == list.end())
                            {
                                tmp.push_back(str);
                                list.insert(str);
                            }
                            else
                            {
                                cout << "READ : " + str + " multiple request ignored" << endl;
                            }
                        }
                        continue;
                    }
                    // Old format: process full content with str_erase
                    // str_erase(val);  // BUG: removes ALL whitespace, merging multiple paths into one
                    istringstream is_old(trimmed);  // Use trimmed (whitespace-trimmed only) to preserve internal spaces as token separators
                    while (is_old >> str && !is_old.fail())
                    {
                        if (str.find("://") == string::npos)
                            str = GFILE_PREFIX + str;
                        if (list.find(str) == list.end())
                        {
                            tmp.push_back(str);
                            list.insert(str);
                        }
                        else
                        {
                            cout << "READ : " + str + " multiple request ignored" << endl;
                        }
                    }
                }
            }
        }
        return tmp;
    }

    vector<string> t_gsetinp::_inputs(const IFMT &fmt)
    {
        vector<string> tmp;
        set<string> list;
        string str;

        for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
        {
            // jdhuang
            IFMT ifmt = str2ifmt(node.name());
            if (ifmt == IFMT::UNDEF)
                continue;
            if (ifmt == fmt)
            {
                string val = node.child_value();
                size_t s = val.find_first_not_of(" \t\n\r");
                if (s == string::npos) continue;
                size_t e = val.find_last_not_of(" \t\n\r");
                string trimmed = val.substr(s, e - s + 1);

                istringstream is(trimmed);
                if (is >> str && !is.fail())
                {
                    if (str == "0") continue;
                    if (str == "1")
                    {
                        while (is >> str && !is.fail())
                        {
                            if (str.find("://") == string::npos)
                                str = GFILE_PREFIX + str;
                            if (list.find(str) == list.end())
                            {
                                tmp.push_back(str);
                                list.insert(str);
                            }
                            else
                            {
                                cout << "READ : " + str + " multiple request ignored" << endl;
                            }
                        }
                        continue;
                    }
                    // Old format
                    str_erase(val);
                    istringstream is_old(val);
                    while (is_old >> str && !is_old.fail())
                    {
                        if (str.find("://") == string::npos)
                            str = GFILE_PREFIX + str;
                        if (list.find(str) == list.end())
                        {
                            tmp.push_back(str);
                            list.insert(str);
                        }
                        else
                        {
                            cout << "READ : " + str + " multiple request ignored" << endl;
                        }
                    }
                }
            }
        }
        return tmp;
    }

    set<string> t_gsetinp::_iformats()
    {
        set<string> tmp;
        for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
        {
            tmp.insert(node.name());
        }
        return tmp;
    }

    void t_gsetinp::check()
    {
        _gmutex.lock();

        // check existence of nodes/attributes
        xml_node parent = _doc.child(XMLKEY_ROOT);
        xml_node node = _default_node(parent, XMLKEY_INP);

        // check supported input formats (see IFMT enum !)
        set<string> ifmt = _iformats();
        set<string>::const_iterator itFMT = ifmt.begin();
        while (itFMT != ifmt.end())
        {
            string fmt = *itFMT;
            if (fmt == "basepath" || fmt == "tbl")
            {
                itFMT++;
                continue;
            }
            try
            {
                str2ifmt(fmt);
            }
            catch (const std::exception &e)
            {
                if (fmt == "basepath" || fmt == "tbl")
                {
                    itFMT++;
                    continue;
                }
                _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
                itFMT++;
                continue;
            }
            // check application-specific output format
            if (_IFMT_supported.find(str2ifmt(fmt)) == _IFMT_supported.end())
            {
                _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
                spdlog::warn(fmt + " inp format not supported by this application!");
            }
            // Skip empty check for basepath itself; check only real input formats
            if (fmt != "basepath")
                check_input(fmt, "your fmt : " + fmt + " is empty");
            itFMT++;
        }

        _default_attr(node, "chk_nav", _chkNavig);
        _default_attr(node, "chk_health", _chkHealth);

        xml_node nodeBNCRTCM = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child("bncrtcm");
        _default_attr(nodeBNCRTCM, "_corrStream", _corrStream);

        _gmutex.unlock();
        return;
    }

    void t_gsetinp::help()
    {
        _gmutex.lock();

        cerr << " <inputs>\n"
             << "   <basepath> file://dir/name </basepath> \t\t <!-- base directory for auto-discovery -->\n"
             << "   <tbl> file://dir/name </tbl> \t\t <!-- fallback directory for model files -->\n"
             << "   <rinexo> file://dir/name </rinexo> \t\t <!-- obs RINEX decoder -->\n"
             << "   <rinexn> file://dir/name </rinexn> \t\t <!-- nav RINEX decoder -->\n"
             << " </inputs>\n";

        cerr << "\t<!-- inputs description:\n"
             << "\t <basepath> dir </basepath>           <!-- base directory for auto file discovery -->\n"
             << "\t <tbl> dir </tbl>                     <!-- fallback directory for model/table files -->\n"
             << "\t <decoder> path1 path2 path3  </decoder>\n"
             << "\t ... \n"
             << "\t where path(i) contains [file,tcp,ntrip]:// depending on the application\n"
             << "\t empty or missing tags trigger auto-discovery when basepath/tbl is set\n"
             << "\t -->\n\n";

        _gmutex.unlock();
        return;
    }

} // namespace
