/**
*
* @verbatim
    History

    @endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

*
* @file        gsetout.h
* @brief       implements output setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#ifndef GSETOUT_H
#define GSETOUT_H

#define XMLKEY_OUT "outputs"   ///< The defination of outputs module in xml file

#define DEFAULT_FILE_VER ""    ///< default version for file format
#define DEFAULT_FILE_UPD 0     ///< default auto update [min] for file saving
#define DEFAULT_FILE_LEN 0     ///< default auto length [min] for file saving
#define DEFAULT_FILE_SMP 0     ///< default auto sample [sec] for file saving
#define DEFAULT_FILE_OFF 0     ///< default file offset [min] for file saving
#define DEFAULT_FILE_SYS "UTC" ///< default file system       for file saving

#include <string>
#include <iostream>

#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"
#include "spdlog/spdlog.h"

using namespace std;
using namespace pugi;
using namespace spdlog;

namespace gnut
{
    enum OFMT
    {
        XXX_OUT,
        LOG_OUT,
        PPP_OUT,
        FLT_OUT,
        FLT_FLOAT_OUT,
        AUG_OUT,
        FLT_PPPRTK_OUT
    };

    class LibGnut_LIBRARY_EXPORT t_gsetout : public virtual t_gsetbase
    {
    public:
        /**@brief defalut constructor */
        t_gsetout();
        /**@brief defalut destructor */
        ~t_gsetout();

        /**
         * @brief change from string to OFMT
         * @param[in] s file format
         * @return OFMT : file format
         */
        static OFMT str2ofmt(const string &s);

        /**
         * @brief change from OFMT to string
         * @param[in] f file format
         * @return string : file format
         */
        static string ofmt2str(const OFMT &f);

        /**@brief settings check */
        void check();

        /**@brief settings help */
        void help();

        // attributes
        /**
         * @brief  get verbosity attribute
         * @return int : verbosity attribute
         */
        int verb();

        /**
         * @brief  get append request
         * @return bool : append request
         */
        bool append();

        // elements
        /**
         * @brief  get format output size
         * @param[in] fmt file format
         * @return int : format output size
         */
        int output_size(const string &fmt);

        /**
         * @brief  get string outputs
         * @param[in] fmt file format
         * @return string : string outputs
         */
        string outputs(const string &fmt);

        /**
         * @brief  get string outputs with placeholder substitution
         * @param[in] fmt file format
         * @param[in] site station/site name
         * @param[in] year 4-digit year
         * @param[in] doy day of year
         * @param[in] mode_str mode string (e.g. "PPP", "FLT", "AUG")
         * @param[in] iono_str iono string (e.g. "DF", "FF", "IF")
         * @param[in] sys_str system string (e.g. "GREC")
         * @return string : string outputs with placeholders replaced
         *
         * Supported placeholders:
         *   $(MODE)    -> mode_str
         *   $(IONO)    -> iono_str
         *   $(SYSTEM)  -> sys_str
         *   $(REC)     -> site
         *   $(STATION) -> site
         *   $(YEAR)    -> 4-digit year
         *   $(DOY)     -> 3-digit DOY (zero-padded)
         *   $(EXT)     -> file extension derived from fmt
         */
        string outputs(const string &fmt, const string &site,
                       int year, int doy,
                       const string &mode_str,
                       const string &iono_str,
                       const string &sys_str);

        /**
         * @brief substitute placeholders in a path string
         */
        static string substitute_placeholders(const string &path,
                                               const string &site,
                                               int year, int doy,
                                               const string &mode_str,
                                               const string &iono_str,
                                               const string &sys_str);

        /*
         * @brief  get string log type
        */
        string log_type();

        /*
         * @brief  get string log name
        */
        string log_name();

        /*
         * @brief  get string log level
        */
        level::level_enum log_level();

        /*
        * @brief  get string log pattern
        */
        string log_pattern();

        /**
         * @brief  get formats
         * @return set<string> : all the outputs
         */
        set<string> oformats();

        /**
         * @brief  get string output version
         * @param[in] fmt file format
         * @return string : string output version
         */
        string version(const string &fmt);

        /**
         * @brief set context for placeholder substitution
         * @param[in] site station name
         * @param[in] year 4-digit year
         * @param[in] doy day of year
         * @param[in] mode mode string
         * @param[in] iono iono string
         * @param[in] sys system string
         */
        void set_context(const string &site, int year, int doy,
                         const string &mode, const string &iono, const string &sys);

        /**
         * @brief  get whether context has been set
         * @return bool : true if context has been set
         */
        bool ctx_set() const;

        /**
         * @brief  get context year
         * @return int : 4-digit year (0 if not set)
         */
        int ctx_year() const;

        /**
         * @brief  get context doy
         * @return int : day of year (0 if not set)
         */
        int ctx_doy() const;

    protected:
        /**
         * @brief  get string output file
         * @param[in] fmt file format
         * @return string : string output file ("AUTO" for auto-generated path)
         */
        string _outputs(const string &fmt);

        /**
         * @brief  generate default output path for auto mode
         * @param[in] fmt file format
         * @param[in] site station name
         * @param[in] year 4-digit year
         * @param[in] doy day of year
         * @param[in] mode_str mode string
         * @param[in] iono_str iono string
         * @param[in] sys_str system string
         * @return string : auto-generated output path
         */
        string _default_output_path(const string &fmt, const string &site,
                                    int year, int doy,
                                    const string &mode_str,
                                    const string &iono_str,
                                    const string &sys_str);

        /**
         * @brief  get all string output file
         * @return set<string> : all string output file
         */
        set<string> _oformats();

        set<OFMT> _OFMT_supported; ///< vector of supported OFMTs (app-specific)
        bool _append;              ///< append mode
        int _verb;                 ///< output verbosity

        // context for placeholder substitution
        string _ctx_site;          ///< station name
        int _ctx_year;             ///< 4-digit year
        int _ctx_doy;              ///< day of year
        string _ctx_mode;          ///< mode string
        string _ctx_iono;          ///< iono string
        string _ctx_sys;           ///< system string
        bool _ctx_set;             ///< whether context has been set

    private:
    };

} // namespace
#endif
