/*
*
* @verbatim
    History

    @endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

*
* @file        gsetrec.h
* @brief       implements receiver object setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#ifndef GSETREC_H
#define GSETREC_H

#define XMLKEY_REC "receiver"

#include <map>
#include <string>
#include <iostream>
#include <memory>

#include "gdata/grec.h"
#include "gset/gsetbase.h"
#include "gset/gtruecrd.h"
#include "gutils/gtypeconv.h"
#include "gutils/gtriple.h"
#include "gmodels/ggpt.h"

#define HSL_UNKNOWN -9999 

using namespace std;

namespace gnut
{

    class LibGnut_LIBRARY_EXPORT t_gsetrec : public virtual t_gsetbase
    {
    public:
        /** @brief constructor */
        t_gsetrec();

        /** @brief destructor */
        ~t_gsetrec();

        /** @brief settings check */
        void check(); 

        /** @brief settings help */
        void help(); 

        /** @brief get crd xyz */
        t_gtriple get_crd_xyz(string s);

        /** @brief get crd xyz from true_crd file (if available) */
        t_gtriple get_true_crd_xyz(string s);

        /** @brief check if true_crd file is loaded */
        bool has_true_crd() const;

        /** @brief get all objects IDs */
        set<string> objects();                                  
        shared_ptr<t_grec> grec(string s, t_spdlog spdlog = 0); 

        /**
         * @brief get the List of recevier names
         * @return set<string> : List of recevier names
         */
        virtual set<string> recs();
        set<string> all_rec();

    protected:
        /** @brief get crd xyz */
        t_gtriple _get_crd_xyz(string s);

        /** @brief get ecc neu */
        t_gtriple _get_ecc_neu(string s);

        /** @brief get crd blh */
        t_gtriple _get_crd_blh(string s); 

        /** @brief Global Pressure Temperature model */
        t_gpt _ggpt; 

        /** @brief get all objects IDs */
        set<string> _objects();

        /** @brief load true_crd file from XML config */
        void _load_true_crd();

        string _rec;      // default receiver name
        string _id;       // receiver id
        string _name_rec; // receiver name
        t_gtime _beg;     // default begin time
        t_gtime _end;     // default end time
        double _X;        // receiver X-coordinate [m]
        double _Y;        // receiver Y-coordinate [m]
        double _Z;        // receiver Z-coordinate [m]
        bool _overwrite;

        shared_ptr<t_gtruecrd> _true_crd_reader;  ///< true_crd file reader
        bool _true_crd_enabled;                   ///< whether true_crd is enabled globally

    private:
    };

} // namespace

#endif
