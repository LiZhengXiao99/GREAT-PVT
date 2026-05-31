
#include "gcfg_ppp.h"
#include "gutils/gfileconv.h"
#include <chrono>
#include <thread>

using namespace std;
using namespace gnut;
using namespace gsins;
using namespace std::chrono;

void catch_signal(int) { cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; }

// MAIN
// ----------
int main(int argc, char** argv)
{
    // Only to cout the Reminder here
    signal(SIGINT, catch_signal);

    // Construct the gset class and init some values in the class
    t_gcfg_ppp gset;
    gset.app("GRAET-PVT", "$Ver: 1.1 $", "$Rev:  $", "(https://github.com/GREAT-WHU/GREAT-PVT)", __DATE__, __TIME__);

    // Get the arguments from the command line
    gset.arg(argc, argv, true, false);

    // Creat and set the log file : ppp.log
    auto log_type = dynamic_cast<t_gsetout *>(&gset)->log_type();
    auto log_level = dynamic_cast<t_gsetout *>(&gset)->log_level();
    auto log_name = dynamic_cast<t_gsetout *>(&gset)->log_name();
    auto log_pattern = dynamic_cast<t_gsetout *>(&gset)->log_pattern();
    spdlog::set_level(log_level);
    spdlog::set_pattern(log_pattern);
    spdlog::flush_on(spdlog::level::err);
    t_grtlog great_log = t_grtlog(log_type, log_level, log_name);
    auto my_logger = great_log.spdlog();

    // Handle command-line observation file (-o): override XML rinexo and rec
    string obsfile = gset.obsfile();
    if (!obsfile.empty())
    {
        // Extract site name from filename (first 4 chars, uppercased)
        string fname = gnut::base_name(obsfile);
        string site = fname.substr(0, min((size_t)4, fname.size()));
        transform(site.begin(), site.end(), site.begin(), ::toupper);

        xml_node config = gset.config_node();

        auto clear_and_set = [](xml_node parent, const char* name, const char* val)
        {
            xml_node node = parent.child(name);
            if (!node) node = parent.append_child(name);
            for (xml_node child = node.first_child(); child; )
            {
                xml_node next = child.next_sibling();
                node.remove_child(child);
                child = next;
            }
            node.append_child(pugi::node_pcdata).set_value(val);
        };

        // Override <inputs><rinexo>
        xml_node inpNode = config.child("inputs");
        if (!inpNode) inpNode = config.append_child("inputs");
        clear_and_set(inpNode, "rinexo", obsfile.c_str());

        // Override <gen><rec>
        xml_node genNode = config.child("gen");
        if (!genNode) genNode = config.append_child("gen");
        clear_and_set(genNode, "rec", site.c_str());

        SPDLOG_LOGGER_INFO(my_logger,
            "Command-line obs file override: site=" + site + " file=" + obsfile);
    }

    // Check the base station
    bool isBase = dynamic_cast<t_gsetgen*>(&gset)->list_base().size();

    // Prepare site list from gset
    set<string> sites = dynamic_cast<t_gsetgen*>(&gset)->recs();

    // Auto-discover missing input files from basepath/tbl
    t_gtime beg_time = dynamic_cast<t_gsetgen*>(&gset)->beg();
    t_gsetinp* gsetinp = dynamic_cast<t_gsetinp*>(&gset);
    if (!gsetinp->basepath().empty() || !gsetinp->tbl().empty())
    {
        int year = beg_time.year();
        int doy = beg_time.doy();
        // If beg=0 (FIRST_TIME), use 0/0 to trigger broad file discovery
        if (beg_time == FIRST_TIME) { year = 0; doy = 0; }
        gsetinp->auto_discover(
            year, doy, beg_time.gwk(), beg_time.dow(), sites);
    }

    // Prepare input files list from gset
    multimap<IFMT, string> inp = gset.inputs_all();

    // Get sample intval from gset. if not, init with the default value
    int sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling());
    if (!sample) sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling_default());

    //--- INITIALIZATIONS ---
    t_gdata* gdata = nullptr;
    t_gnavde* gde = new t_gnavde;
    t_gpoleut1* gerp = new t_gpoleut1;
    t_gallobs* gobs = new t_gallobs();  gobs->spdlog(my_logger); gobs->gset(&gset);
    t_gallprec* gorb = new t_gallprec(); gorb->spdlog(my_logger);
    t_gallpcv* gpcv = nullptr; if (gset.input_size("atx") > 0) { gpcv = new t_gallpcv;  gpcv->spdlog(my_logger); }
    t_gallotl* gotl = nullptr; if (gset.input_size("blq") > 0) { gotl = new t_gallotl;  gotl->spdlog(my_logger); }
    t_gallbias* gbia = nullptr; if (gset.input_size("biasinex") > 0 ||gset.input_size("bias") > 0) { gbia = new t_gallbias; gbia->spdlog(my_logger); }
    t_gallobj* gobj = new t_gallobj(my_logger, gpcv, gotl); gobj->spdlog(my_logger);
    t_gupd* gupd = nullptr; if (gset.input_size("upd") > 0) { gupd = new t_gupd;  gupd->spdlog(my_logger); }
    t_gifcb* gifcb = nullptr;  if (gset.input_size("ifcb") > 0) { gifcb = new t_gifcb;  gifcb->spdlog(my_logger); }

    // vgppp for the process of ppp with filter
    vector<t_gpvtflt*> vgpvt;

    // runepoch for the time costed each epoch 
    t_gtime runepoch(t_gtime::GPS);

    // lstepoch for the time of all epoches 
    t_gtime lstepoch(t_gtime::GPS);

    // check precision products
    if (gset.input_size("sp3") == 0 && gset.input_size("rinexc") == 0)
    {
        gorb->use_clknav(true);
        gorb->use_posnav(true);
    }
    else if (gset.input_size("sp3") > 0 && gset.input_size("rinexc") == 0)
    {
        gorb->use_clksp3(true);
    }

    // SET OBJECTS
    set<string>::const_iterator itOBJ;
    set<string> obj = dynamic_cast<t_gsetrec*>(&gset)->objects();
    for (itOBJ = obj.begin(); itOBJ != obj.end(); ++itOBJ) 
    {
        string name = *itOBJ;
        shared_ptr<t_grec> rec = dynamic_cast<t_gsetrec*>(&gset)->grec(name, my_logger);
        gobj->add(rec);
    }

    // Multi gcoder for multi-thread decoding data
    vector<t_gcoder*> gcoder;

    // Multi gior for multi-thread receiving data
    vector<t_gio*> gio;

    // multi-thread 
    vector<thread> gthread;

    t_gio* tgio = 0;
    t_gcoder* tgcoder = 0;

    if (!isBase)
    {
        // CHECK INPUTS, sp3+rinexc+rinexo, Necessary data
        if (gset.input_size("sp3") == 0 &&
            gset.input_size("rinexc") == 0 &&
            gset.input_size("rinexo") == 0
            ) 
        {
            SPDLOG_LOGGER_INFO(my_logger, "Error: incomplete input: rinexo + rinexc + sp3");
            gset.usage();
        }
    }

    // DATA READING
    multimap<IFMT, string>::const_iterator itINP = inp.begin();
    for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
    {
        // Get the file format/path, which will be used in decoder
        IFMT   ifmt(itINP->first);
        string path(itINP->second);
        string id("ID" + int2str(i));

        // For different file format, we prepare different data container and decoder for them.
        if (ifmt == IFMT::RINEXO_INP) { gdata = gobs; tgcoder = new t_rinexo(&gset, "", 4096); }
        else if (ifmt == IFMT::SP3_INP) { gdata = gorb; tgcoder = new t_sp3(&gset, "", 8172); }
        else if (ifmt == IFMT::RINEXC_INP) { gdata = gorb; tgcoder = new t_rinexc(&gset, "", 4096); }
        else if (ifmt == IFMT::RINEXN_INP) { gdata = gorb; tgcoder = new t_rinexn(&gset, "", 4096); }
        else if (ifmt == IFMT::ATX_INP) { gdata = gpcv; tgcoder = new t_atx(&gset, "", 4096); }
        else if (ifmt == IFMT::BLQ_INP) { gdata = gotl; tgcoder = new t_blq(&gset, "", 4096); }
        else if (ifmt == IFMT::UPD_INP) { gdata = gupd; tgcoder = new t_upd(&gset, "", 4096); }
        else if (ifmt == IFMT::BIASINEX_INP) { gdata = gbia; tgcoder = new t_biasinex(&gset, "", 20480); }
        else if (ifmt == IFMT::BIAS_INP) { gdata = gbia; tgcoder = new t_biabernese(&gset, "", 20480); }
        else if (ifmt == IFMT::DE_INP) { gdata = gde; tgcoder = new t_dvpteph405(&gset, "", 4096); }
        else if (ifmt == IFMT::EOP_INP) { gdata = gerp; tgcoder = new t_poleut1(&gset, "", 4096); }
        else if (ifmt == IFMT::IFCB_INP) { gdata = gifcb; tgcoder = new t_ifcb(&gset, "", 4096); }
        else 
        {
            SPDLOG_LOGGER_INFO(my_logger, "Error: unrecognized format " + int2str(int(ifmt)));
            gdata = 0;
        }

        // Check the file path
        if (path.substr(0, 7) == "file://") 
        {
            SPDLOG_LOGGER_INFO(my_logger, "path is file!");
            tgio = new t_gfile(my_logger);
            tgio->spdlog(my_logger);
            tgio->path(path);
        }

        // READ DATA FROM FILE
        if (tgcoder) 
        {
            // Put the file into gcoder
            tgcoder->clear();
            tgcoder->path(path);
            tgcoder->spdlog(my_logger);

            // Put the data container into gcoder
            tgcoder->add_data(id, gdata);
            tgcoder->add_data("OBJ", gobj); 

            // Put the gcoder into the gio
            // Note, gcoder contain the gdata and gio contain the gcoder
            tgio->coder(tgcoder);

            runepoch = t_gtime::current_time(t_gtime::GPS);

            // Read the data from file here
            tgio->run_read();
            lstepoch = t_gtime::current_time(t_gtime::GPS);

            // Write the information of reading process to log file
            SPDLOG_LOGGER_INFO(my_logger, "READ: " + path + " time: " + dbl2str(lstepoch.diff(runepoch)) + " sec");

            // Delete 
            delete tgio;
            delete tgcoder;

        }
    }
    // Determine actual observation time range when beg/end are set to 0
    t_gtime actual_beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
    t_gtime actual_end = dynamic_cast<t_gsetgen*>(&gset)->end();
    if (actual_beg == FIRST_TIME && gobs)
    {
        actual_beg = LAST_TIME;
        for (const string& s : gobs->stations())
        {
            t_gtime t = gobs->beg_obs(s);
            if (t != FIRST_TIME && t < actual_beg) actual_beg = t;
        }
        if (actual_beg == LAST_TIME) actual_beg = FIRST_TIME;
    }
    if (actual_end == LAST_TIME && gobs)
    {
        actual_end = FIRST_TIME;
        for (const string& s : gobs->stations())
        {
            t_gtime t = gobs->end_obs(s);
            if (t != LAST_TIME && t > actual_end) actual_end = t;
        }
        if (actual_end == FIRST_TIME) actual_end = LAST_TIME;
    }
    // Update beg_time for output placeholder substitution
    if (beg_time == FIRST_TIME && actual_beg != FIRST_TIME)
    {
        beg_time = actual_beg;
    }

    // set antennas for satllites (must be before PCV assigning)
    gobj->read_satinfo(beg_time);

    // assigning PCV pointers to objects
    gobj->sync_pcvs();

    // add all data
    t_gallproc* data = new t_gallproc();
    if (gobs)data->Add_Data(t_gdata::type2str(gobs->id_type()), gobs);
    if (gorb)data->Add_Data(t_gdata::type2str(gorb->id_type()), gorb);
    if (gobj)data->Add_Data(t_gdata::type2str(gobj->id_type()), gobj);
    if (gbia)data->Add_Data(t_gdata::type2str(gbia->id_type()), gbia);
    if (gotl)data->Add_Data(t_gdata::type2str(gotl->id_type()), gotl);
    if (gde)data->Add_Data(t_gdata::type2str(gde->id_type()), gde);
    if (gerp)data->Add_Data(t_gdata::type2str(gerp->id_type()), gerp);
    if (gupd && dynamic_cast<t_gsetamb*>(&gset)->fix_mode() != FIX_MODE::NO && !isBase)
    {
        data->Add_Data(t_gdata::type2str(gupd->id_type()), gupd);
    }
    int frequency = dynamic_cast<t_gsetproc*>(&gset)->frequency();
    set<string> system = dynamic_cast<t_gsetgen*>(&gset)->sys();
    if (frequency >= 3 && system.find("GPS") != system.end() && !isBase)
    {
        if (gifcb)data->Add_Data(t_gdata::type2str(gifcb->id_type()), gifcb);
    }
    
    // If no sites configured (e.g., rec=0), use all stations from observations
    if (!isBase && sites.empty() && gobs)
    {
        sites = gobs->stations();
    }

    // ============================================
    // CRITICAL INPUT VALIDATION
    // Abort early if mandatory data is missing to avoid
    // segfaults downstream (e.g. null _trs2crs_2000).
    // ============================================
    bool critical_missing = false;

    // 1. Observations
    if (!gobs || gobs->stations().empty())
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "CRITICAL: No observation data loaded! Check rinexo input.");
        critical_missing = true;
    }

    // 2. Ephemeris (SP3 / CLK / broadcast)
    if (!gorb || gorb->satellites().empty())
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "CRITICAL: No ephemeris data loaded! Check sp3 / rinexc / rinexn input.");
        critical_missing = true;
    }

    // 3. ERP / PoleUT1 (required by t_gprecisebias::_update_rot_matrix)
    if (!gerp || gerp->isEmpty())
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "CRITICAL: No ERP/poleut1 data loaded! Check eop / poleut1 input (also in tbl/).");
        critical_missing = true;
    }

    // 4. Ocean load (BLQ) – only if explicitly configured
    if (gset.input_size("blq") > 0 && (!gotl))
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "CRITICAL: BLQ/oceanload configured but not loaded!");
        critical_missing = true;
    }

    // 5. ATX – only if explicitly configured
    if (gset.input_size("atx") > 0 && (!gpcv))
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "CRITICAL: ATX configured but not loaded!");
        critical_missing = true;
    }

    if (critical_missing)
    {
        SPDLOG_LOGGER_CRITICAL(my_logger, "Missing critical input files. Aborting.");
        return -1;
    }

    // Record current time
    auto tic_start = system_clock::now();

    // PVT PROCESSING - loop over sites from settings
    int i = 0;
    if (isBase)
    {
        vector<string> list_rover = gset.list_rover();  
        sites = set<string>(list_rover.begin(), list_rover.end());
    }
    int nsite = sites.size();
    set<string>::iterator it = sites.begin();

    // Prepare context strings for output placeholder substitution
    string sys_str;
    if (system.find("GPS") != system.end()) sys_str += "G";
    if (system.find("GLO") != system.end()) sys_str += "R";
    if (system.find("GAL") != system.end()) sys_str += "E";
    if (system.find("BDS") != system.end()) sys_str += "C";
    if (system.find("QZS") != system.end()) sys_str += "J";
    if (system.find("SBS") != system.end()) sys_str += "S";
    if (sys_str.empty()) sys_str = "GREC";

    t_gsetproc* gsetproc = dynamic_cast<t_gsetproc*>(&gset);
    OBSCOMBIN obs_combin = gsetproc->obs_combin();
    string mode_str;
    if (gsetproc->pos_kin()) mode_str = "KIN";
    else if (gsetproc->crd_est() == CONSTRPAR::FIX) mode_str = "FIX";
    else mode_str = "EST";

    string iono_str;
    if (obs_combin == OBSCOMBIN::IONO_FREE)
        iono_str = "IF";
    else
        iono_str = "UC";

    while (i < nsite)
    {
        string site_base = "";
        string site = *it;
        // Check site data
        if (isBase)
        {
            site_base = (gset.list_base())[i];
            site = (gset.list_rover())[i];
            if (gobs->beg_obs(site_base) == LAST_TIME || gobs->end_obs(site_base) == FIRST_TIME ||
                 site_base.empty() || gobs->isSite(site_base) == false) 
            {
                SPDLOG_LOGGER_INFO(my_logger, "No two site/data for processing!");
                i++;
                continue;
            }
        }
        if (gobs->beg_obs(site) == LAST_TIME ||
            gobs->end_obs(site) == FIRST_TIME ||
            site.empty() || gobs->isSite(site) == false) 
        {
            SPDLOG_LOGGER_INFO(my_logger, "No site/data for processing!");
            if (!isBase) it++;
            i++;
            continue;
        }

        // Set output context for placeholder substitution
        dynamic_cast<t_gsetout*>(&gset)->set_context(site, beg_time.year(), beg_time.doy(),
                                                     mode_str, iono_str, sys_str);

        // Add site data
        vgpvt.push_back(0); int idx = vgpvt.size() - 1;
        vgpvt[idx] = new t_gpvtflt(site, site_base, &gset, my_logger, data);
        if (dynamic_cast<t_gsetamb*>(&gset)->fix_mode() != FIX_MODE::NO && !isBase)
        {
            vgpvt[idx]->Add_UPD(gupd);
        }

        SPDLOG_LOGGER_INFO(my_logger, "PVT processing started ");
        SPDLOG_LOGGER_INFO(my_logger, actual_beg.str_ymdhms("  beg: ") + actual_end.str_ymdhms("  end: "));

        runepoch = t_gtime::current_time(t_gtime::GPS);

        // The main processing code : processBatch
        vgpvt[idx]->processBatch(actual_beg, actual_end, true);

        // The time when process ends
        lstepoch = t_gtime::current_time(t_gtime::GPS);

        // Write the log file
        SPDLOG_LOGGER_INFO(my_logger, site_base + site + "PVT processing finished : duration  " + dbl2str(lstepoch.diff(runepoch)) + " sec");

        if (!isBase) it++;
        i++;
    }

    //Delete pointer
    for (size_t i = 0; i < gio.size(); ++i) { delete gio[i]; }; gio.clear();
    for (size_t i = 0; i < gcoder.size(); ++i) { delete gcoder[i]; }; gcoder.clear();
    for (unsigned int i = 0; i < vgpvt.size(); ++i) { if (vgpvt[i])  delete vgpvt[i]; }
    if (gobs) delete gobs;
    if (gpcv) delete gpcv;
    if (gotl) delete gotl;
    if (gobj) delete gobj;
    if (gorb) delete gorb;
    if (gbia) delete gbia;
    if (gde)  delete gde;
    if (gerp)  delete gerp;
    if (gupd)  delete gupd;
    if (gifcb)  delete gifcb;
    if (data)  delete data;

    // Record the current time and the time spent
    auto tic_end = system_clock::now();
    auto duration = duration_cast<microseconds>(tic_end - tic_start);
    cout << "Spent" << double(duration.count()) * microseconds::period::num / microseconds::period::den << " seconds." << endl;

    return 0;
}