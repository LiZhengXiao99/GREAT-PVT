#!/dsk/lipan/programs/anaconda2/bin/python
# coding:utf-8

import os
import re
import shutil
import datetime
import math
import copy
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from functools import partial
from typing import List
from itertools import groupby
import multiprocessing
import time
import sys
from collections import Counter

# 字体设置
font_path ="/home/lzx/code/GREAT-PVT/script/arial.ttf"
font_prop =fm.FontProperties(fname=font_path)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False
mm = 1/25.4 # inch 和 毫米的转换

# ------------------------------------------------------------------
# Borrowed from plot_pos_error.py: true coordinate loading & ENU computation
# ------------------------------------------------------------------
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)

def ecef_to_geodetic(x, y, z):
    a = WGS84_A
    e2 = WGS84_E2
    lon = math.atan2(y, x)
    p = math.sqrt(x * x + y * y)
    lat = math.atan2(z, p * (1 - e2))
    for _ in range(10):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sin_lat * sin_lat)
        h = p / math.cos(lat) - N
        lat_new = math.atan2(z, p * (1 - e2 * (N / (N + h))))
        if abs(lat_new - lat) < 1e-12:
            lat = lat_new
            break
        lat = lat_new
    sin_lat = math.sin(lat)
    N = a / math.sqrt(1 - e2 * sin_lat * sin_lat)
    h = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), h

def ecef_to_enu_matrix(lat_deg, lon_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    sb, cb = math.sin(lon), math.cos(lon)
    R = np.array([[-sb, cb, 0],
                  [-sl*cb, -sl*sb, cl],
                  [cl*cb, cl*sb, sl]])
    return R

def load_true_coords(true_file):
    tbl = {}
    if not os.path.isfile(true_file):
        raise FileNotFoundError("True coordinate file not found: %s" % true_file)
    if true_file.lower().endswith('.snx'):
        coords_tmp = {}
        with open(true_file, 'r', encoding='utf-8', errors='ignore') as f:
            in_est = False
            for line in f:
                if line.startswith('+SOLUTION/ESTIMATE'):
                    in_est = True
                    continue
                if line.startswith('-SOLUTION/ESTIMATE'):
                    break
                if not in_est:
                    continue
                if not line.strip() or line.lstrip().startswith('*'):
                    continue
                parts = line.strip().split()
                if len(parts) < 9:
                    continue
                ptype = parts[1].upper()
                site = parts[2].strip().upper()
                if ptype not in ('STAX', 'STAY', 'STAZ'):
                    continue
                try:
                    val = float(parts[8])
                except Exception:
                    continue
                rec = coords_tmp.setdefault(site, {})
                if ptype == 'STAX':
                    rec['x'] = val
                elif ptype == 'STAY':
                    rec['y'] = val
                elif ptype == 'STAZ':
                    rec['z'] = val
        for site, rec in coords_tmp.items():
            if not all(k in rec for k in ('x', 'y', 'z')):
                continue
            x, y, z = rec['x'], rec['y'], rec['z']
            lat_deg, lon_deg, hgt_m = ecef_to_geodetic(x, y, z)
            tbl[site] = (x, y, z, lat_deg, lon_deg, hgt_m)
        if not tbl:
            raise ValueError("No station ECEF coordinates parsed from SINEX file: %s" % true_file)
        return tbl
    with open(true_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line.startswith('%'):
                continue
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            site = parts[0].strip().upper()
            try:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                lat = float(parts[4])
                lon = float(parts[5])
                hgt = float(parts[6])
            except Exception:
                continue
            tbl[site] = (x, y, z, lat, lon, hgt)
    return tbl

def compute_enu_error(ecef_arr, true_ecef, true_llh):
    e_ref = np.array(true_ecef, dtype=float)
    lat_deg, lon_deg, hgt = true_llh
    R = ecef_to_enu_matrix(lat_deg, lon_deg)
    diff = ecef_arr - e_ref
    enu = np.dot(R, diff)
    return enu

# 全局常量
CVG_CONTIOUS_EPNUM = 10
CVG_ERROR_THRES_3D = 0.1
CVG_ERROR_THRES_HRZ = 0.1
CVG_ERROR_THRES_HGT = 0.1


class denu_base_t:
    def __init__(self):
        self.tm = datetime.datetime.now()
        self.ns = [0, 0, 0, 0, 0]
        self.denu = [0.0, 0.0, 0.0]
        self.err_hrz = 0.0
        self.err_3d = 0.0


class error_t:
    def __init__(self):
        self.dir0 = ''
        self.stname = ''
        self.fname = ''
        self.ses_id = 0
        self.lst_denu_bases = []


class cvg_prc_opt_t:
    def __init__(self):
        self.interval_sec = 0
        self.is_sitebysite_output = False
        self.init_freq_min = 0
        self.dirres = ''
        self.start_time_str = ''
        self.end_time_str = ''


class cvg_res_1direction_t:
    def __init__(self):
        self.lst_site = []
        self.lst_ave_cvg = []
        self.nbad = 0
        self.nall = 0
        self.ave_cvg = 0.0


class pair_dir_file_mode_t:
    def __init__(self):
        self.dir0 = ''
        self.name = ''
        self.md_src = ''


def cal_cvg_3d_hrz(rEps):
    cvg_3d, cvg_hrz = 99999.9, 99999.9
    if len(rEps) <= 0:
        return cvg_3d, cvg_hrz

    lst_tms = []
    nzero, nep = 0, 0
    for r1 in rEps:
        nep += 1
        if r1.err_3d <= 0.0:
            nzero += 1
        else:
            nzero = 0

        if r1.err_3d <= CVG_ERROR_THRES_3D:
            lst_tms.append(r1.tm)
        else:
            lst_tms.clear()

        if len(lst_tms) >= CVG_CONTIOUS_EPNUM:
            if nzero >= nep - 1 and nzero >= CVG_CONTIOUS_EPNUM - 1 and nzero >= 2:
                lst_tms.clear()
                continue
            else:
                dtmp = (lst_tms[0] - rEps[0].tm).total_seconds()
                cvg_3d = (dtmp / 60.0)
                break

    lst_tms = []
    nzero, nep = 0, 0
    for r1 in rEps:
        nep += 1
        if r1.err_hrz <= 0.0:
            nzero += 1
        else:
            nzero = 0

        if r1.err_hrz <= CVG_ERROR_THRES_HRZ:
            lst_tms.append(r1.tm)
        else:
            lst_tms.clear()

        if len(lst_tms) >= CVG_CONTIOUS_EPNUM:
            if nzero >= nep - 1 and nzero >= CVG_CONTIOUS_EPNUM - 1 and nzero >= 2:
                lst_tms.clear()
                continue
            else:
                d1 = (lst_tms[0] - rEps[0].tm).total_seconds()
                cvg_hrz = (d1 / 60.0)
                break

    return cvg_3d, cvg_hrz


def cal_cvg_enu(rEps):
    cvg = [99999.9, 99999.9, 99999.9]
    if len(rEps) <= 0:
        return cvg

    err_thres = [CVG_ERROR_THRES_HRZ, CVG_ERROR_THRES_HRZ, CVG_ERROR_THRES_HGT]

    for i in range(3):
        lst_tms = []
        nzero, nep = 0, 0

        for r1 in rEps:
            nep += 1
            if abs(r1.denu[i]) <= 0.0:
                nzero += 1
            else:
                nzero = 0

            if abs(r1.denu[i]) <= err_thres[i]:
                lst_tms.append(r1.tm)
            else:
                lst_tms = []

            if len(lst_tms) >= CVG_CONTIOUS_EPNUM:
                if nzero >= nep - 1 and nzero >= CVG_CONTIOUS_EPNUM - 1 and nzero >= 2:
                    lst_tms = []
                    continue

                d1 = (lst_tms[0] - rEps[0].tm).total_seconds()
                cvg[i] = d1 / 60.0
                break

        if nzero >= len(rEps) - 2 and nzero >= 3:
            cvg[i] = 99999.9

    return cvg


def cal_ave_sat(rEps):
    ns_ave = [0.0 for _ in range(5)]
    if len(rEps) <= 0:
        return ns_ave

    for i in range(5):
        lst = [r1.ns[i] for r1 in rEps]
        ns_ave[i] = sum(lst) / len(lst)

    return ns_ave


def check_partion(nsite, nsub):
    npart = int(nsite / nsub)
    nfac = nsite - npart * nsub

    if nfac * 4 <= nsub:
        return True, npart, nfac, 0
    elif nfac * 4 >= 3 * nsub:
        return True, npart, nfac, -1
    return False, npart, nfac, 999


def divide_ep_ses_with_init_freq(err_1site: error_t, opt: cvg_prc_opt_t):
    if opt.init_freq_min <= 0:
        return [err_1site]

    lst_nodes, lst_seses = [], []

    r0 = err_1site.lst_denu_bases[0]
    lst_nodes.append(0)

    mi_ = r0.tm.hour * 60.0 + r0.tm.minute
    ses = int(mi_ / opt.init_freq_min + 0.5)
    lst_seses.append(ses)

    ix = 0
    for r in err_1site.lst_denu_bases[1:]:
        ix += 1
        mi_ = r.tm.hour * 60 + r.tm.minute

        is_new_ses1 = 0 == mi_ % opt.init_freq_min and r.tm.second < 0.01
        is_new_ses2 = (r.tm - r0.tm).total_seconds() >= opt.init_freq_min * 60.0

        if is_new_ses1 or is_new_ses2:
            r0 = r
            ses = int(mi_ / opt.init_freq_min + 0.5)
            lst_nodes.append(ix)
            lst_seses.append(ses)

    if ix > lst_nodes[-1] + 1:
        lst_nodes.append(ix)
        lst_seses.append(ses)

    nnodes = len(lst_nodes)

    nskip = 0
    nused = 0

    lst_res = []
    for i in range(nnodes - 1):
        ix0, ix1 = lst_nodes[i], lst_nodes[i + 1]
        if ix1 - ix0 < 3:
            continue

        err_ = error_t()
        err_.dir0 = err_1site.dir0
        err_.fname = err_1site.fname
        err_.stname = err_1site.stname
        err_.ses_id = lst_seses[i]
        err_.lst_denu_bases = err_1site.lst_denu_bases[ix0:ix1]

        if '-GREC23' in err_1site.fname:
            lst_ns_E, lst_ns_C = [], []
            for base in err_.lst_denu_bases:
                lst_ns_C.append(base.ns[-1])
                lst_ns_E.append(base.ns[2])
            nsmax_E, nsmax_C = max(lst_ns_E), max(lst_ns_C)
            if nsmax_E <= 2 or nsmax_C <= 4:
                print('skip [%s]-[%02d] due to no E/C sat used [%2d %2d]' % (err_.fname, err_.ses_id, nsmax_E, nsmax_C))
                nskip += 1
                continue
        nused += 1

        lst_res.append(err_)
    if nskip > 0:
        print('nused/nskip = %d/%d' % (nused, nskip))
    return lst_res


def get_modestr_from_folder(dir0):
    substop1 = 'result-'
    substop2 = 'result_'
    id00 = dir0.find(substop1, 1)
    id11 = dir0.find(substop2, 1)
    id0, substop = -1, 'xxx'
    if id00 >= 0:
        id0, substop = id00, substop1
    elif id11 >= 0:
        id0, substop = id11, substop2

    if id0 < 0:
        return ''
    md0 = dir0[id0 + len(substop):].strip()
    return md0


def get_modestr_from_filename(nm):
    sub = '-PPP_'
    id0 = nm.find(sub)
    if id0 >= 0:
        md1 = nm[id0 + 1:].replace('.ppprtk', '').replace('.pos', '').replace('.flt', '')
        return md1
    # For .flt files with suffixes like _float, _pppar, _ppprtk
    if nm.endswith('.flt'):
        base = nm[:-4]
        for suffix in ['_float', '_pppar', '_ppprtk']:
            if base.endswith(suffix):
                return suffix[1:]
        return 'flt'
    print('no substr [%s] in [%s]' % (sub, nm))
    input('not a valid file ' + nm)


def simplify_mode_str_list_remove_cmn_str(lst_modes_all):
    if len(lst_modes_all) <= 1:
        return lst_modes_all

    lst_len = [len(s0) for s0 in lst_modes_all]
    len0 = min(lst_len)

    ipos0, ipos9 = 0, -999

    for i in range(0, len0):
        lst_ch = [s0[i] for s0 in lst_modes_all]
        lst_ch = list(set(lst_ch))
        if 1 == len(lst_ch):
            ipos0 = i
            continue
        else:
            break

    for i in range(-1, -len0, -1):
        lst_ch = [s0[i] for s0 in lst_modes_all]
        lst_ch = list(set(lst_ch))
        if 1 == len(lst_ch):
            ipos9 = i
            continue
        else:
            break

    if -999 == ipos9:
        lst_new = [s0[ipos0 + 1:] for s0 in lst_modes_all]
    else:
        lst_new = [s0[ipos0 + 1:ipos9] for s0 in lst_modes_all]

    return lst_new


def rename_mode_str(md_str):
    md = md_str
    return md


def get_pos_file_name_list(dir0):
    lst_nms = []

    lst_tmp = os.listdir(dir0)
    lst_tmp.sort()
    for nm in lst_tmp:
        if (not nm.endswith('.ppprtk') and not nm.endswith('.pos') and not nm.endswith('.flt')) or '-augoff.' in nm:
            continue

        fpath_f = os.path.join(dir0, nm)
        if not os.path.exists(fpath_f):
            print('path [%s] not exists.' % fpath_f)
            continue

        sz_kb = os.path.getsize(fpath_f) / 1024.0
        if sz_kb < 5:
            print('[%s] is too short. only [%4.1f]kb.' % (fpath_f, sz_kb))
            continue
        lst_nms.append(nm)

    return lst_nms


def get_percentile_x_y_series(perct, dicts_lst_cvg):
    lst_x, lst_y, lst_nsmp = [], [], []

    for tm, lst in dicts_lst_cvg.items():
        if len(lst) < 5:
            continue

        lst.sort()
        val_90 = np.percentile(lst, perct)

        lst_x.append(tm / 60.0)
        lst_y.append(val_90)
        lst_nsmp.append(len(lst))

    return [lst_x, lst_y, lst_nsmp]


def draw_cvg_series_allmodes_1direction(lst_x, lst_y, lst_md, perct, str_direction, pth_jpg):
    """单独绘制每个方向的图（保留原功能）"""
    plt.figure(figsize=(129*mm,100*mm))
    
    x9, xmax_set = -999.0, 90

    lst_clr = ['C0', 'C1', 'C2', 'C4', 'C5', 'C3', 'C6', 'C7', 'C8', 'C9']
    for x, y, md, clr in zip(lst_x, lst_y, lst_md, lst_clr):
        if not x:
            continue
        x9 = max(x9, max(x))
        plt.plot(x, y, label=md, color=clr, linewidth=1.5)

    if x9 < 0:
        return
    
    xmax_set = 30 if x9 <= 30 else 45 if x9 <= 45 else 60 if x9 <= 60 else 90 if x9 <= 90 else x9

    plt.legend(loc='upper right', fontsize=8, prop=font_prop)

    plt.ylim(ymin=0, ymax=0.401)
    plt.xlim(0, 60)
    plt.ylabel('%s Error(m)' % str_direction, fontsize=9, fontproperties=font_prop)
    plt.xlabel('Time(min)', fontsize=9, fontproperties=font_prop)
    plt.title('%d-Percentile' % perct, fontsize=10, fontproperties=font_prop)

    plt.hlines(y=0.2, xmin=0, xmax=x9, color='gray', linestyles=':', linewidth=0.8)
    plt.hlines(y=0.1, xmin=0, xmax=x9, color='gray', linestyles=':', linewidth=0.8)
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig(pth_jpg, dpi=400, bbox_inches='tight')
    plt.close('all')


def draw_cvg_series_allmodes_1direction_combined(lst_lst_x_dict, lst_lst_y_dict, lst_md, perct, str_direction, pth_jpg):
    """将同一百分位数下的East/North/Vertical/3D绘制在一张图上"""
    plt.figure(figsize=(129*mm,100*mm))
    
    # 定义子图布局
    directions = ['East', 'North', 'Vertical', '3D']
    colors_list = ['C0', 'C1', 'C2', 'C4', 'C5', 'C3', 'C6', 'C7', 'C8', 'C9']  
    
    for idx, (dir_name) in enumerate(directions):
        ax = plt.subplot(2, 2, idx + 1)
        
        if dir_name in lst_lst_x_dict:
            lst_x = lst_lst_x_dict[dir_name]
            lst_y = lst_lst_y_dict[dir_name]
            
            x9 = -999.0
            for x, y, md, color in zip(lst_x, lst_y, lst_md, colors_list):
                if not x:
                    continue
                x9 = max(x9, max(x))
                # 使用该子图对应的颜色绘制所有策略
                ax.plot(x, y, label=md, color=color, linewidth=1.5)
            
            if x9 < 0:
                continue
            
            # 只在第一个子图显示legend
            if idx == 1:
                ax.legend(loc='upper right', fontsize=7)
            
            ax.set_ylim(ymin=0, ymax=0.401)
            ax.set_xlim(0, 60)
            
            # 纵坐标标题：只在子图1和3显示
            if idx == 0 or idx == 2:
                ax.set_ylabel('Error (m)', fontsize=8)
            
            # 横坐标标题：只在子图3和4显示
            if idx == 2 or idx == 3:
                ax.set_xlabel('Time (min)', fontsize=8)
            
            ax.set_title('%s - %d-Percentile' % (dir_name, perct), fontsize=8)
            ax.hlines(y=0.2, xmin=0, xmax=x9, color='gray', linestyles=':', linewidth=0.8)
            ax.hlines(y=0.1, xmin=0, xmax=x9, color='gray', linestyles=':', linewidth=0.8)
            ax.grid(True, linestyle=':', alpha=0.6)
            
            # 设置刻度标签字体
            for label in ax.get_xticklabels():
                label.set_fontproperties(font_prop)
                label.set_fontsize(8)
            for label in ax.get_yticklabels():
                label.set_fontproperties(font_prop)
                label.set_fontsize(8)
    
    # 调整子图间距，减小间隔
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    plt.tight_layout()
    plt.savefig(pth_jpg, dpi=400, bbox_inches='tight')
    plt.close('all')


def read_err_series_sub(pth_series, lines, init_freq_min, target_site, is_bad_skip):
    is_good = True
    is_tgt = len(target_site) >= 2

    nm_site = ''
    lst_tm_dif, lst_err_hrz, lst_err_e, lst_err_n, lst_err_up, lst_err_3d = [], [], [], [], [], []

    for ln in lines:
        if len(ln) < 5:
            continue
        if ln.startswith('##'):
            if '##%' in ln and '[folder]' in ln:
                continue

            if '[folder]' in ln:
                input('not a good line. should be skip: ' + ln)
            nm_file = ln.split()[1]
            nm_site = nm_file[0:4]
            is_good = True
            if 'bad' == ln.split()[-1]:
                is_good = False
            continue
        if is_tgt:
            if nm_site.upper() != target_site.upper():
                continue
        if is_bad_skip:
            if not is_good:
                continue

        ss = ln.split()

        tm_dif, e, n, err_up, err_hrz, err_3d = int(float(ss[0])), abs(float(ss[1])), abs(float(ss[2])), abs(
            float(ss[-3])), float(ss[-2]), float(ss[-1])
        lst_tm_dif.append(tm_dif)
        lst_err_e.append(e)
        lst_err_n.append(n)
        lst_err_up.append(err_up)
        lst_err_hrz.append(err_hrz)
        lst_err_3d.append(err_3d)

        if tm_dif > init_freq_min * 60.0:
            print('[%s] tm_dif=%.1f > init_freq_min=%.1f' % (pth_series, tm_dif / 60.0, init_freq_min))
            print('line : ' + ln.strip())
            continue

        if tm_dif > 15400:
            print('skip [%s] because tm_dif:%.1f > 15400' % (ln, tm_dif))
            continue
    return lst_tm_dif, lst_err_hrz, lst_err_e, lst_err_n, lst_err_up, lst_err_3d


def read_err_series(pth_cvg, is_bad_skip, init_freq_min=10000, target_site=''):
    with open(pth_cvg) as fp:
        lines = fp.readlines()

        lst_tm_dif, lst_err_hrz, lst_err_e, lst_err_n, lst_err_up, lst_err_3d = read_err_series_sub(pth_cvg, lines,
                                                                                                    init_freq_min,
                                                                                                    target_site,
                                                                                                    is_bad_skip)

        dicts_lst_errhrz, dicts_lst_errE, dicts_lst_errN, dicts_lst_errup, dicts_lst_err3d = dict(), dict(), dict(), dict(), dict()
        lst_tm_dif_no_rpt = list(set(lst_tm_dif))
        lst_tm_dif_no_rpt.sort()

        for tmdif in lst_tm_dif_no_rpt:
            dicts_lst_errE[tmdif], dicts_lst_errN[tmdif], dicts_lst_errhrz[tmdif], dicts_lst_errup[tmdif], dicts_lst_err3d[tmdif] = [], [], [], [], []

        for tmdif, errE, errN, err_hrz, err_up, err_3d in zip(lst_tm_dif, lst_err_e, lst_err_n, lst_err_hrz,
                                                              lst_err_up, lst_err_3d):
            dicts_lst_errE[tmdif].append(errE)
            dicts_lst_errN[tmdif].append(errN)
            dicts_lst_errhrz[tmdif].append(err_hrz)
            dicts_lst_errup[tmdif].append(err_up)
            dicts_lst_err3d[tmdif].append(err_3d)

        return dicts_lst_errE, dicts_lst_errN, dicts_lst_errhrz, dicts_lst_errup, dicts_lst_err3d


def get_doy_from_file_name(nm: str) -> int:
    """从文件名中提取年积日(DOY),文件名第13-19位为年份+年积日(例如2023070)"""
    if len(nm) >= 19:
        try:
            year_doy_str = nm[12:19]
            if len(year_doy_str) == 7 and year_doy_str.isdigit():
                doy = int(year_doy_str[1:7])
                return doy
        except:
            pass
    return 0


def read_draw_cvg_series_multi_files(lst_pth_cvg, lst_mds, dirout, dir_cvg_pert, dir_cvg_series_each_site, is_bad_skip,
                                     init_freq_min=0, target_site='', perct=68):
    print('call read_draw_cvg_series_multi_files perct=%d' % perct)

    doy_str = ''
    if len(target_site) > 2 and len(lst_pth_cvg) > 0:
        pth_cvg_time = lst_pth_cvg[0]
        pth_cvg_series = pth_cvg_time.replace('cvg_time-mode-', 'cvg_series-mode-')
        if os.path.exists(pth_cvg_series):
            try:
                with open(pth_cvg_series, 'r') as fp:
                    for ln in fp:
                        if ln.startswith('###') and target_site.upper() in ln.upper():
                            parts = ln.split()
                            if len(parts) >= 2:
                                filename = parts[1]
                                doy = get_doy_from_file_name(filename)
                                if doy > 0:
                                    doy_str = '-DOY%03d' % doy
                                    break
            except:
                pass

    lst_md_new = []
    # 存储各个方向的数据，用于合并绘图
    dict_lst_x = {'East': [], 'North': [], 'Vertical': [], '3D': []}
    dict_lst_y = {'East': [], 'North': [], 'Vertical': [], '3D': []}
    
    # 存储各个方向的数据，用于单独绘图
    dict_single_x = {'East': [], 'North': [], 'Horizontal': [], 'Vertical': [], '3D': []}
    dict_single_y = {'East': [], 'North': [], 'Horizontal': [], 'Vertical': [], '3D': []}

    for pth_cvg, md in zip(lst_pth_cvg, lst_mds):
        md = rename_mode_str(md)
        lst_md_new.append(md)

        dicts_lst_errE, dicts_lst_errN, dicts_lst_errhrz, dicts_lst_errup, dicts_lst_err3d = read_err_series(pth_cvg,is_bad_skip,init_freq_min,target_site)

        lst_tmp_dict = [dicts_lst_errE, dicts_lst_errN, dicts_lst_errhrz, dicts_lst_errup, dicts_lst_err3d]
        lst_x_, lst_y_, lst_nsmp_ = [], [], []

        for tmp in lst_tmp_dict:
            x, y, n = get_percentile_x_y_series(perct, tmp)
            lst_x_.append(x)
            lst_y_.append(y)
            lst_nsmp_.append(n)

        lst_nsmp = lst_nsmp_[0]

        # 存储各方向的数据用于合并绘图
        dict_lst_x['East'].append(lst_x_[0])
        dict_lst_x['North'].append(lst_x_[1])
        dict_lst_x['Vertical'].append(lst_x_[3])
        dict_lst_x['3D'].append(lst_x_[4])

        dict_lst_y['East'].append(lst_y_[0])
        dict_lst_y['North'].append(lst_y_[1])
        dict_lst_y['Vertical'].append(lst_y_[3])
        dict_lst_y['3D'].append(lst_y_[4])
        
        # 存储各方向的数据用于单独绘图
        dict_single_x['East'].append(lst_x_[0])
        dict_single_x['North'].append(lst_x_[1])
        dict_single_x['Horizontal'].append(lst_x_[2])
        dict_single_x['Vertical'].append(lst_x_[3])
        dict_single_x['3D'].append(lst_x_[4])
        
        dict_single_y['East'].append(lst_y_[0])
        dict_single_y['North'].append(lst_y_[1])
        dict_single_y['Horizontal'].append(lst_y_[2])
        dict_single_y['Vertical'].append(lst_y_[3])
        dict_single_y['3D'].append(lst_y_[4])

        if len(lst_x_[0]) == len(lst_x_[1]) and len(lst_x_[1]) == len(lst_x_[2]) and len(lst_x_[1]) == len(
                lst_x_[3]) and len(lst_x_[1]) == len(lst_x_[4]):

            nm0 = os.path.split(pth_cvg)[-1]
            nm0 = os.path.splitext(nm0)[0]

            if len(target_site) > 2:
                dir_tmp = dir_cvg_series_each_site
                nm0 = nm0 + '-' + target_site
            else:
                dir_tmp = dir_cvg_pert
            pth_cvg_perct = os.path.join(dir_tmp, '%s-P%.0f.perct' % (nm0, perct))

            with open(pth_cvg_perct, 'w') as fp:
                fp.write('%7s  %9s  %9s  %9s  %9s  %9s  %9s \n' % ('time', 'E', 'N', 'HRZ', 'U', '3D', 'nsmp'))
                for th, ve, vn, vh, vu, v3, nsmp in zip(lst_x_[0], lst_y_[0], lst_y_[1], lst_y_[2], lst_y_[3],
                                                        lst_y_[4], lst_nsmp):
                    fp.write('%7.1f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9d\n' % (th, ve, vn, vh, vu, v3, nsmp))

    # 绘制合并图
    if 0 < len(lst_pth_cvg):
        dir_tmp = dirout
        if len(target_site) > 2:
            dir_tmp = dir_cvg_series_each_site
            nm_combined = 'allmodes-%s%s-P%d-combined.png' % (target_site, doy_str, perct)
        else:
            nm_combined = 'allmodes-P%d-combined.png' % perct

        if not os.path.exists(dir_tmp):
            print('makedirs [%s] for draw_cvg_series_' % dir_tmp)
            os.makedirs(dir_tmp)

        pth_jpg = os.path.join(dir_tmp, nm_combined)
        draw_cvg_series_allmodes_1direction_combined(dict_lst_x, dict_lst_y, lst_md_new, perct, 'combined', pth_jpg)
        
        # 绘制单独各方向的图
        str_directions = ['East', 'North', 'Horizontal', 'Vertical', '3D']
        for str_dir in str_directions:
            if len(target_site) > 2:
                nm_single = 'allmodes-%s%s-P%d-%s.png' % (target_site, doy_str, perct, str_dir.lower())
            else:
                nm_single = 'allmodes-P%d-%s.png' % (perct, str_dir.lower())
            
            pth_single_jpg = os.path.join(dir_tmp, nm_single)
            draw_cvg_series_allmodes_1direction(dict_single_x[str_dir], dict_single_y[str_dir], lst_md_new, perct, str_dir, pth_single_jpg )


def rm_mk_dir(dir0):
    if os.path.exists(dir0):
        shutil.rmtree(dir0)
    os.makedirs(dir0)


def check_dirlst(dirlst):
    for dir0 in dirlst:
        if not os.path.exists(dir0):
            print('[check_dirlst] dir [%s] is not exists. skip' % (dir0))
            return False
        elif not os.path.isdir(dir0):
            print('[check_dirlst] dir [%s] is not dir. skip' % (dir0))
            return False
    return True


def replace_special_mode_substring(md1):
    md1 = md1.replace('0=Flo', 'Flo').replace('1=AR', 'AR')
    for t in range(2, 100):
        md1 = md1.replace('%d=' % t, '')
    return md1


def get_dir_file_pair_md_list(dirlst, lst_md_input):
    lst_dir_name_md = []

    for idx, dir0 in enumerate(dirlst):
        dir0 = dir0.rstrip('/')
        dir0 = dir0.rstrip('\\')

        lst_nms = get_pos_file_name_list(dir0)
        print('[%3d] files in [%s]' % (len(lst_nms), dir0))

        md0 = get_modestr_from_folder(dir0)

        for nm in lst_nms:
            if lst_md_input:
                md = lst_md_input[idx]
            else:
                md1 = get_modestr_from_filename(nm)
                md = md0 + '-' + md1

            inf = pair_dir_file_mode_t()
            inf.dir0 = dir0
            inf.name = nm
            inf.md_src = md

            lst_dir_name_md.append(inf)

    return lst_dir_name_md


def divide_sub_sites(nsite, lst_site, lst_cvg):
    lst_lst_site, lst_lst_cvg = [], []

    for nsub in [20, 24, 25, 16, 15, 28, 30, 100, 200, 150, 99999999, nsite]:
        is_gd, npart, nfac, istat = check_partion(nsite, nsub)

        if not is_gd:
            continue

        ix9 = 0
        for i in range(npart):
            ix0 = i * nsub
            ix9 = ix0 + nsub

            lst_lst_site.append(lst_site[ix0:ix9])
            lst_lst_cvg.append(lst_cvg[ix0:ix9])

        if 0 == istat:
            if npart > 0:
                lst_lst_site[-1].extend(lst_site[ix9:])
                lst_lst_cvg[-1].extend(lst_cvg[ix9:])
            else:
                lst_lst_site.append(lst_site[ix9:])
                lst_lst_cvg.append(lst_cvg[ix9:])
        elif -1 == istat:
            lst_lst_site.append(lst_site[ix9:])
            lst_lst_cvg.append(lst_cvg[ix9:])
        break

    return lst_lst_site, lst_lst_cvg


def draw_cvg_hbar_sub(str_direction, cvg_ave_all, tup_lst_site_cvg_pth):
    lst_site_, lst_cvg_, pth = tup_lst_site_cvg_pth

    lst_site = [s for s in lst_site_]
    lst_cvg = [c for c in lst_cvg_]
    lst_site.reverse()
    lst_cvg.reverse()

    fig, ax = plt.subplots(figsize=(129/25.4, 80.625/25.4), dpi=600)
    ax_rect = ax.barh(lst_site, lst_cvg)

    def autolabel(rects, xlimmax):
        for rect in rects:
            width = rect.get_width()
            s0 = '%.1f' % width
            if width > xlimmax:
                width = xlimmax
            ax.annotate(s0,
                        xy=(width, rect.get_y() + rect.get_height() * 0 / 3),
                        xytext=(15, 0),
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontproperties=font_prop)

    max_x, max_y = 0.0, 0.0
    for rect in ax_rect:
        width = rect.get_width()
        if width > max_x:
            max_x = width

        tmp = rect.get_y()
        if tmp > max_y:
            max_y = tmp

    fig.tight_layout()

    cvg9 = max(lst_cvg)
    x9 = 40 if cvg9 <= 40 else 60
    cvg_ave = sum(lst_cvg) / len(lst_cvg)

    autolabel(ax_rect, x9)

    s0 = 'AVE    : %.1f' % cvg_ave
    if cvg_ave_all > 0.0:
        s0 += '\nAVE_ALL: %.1f' % cvg_ave_all
    ax.text(x9 * 3 / 4, max_y * 4 / 7, s0, color='r', fontproperties=font_prop)

    ax.set_xlim(xmax=x9, xmin=0)
    ax.set_ylim(ymin=-1, ymax=len(lst_site))
    ax.set_title('CVG-%s/min' % str_direction.upper(), fontproperties=font_prop)
    ax.set_xlabel('Time (min)', fontproperties=font_prop)
    ax.set_ylabel('Site', fontproperties=font_prop)
    for label in ax.get_yticklabels():
        label.set_fontproperties(font_prop)
    for label in ax.get_xticklabels():
        label.set_fontproperties(font_prop)
    fig.tight_layout()
    plt.savefig(pth, bbox_inches='tight', dpi=600)
    plt.close('all')


def draw_cvg_hbar(res_cvg_sum: cvg_res_1direction_t, str_mode, str_direction, dirres):
    nsub = 30
    nsite = len(res_cvg_sum.lst_site)

    if nsite <= 0:
        print('[%s-%s] no valid site for hbar. skip' % (str_mode, str_direction))
        return

    if nsite <= nsub:
        lst_lst_site = [res_cvg_sum.lst_site]
        lst_lst_cvg = [res_cvg_sum.lst_ave_cvg]
    else:
        lst_lst_site, lst_lst_cvg = divide_sub_sites(nsite, res_cvg_sum.lst_site, res_cvg_sum.lst_ave_cvg)

    lst_tup_info = []
    n = 0
    for lst_s, lst_c in zip(lst_lst_site, lst_lst_cvg):
        n += 1
        pth0 = os.path.join(dirres, 'site-cvg-bar-%s-%s-sub%d.png' % (str_mode, str_direction, n))
        lst_tup_info.append((lst_s, lst_c, pth0))

    pt = partial(draw_cvg_hbar_sub, str_direction, res_cvg_sum.ave_cvg)

    with multiprocessing.Pool() as p:
        p.map(pt, lst_tup_info)


def check_is_pos_line(ln0):
    if ('/' in ln0 and '.' in ln0 and ln0.startswith(('19', '20'))):
        return True
    else:
        return False


def add_denu_base_to_list(ep, err_1site):
    ep.err_3d = math.sqrt(ep.denu[0] * ep.denu[0] + ep.denu[1] * ep.denu[1] + ep.denu[2] * ep.denu[2])
    ep.err_hrz = math.sqrt(ep.denu[0] * ep.denu[0] + ep.denu[1] * ep.denu[1])
    err_1site.lst_denu_bases.append(ep)


def parse_time_bound(tm_str, ref_dt=None):
    """解析时间边界字符串。
    支持完整日期时间格式: YYYY/MM/DD HH:MM:SS 或 YYYY/MM/DD HH:MM:SS.sss
    支持仅时间格式: HH:MM:SS (自动使用参考日期 ref_dt 补全)
    """
    if not tm_str or len(tm_str.strip()) < 3:
        return None
    tm_str = tm_str.strip()
    if '/' in tm_str:
        # 完整日期时间
        fmt = '%Y/%m/%d %H:%M:%S'
        if '.' in tm_str.split()[-1]:
            fmt = '%Y/%m/%d %H:%M:%S.%f'
        return datetime.datetime.strptime(tm_str, fmt)
    else:
        # 仅时间，需要参考日期
        if ref_dt is None:
            return None
        tm = datetime.datetime.strptime(tm_str, '%H:%M:%S').time()
        return datetime.datetime.combine(ref_dt.date(), tm)


def read_1_posfile(fpath: str, interval_sec: int = 0, start_time_str: str = '', end_time_str: str = '') -> error_t:
    err_1site = error_t()

    nep, nzero = 0, 0
    start_time = None
    end_time = None
    bounds_parsed = not (start_time_str or end_time_str)

    with open(fpath) as fp:
        lines = fp.readlines()

        for ln in lines:
            if 160 > len(ln):
                continue
            if not (':' in ln and '/' in ln):
                continue
            if ln.startswith('% '):
                continue
            ep = denu_base_t()
            ss = ln.split()
            if 18 > len(ss):
                continue
            try:
                ep.tm = datetime.datetime.strptime(ln[0:23].strip(), '%Y/%m/%d %H:%M:%S.%f')
            except ValueError:
                continue

            # 利用首个历元日期解析仅含时间的边界字符串
            if not bounds_parsed:
                start_time = parse_time_bound(start_time_str, ep.tm)
                end_time = parse_time_bound(end_time_str, ep.tm)
                bounds_parsed = True

            # 应用起止时间过滤
            if start_time and ep.tm < start_time:
                continue
            if end_time and ep.tm > end_time:
                continue

            if interval_sec > 0:
                sec = ep.tm.second + ep.tm.minute * 60
                if 0 != sec % interval_sec:
                    continue

            if fpath.endswith('.pos'):
                if len(ss) < 20:
                    continue
                for i in range(0, 3):
                    ep.denu[i] = float(ss[17 + i])
            else:
                if '0=Flo' in fpath:
                    for i in range(0, 3):
                        ep.denu[i] = float(ss[2 + i])
                else:
                    for i in range(0, 3):
                        ep.denu[i] = float(ss[i + 5])
            ep.err_3d = math.sqrt(ep.denu[0] * ep.denu[0] + ep.denu[1] * ep.denu[1] + ep.denu[2] * ep.denu[2])
            ep.err_hrz = math.sqrt(ep.denu[0] * ep.denu[0] + ep.denu[1] * ep.denu[1])
            for i in range(0, 5):
                ep.ns[i] = int(ss[8 + i])

            err_1site.lst_denu_bases.append(ep)

            nep += 1
            if ep.err_3d <= 0.000001:
                nzero += 1

    if nep - nzero <= 3:
        err_1site.lst_denu_bases.clear()
        print('file [%s] zero enu error' % fpath)

    return err_1site


def read_1_fltfile(fpath, true_tbl, interval_sec=0, start_time_str='', end_time_str=''):
    """Read GREAT-PVT .flt file (ECEF XYZ) and compute ENU errors using true coordinates."""
    err_1site = error_t()

    basename = os.path.basename(fpath)
    stname = basename[0:4].upper()

    tcoord = true_tbl.get(stname)
    if tcoord is None:
        name_noext = os.path.splitext(basename)[0]
        tokens = re.split(r"[^A-Za-z0-9]+", name_noext)
        tokens = [t.upper() for t in tokens if t]
        for t in tokens:
            if t in true_tbl:
                tcoord = true_tbl[t]
                stname = t
                break
    if tcoord is None:
        print('WARNING: cannot find true coordinates for site [%s] in true table. skip [%s]' % (stname, fpath))
        return err_1site

    true_ecef = (tcoord[0], tcoord[1], tcoord[2])
    true_llh = (tcoord[3], tcoord[4], tcoord[5])

    # Parse time bounds (seconds of day) if provided as HH:MM:SS
    start_sec = -1.0
    end_sec = -1.0
    if start_time_str:
        try:
            t = datetime.datetime.strptime(start_time_str.strip(), '%H:%M:%S')
            start_sec = t.hour * 3600 + t.minute * 60 + t.second
        except Exception:
            pass
    if end_time_str:
        try:
            t = datetime.datetime.strptime(end_time_str.strip(), '%H:%M:%S')
            end_sec = t.hour * 3600 + t.minute * 60 + t.second
        except Exception:
            pass

    nep, nzero = 0, 0
    with open(fpath) as fp:
        for ln in fp:
            if not ln.strip() or ln.strip().startswith('#'):
                continue
            ss = ln.split()
            if len(ss) < 19:
                continue
            try:
                sow = float(ss[0])
                x = float(ss[1])
                y = float(ss[2])
                z = float(ss[3])
                nsat = int(ss[13])
            except ValueError:
                continue

            if interval_sec > 0:
                if int(sow) % interval_sec != 0:
                    continue

            # Apply time bounds in seconds-of-week space if parsed
            if start_sec >= 0 and sow < start_sec:
                continue
            if end_sec >= 0 and sow > end_sec:
                continue

            ecef_arr = np.array([x, y, z], dtype=float)
            enu = compute_enu_error(ecef_arr, true_ecef, true_llh)

            ep = denu_base_t()
            # Use a fixed base date + SOW seconds for relative time difference calculations
            ep.tm = datetime.datetime(2020, 1, 1) + datetime.timedelta(seconds=sow)
            ep.denu[0] = float(enu[0])   # East
            ep.denu[1] = float(enu[1])   # North
            ep.denu[2] = float(enu[2])   # Up
            ep.err_3d = math.sqrt(enu[0]**2 + enu[1]**2 + enu[2]**2)
            ep.err_hrz = math.sqrt(enu[0]**2 + enu[1]**2)
            ep.ns[0] = nsat
            for i in range(1, 5):
                ep.ns[i] = 0

            err_1site.lst_denu_bases.append(ep)
            nep += 1
            if ep.err_3d <= 0.000001:
                nzero += 1

    if nep - nzero <= 3:
        err_1site.lst_denu_bases.clear()
        print('file [%s] zero enu error' % fpath)

    return err_1site


def write_ave_cvg_sum_file(pth_cvg_ave_sum, lst_cvg_sum_1dir: List[cvg_res_1direction_t]):
    lst_all_site = []
    for i in range(5):
        lst_all_site.extend(lst_cvg_sum_1dir[i].lst_site)
    lst_all_site = list(set(lst_all_site))
    lst_all_site.sort()

    dict_site_cvg = dict()
    for st in lst_all_site:
        dict_site_cvg[st] = [999.9 for _ in range(5)]
    for i in range(5):
        for st, cvg in zip(lst_cvg_sum_1dir[i].lst_site, lst_cvg_sum_1dir[i].lst_ave_cvg):
            dict_site_cvg[st][i] = cvg

    with open(pth_cvg_ave_sum, 'w') as fp:
        fp.write('%10s %9s %9s %9s %9s %9s %9s %9s\n' % ('', 'e', 'n', 'u', 'hrz', '3d', 'nall', 'nbad='))

        s0 = '  AVE-CVG '
        str_ave_cvg = ''
        for i in range(5):
            str_ave_cvg += ' %9.2f' % lst_cvg_sum_1dir[i].ave_cvg
        str_ave_cvg += ' %9d' % (lst_cvg_sum_1dir[-1].nall)
        for i in range(5):
            str_ave_cvg += ' %4d' % lst_cvg_sum_1dir[i].nbad
        s0 += str_ave_cvg

        fp.write(s0 + '\n')

        for st, lst_cvg in dict_site_cvg.items():
            s0 = '  %7s ' % st
            for i in range(5):
                if lst_cvg[i] > 998.0:
                    s0 += '   ---    '
                else:
                    s0 += ' %9.2f' % lst_cvg[i]

            fp.write(s0 + '\n')
    return str_ave_cvg


def write_pos_session_res(lst_res, fp_cvg, fp_series):
    for err in lst_res:
        cvg_enu = cal_cvg_enu(err.lst_denu_bases)
        cvg_3d, cvg_hrz = cal_cvg_3d_hrz(err.lst_denu_bases)
        ave_sat = cal_ave_sat(err.lst_denu_bases)

        str_good = 'g'
        if cvg_3d > 999:
            d3d_max = 0.0
            for base in err.lst_denu_bases[-CVG_CONTIOUS_EPNUM:]:
                if abs(base.denu[-1]) > d3d_max:
                    d3d_max = abs(base.denu[-1])
            if d3d_max > 0.2 and 1:
                str_good = 'bad'

        s0 = '%45s  %4d  %5d   %5.2f %5.2f %5.2f %5.2f %5.2f    %7.1f %7.1f %7.1f  %7.1f %7.1f %s' % (
            err.fname, err.ses_id, len(err.lst_denu_bases), ave_sat[0], ave_sat[1], ave_sat[2], ave_sat[3],
            ave_sat[4],
            cvg_enu[0], cvg_enu[1], cvg_enu[2], cvg_hrz, cvg_3d, str_good)
        fp_cvg.write(s0 + '\n')

        fp_series.write('### %40s  %3d  %s\n' % (err.fname, err.ses_id, str_good))
        for r in err.lst_denu_bases:
            s0 = ' %7.1f   %7.3f %7.3f %7.3f   %7.3f %7.3f\n' % (
                (r.tm - err.lst_denu_bases[0].tm).total_seconds(), r.denu[0], r.denu[1], r.denu[2], r.err_hrz,
                r.err_3d)
            fp_series.write(s0)

    fp_cvg.flush()


def get_title_cvg():
    title_cvg = '%45s %5s  %5s   %5s %5s %5s %5s %5s    %7s %7s %7s  %7s %7s / min  \n' % (
    "file       ", "ses", "nep", "ave_G", 'ave_R', 'ave_E', 'ave_C', 'ave_J',
    'cvg-E', 'cvg-N', 'cvg-U', 'cvg_H', 'cvg-3D')
    return title_cvg


def sum_cvg_one_mode_one_direction(str_pth, ix_pos, first_conv_period):
    lst = []
    with open(str_pth) as fp:
        lines = fp.readlines()
        for ln in lines[1:]:
            if ln.startswith("##"):
                continue
            columns = ln.strip().split()
            if first_conv_period and int(columns[1]) == 0:
                continue
            ss = ln.split()
            lst.append((ss[0][0:4], float(ss[ix_pos])))

    lst.sort(key=lambda x: x[0])

    res_cvg_sum = cvg_res_1direction_t()
    res_cvg_sum.nbad = 0
    res_cvg_sum.nall = len(lst)

    for site, items in groupby(lst, key=lambda x: x[0]):
        items = list(items)
        tmp1 = []
        for i in items:
            if i[1] > 200.0:
                res_cvg_sum.nbad += 1
            else:
                tmp1.append(i[1])

        if 1 <= len(tmp1):
            cvg = sum(tmp1) / len(tmp1)
        else:
            continue

        res_cvg_sum.lst_ave_cvg.append(cvg)
        res_cvg_sum.lst_site.append(site)

    nsite = len(res_cvg_sum.lst_site)
    res_cvg_sum.ave_cvg = -999.999
    if nsite > 0:
        res_cvg_sum.ave_cvg = sum(res_cvg_sum.lst_ave_cvg) / nsite

    return res_cvg_sum


def sum_write_cvg_all_modes_multi_files(lst_modes, lst_fp_path, dirres, first_conv_period):
    print('call sum_write_cvg_all_modes_multi_files')
    lst_pos = [8, 9, 10, 11, 12]
    lst_str_dir = ['e', 'n', 'up', 'hrz', '3D']

    lst_ave_sum_str = []
    for md, pth in zip(lst_modes, lst_fp_path):
        lst_cvg_sum_1dir = []

        for ix_pos, str_direction in zip(lst_pos, lst_str_dir):
            cvg_sum_1dir = sum_cvg_one_mode_one_direction(pth, ix_pos, first_conv_period)
            draw_cvg_hbar(cvg_sum_1dir, md, str_direction, dirres)
            lst_cvg_sum_1dir.append(cvg_sum_1dir)

        pth_cvg_ave_sum = os.path.join(dirres, 'ave-cvg-sum-%s.txt' % md)
        s0 = write_ave_cvg_sum_file(pth_cvg_ave_sum, lst_cvg_sum_1dir)
        lst_ave_sum_str.append([md, s0])

    pth = os.path.join(dirres, 'ALL-modes-ave-cvg-sum.txt')
    with open(pth, 'w') as fp:
        lst_len = [len(s[0]) for s in lst_ave_sum_str]
        lenmax = max(lst_len) + 4

        fp.write(
            '%s %9s %9s %9s %9s %9s %9s %9s\n' % ('mode'.center(lenmax), 'e', 'n', 'u', 'hrz', '3d', 'nall', 'nbad='))

        for s in lst_ave_sum_str:
            fp.write('%s %s\n' % (s[0].center(lenmax), s[1]))


def get_pos_result_ses_by_ses_1file(opt: cvg_prc_opt_t, true_tbl, tup_dir_nm):
    fdir, fnm = tup_dir_nm
    fpath_f = os.path.join(fdir, fnm)

    if not os.path.exists(fpath_f):
        print('path [%s] not exists.' % fpath_f)
        return

    sz_kb = os.path.getsize(fpath_f) / 1024.0
    if sz_kb < 5:
        print('[%s] is too short. only [%4.1f]kb.' % (fpath_f, sz_kb))
        return

    if fnm.endswith('.flt'):
        if not true_tbl:
            print('WARNING: no true table available, skip .flt [%s]' % fpath_f)
            return []
        err_1site = read_1_fltfile(fpath_f, true_tbl, opt.interval_sec, opt.start_time_str, opt.end_time_str)
    else:
        err_1site = read_1_posfile(fpath_f, opt.interval_sec, opt.start_time_str, opt.end_time_str)

    err_1site.dir0 = fdir
    err_1site.fname = fnm
    err_1site.stname = fnm[0:4]
    err_1site.ses_num = 999

    if len(err_1site.lst_denu_bases) < 3:
        print('[%s] no enough epoch. skip' % fnm)
        return []

    lst_err_res = divide_ep_ses_with_init_freq(err_1site, opt)

    return lst_err_res


def gen_cvg_statictics_file_1mode_multiprocess(strmd, lst_dir_names,
                                               opt: cvg_prc_opt_t, true_tbl):
    print('call gen_cvg_statictics_file_1mode_multiprocess for [%s]' % strmd)

    lst_lst_err_res = []
    with multiprocessing.Pool() as p:
        pt = partial(get_pos_result_ses_by_ses_1file, opt, true_tbl)
        lst_lst_err_res = p.map(pt, lst_dir_names)

    pth_series = os.path.join(opt.dirres, 'cvg_series-mode-%s.txt' % strmd)
    pth_cvg = os.path.join(opt.dirres, 'cvg_time-mode-%s.txt' % strmd)

    with open(pth_series, 'w') as fp_series, open(pth_cvg, 'w') as fp_cvg:
        s0hd = '%%               FILE                                     SES   EPNUM    NG    NR    NE    NJ    NC       cvg-E   cvg-N   cvg-U    cvg-H   cvg-3D/min'
        fp_cvg.write(s0hd + '\n')

        dir0 = ''

        lst_site = []
        for lst_err, dirnm in zip(lst_lst_err_res, lst_dir_names):
            if not lst_err:
                continue

            if dir0 != dirnm[0]:
                fp_cvg.write('##%%%% [folder]: %s\n' % dirnm[0])
                fp_series.write('##%%%% [folder]: %s\n' % dirnm[0])
                dir0 = dirnm[0]

            write_pos_session_res(lst_err, fp_cvg, fp_series)
            for err in lst_err:
                lst_site.append(err.stname)

    lst_site = list(set(lst_site))
    lst_site.sort()

    return pth_cvg, pth_series, lst_site


def run_cvg_analyze_(dict_md_files, opt: cvg_prc_opt_t, first_conv_period, true_tbl):
    print('call run_cvg_analyze_')

    lst_sites = []
    lst_modes, lst_pth_cvg_each_mode, lst_pth_series_each_mode = [], [], []
    for strmd, lst_dir_nm_pair in dict_md_files.items():
        pth_cvg, pth_series, lst_st = gen_cvg_statictics_file_1mode_multiprocess(strmd, lst_dir_nm_pair, opt, true_tbl)
        lst_modes.append(strmd)
        lst_pth_cvg_each_mode.append(pth_cvg)
        lst_pth_series_each_mode.append(pth_series)
        lst_sites.extend(lst_st)

    lst_sites = list(set(lst_sites))
    lst_sites.sort()

    dirres_site_cvg = os.path.join(opt.dirres, 'ave-cvg')
    if not os.path.exists(dirres_site_cvg):
        os.makedirs(dirres_site_cvg)
    sum_write_cvg_all_modes_multi_files(lst_modes, lst_pth_cvg_each_mode, dirres_site_cvg, first_conv_period)

    dirout, _ = os.path.split(lst_pth_series_each_mode[0])

    dir_cvg_pert = os.path.join(dirout, 'cvg_perct_txt')
    if not os.path.exists(dir_cvg_pert):
        os.makedirs(dir_cvg_pert)
    dir_cvg_series_each_site = os.path.join(dirout, 'cvg_series_each_site')
    if not os.path.exists(dir_cvg_series_each_site):
        os.makedirs(dir_cvg_series_each_site)

    is_bad_skip = False

    lst_perct = [50, 68, 90, 95]

    # 所有测站整体统计
    for pct in lst_perct:
        read_draw_cvg_series_multi_files(lst_pth_series_each_mode,lst_modes,dirout,dir_cvg_pert,dir_cvg_series_each_site,is_bad_skip,opt.init_freq_min,'',pct)

    # 按测站分别统计输出（开关由 opt.is_sitebysite_output 控制）
    if opt.is_sitebysite_output and lst_sites:
        print('site-by-site output enabled, total sites:', len(lst_sites))
        for st in lst_sites:
            for pct in lst_perct:
                read_draw_cvg_series_multi_files(lst_pth_series_each_mode,lst_modes,dirout,dir_cvg_pert,dir_cvg_series_each_site,is_bad_skip,opt.init_freq_min,st,pct)


def simplify_all_mode_str(lst_dir_name_md_inf, lst_md_input):
    lst_mode_no_rpt_src = []
    for tmp in lst_dir_name_md_inf:
        if tmp.md_src not in lst_mode_no_rpt_src:
            lst_mode_no_rpt_src.append(tmp.md_src)

    is_smp_need = False
    if len(lst_mode_no_rpt_src) > 1 and not lst_md_input:
        is_smp_need = True

    lst_modes_all = [tmp.md_src for tmp in lst_dir_name_md_inf]

    if is_smp_need:
        lst_modes_smp_all = simplify_mode_str_list_remove_cmn_str(lst_modes_all)
    else:
        lst_modes_smp_all = lst_modes_all

    for idx, md in enumerate(lst_modes_smp_all):
        lst_modes_smp_all[idx] = rename_mode_str(md)

    lst_mode_smp_no_rpt = []
    for tmp in lst_modes_smp_all:
        if tmp not in lst_mode_smp_no_rpt:
            lst_mode_smp_no_rpt.append(tmp)

    lst_mode_smp_no_rpt.sort()

    dict_md_files = dict()
    for tmp in lst_mode_smp_no_rpt:
        tmp1 = replace_special_mode_substring(tmp)
        dict_md_files[tmp1] = []

    for md, inf in zip(lst_modes_smp_all, lst_dir_name_md_inf):
        tmp1 = replace_special_mode_substring(md)
        dict_md_files[tmp1].append((inf.dir0, inf.name))

    return dict_md_files


def run_cvg_analyze_dirlst(dirlst, opt: cvg_prc_opt_t, first_conv_period, lst_md_input=[], true_file=''):
    start = time.time()
    print('call run_cvg_analyze_dirlst')

    if not check_dirlst(dirlst):
        return

    rm_mk_dir(opt.dirres)

    # Load true coordinates if provided
    true_tbl = {}
    if true_file and os.path.isfile(true_file):
        true_tbl = load_true_coords(true_file)
        print('Loaded true coordinates for %d sites from [%s]' % (len(true_tbl), true_file))
    else:
        print('No true coordinate file provided. .flt files will be skipped.')

    lst_dir_name_md_inf = get_dir_file_pair_md_list(dirlst, lst_md_input)

    if not lst_dir_name_md_inf:
        print('no valid pos files found. return')
        return

    dict_md_files = simplify_all_mode_str(lst_dir_name_md_inf, lst_md_input)
    run_cvg_analyze_(dict_md_files, opt, first_conv_period, true_tbl)

    end = time.time()
    print('time spend : %.2f' % (end - start))


def trans_ppprtk2pos_1dir(dir0, md, dirtemp, is_pppar=False):
    ext = '.ppprtk'
    if is_pppar:
        ext = '.pppar'

    for nm in os.listdir(dir0):
        if not nm.endswith(ext):
            continue
        if '-augoff' in nm:
            print(f'skip [{nm}] due to no aug sol')
            continue

        fpth_prtk = os.path.join(dir0, nm)
        fpth_pos = os.path.join(dir0, nm.replace(ext, '.pos'))
        if not os.path.exists(fpth_pos):
            input('ppprtk file [%s] not exists' % fpth_pos)

        nm_temp = nm.split('-PPP_')[0] + '-PPP_' + md + '.pos'
        fpth_temp = os.path.join(dirtemp, nm_temp)
        with open(fpth_temp, 'w') as fp:
            lines0 = open(fpth_prtk).readlines()
            lines1 = open(fpth_pos).readlines()

            fp.write(lines1[0])

            idx1 = 0

            for idx0, ln0 in enumerate(lines0):
                if not check_is_pos_line(ln0):
                    continue
                s0tm, s0err = ln0[0:22].rstrip(), ln0[49:72].rstrip()

                for i, ln1 in enumerate(lines1[idx1:], idx1):
                    if ln1.startswith(s0tm):
                        ss1 = ln1.split()
                        if len(ss1) >= 17:
                            s1pre = ' '.join(ss1[:17])
                        else:
                            s1pre = ln1[:149]
                        idx1 = i + 1

                        fp.write(s1pre + ' ' + s0err + '   1.000' + '\n')
                        break
                    if i - idx1 > 20:
                        break


def rename_float_pos_1dir(dir0, dirtemp):
    for nm in os.listdir(dir0):
        if not nm.endswith('.ppprtk'):
            continue
        fpth_src = os.path.join(dir0, nm)

        nm_temp = nm.split('-PPP_')[0] + '-PPP_0=Flo.ppprtk'
        fpth_dst = os.path.join(dirtemp, nm_temp)

        shutil.copy(fpth_src, fpth_dst)


def run_cvg_analyze_4ppprtk(lst_dir_md, opt, dirtemp):
    lst_dir, lst_md = [], []
    for tmp in lst_dir_md:
        if not os.path.exists(tmp[0]):
            continue
        lst_dir.append(tmp[0])
        lst_md.append(tmp[1])

    if not lst_dir:
        print('no valid ppprtk result dir. return')
        return

    rm_mk_dir(dirtemp)

    count = Counter(lst_md).most_common(1)
    md_most = count[0][0]
    print(md_most)

    for dir0, md in zip(lst_dir, lst_md):
        if md != md_most:
            continue
        rename_float_pos_1dir(dir0, dirtemp)
        trans_ppprtk2pos_1dir(dir0, '1=AR', dirtemp, True)

    for dir0, md in zip(lst_dir, lst_md):
        trans_ppprtk2pos_1dir(dir0, md, dirtemp)

    run_cvg_analyze_dirlst([dirtemp], opt, first_conv_period, true_file=true_file)


if __name__ == "__main__":
    opt = cvg_prc_opt_t()
    opt.init_freq_min = 60 * 3
    enable_site_by_site_output = False

    opt.is_sitebysite_output = enable_site_by_site_output
    opt.interval_sec = 30

    # 设置解析数据的起止时间，支持格式:
    #   1) 完整日期时间: '2023/03/11 01:00:00'
    #   2) 仅时间: '01:00:00' (自动使用该文件首个历元的日期补全)
    # 留空字符串表示不限制
    opt.start_time_str = ''
    opt.end_time_str = ''

    dir_00 = r'/home/lzx/code/GREAT-PVT/data/2023305/products/result_stec'
    dir_01 = r'/home/lzx/code/GREAT-PVT/data/2023305/products/result_zwd'
    dir_02 = r'/home/lzx/code/GREAT-PVT/data/2023305/products/result_stec_zwd'
    dir_03 = r'/home/lzx/code/GREAT-PVT/data/2023305/products/result_KIN_UC_GE'
    # dir_04 = r'/home/zxli/code/GKit_uducppp-master/data/GA_2023132_ppprtk_20/result-uduc_dim'
    # dir_05 = r'/home/zxli/code/GKit_uducppp-master/data/GA_2023132_ppprtk_20/result-uduc_rdcb_single_sigma005'
    dirdst = r'/home/lzx/code/GREAT-PVT/data/2023305/compare'

    #dir_00 = r'/stars/multi-freq-PPPRTK/GA-2023070-2023079/All_result/select_ppprtk/result-uduc-usr-GEC-2f-355-200/'
    #dir_01 = r'/stars/multi-freq-PPPRTK/GA-2023070-2023079/All_result/select_ppprtk/result-uduc-usr-GEC-s2u5-355-200/'
    #dir_02 = r'/stars/multi-freq-PPPRTK/GA-2023070-2023079/All_result/select_ppprtk/result-uduc-usr-GEC-s5u2-355-200/'
    #dir_03 = r'/stars/multi-freq-PPPRTK/GA-2023070-2023079/All_result/select_ppprtk/result-uduc-usr-GEC-5f-355-200/'
    #dirdst = r'/stars/multi-freq-PPPRTK/GA-2023070-2023079/All_result/select_ppprtk/picture/cvg_tests/'

    lst_dir = [dir_00,dir_01,dir_02,dir_03]
    lst_md = ['stec', 'zwd','stec_zwd','stec_zwd_noZwdInit']
    # lst_dir = [dir_00]
    # lst_md = ['stec_zwd']
    opt.dirres = dirdst

    first_conv_period = False
    
    #lst_md = ['S2U2','S2U5','S5U2','S5U5']

    # 真值坐标文件路径（用于 .flt 文件 ENU 误差计算），默认从工作目录加载
    true_file = os.path.abspath('true_crd.true_crd')

    # 检测输入目录是否包含 .flt 文件（且不含 .ppprtk）
    has_flt = any(any(f.endswith('.flt') for f in os.listdir(d)) for d in lst_dir if os.path.isdir(d))
    has_ppprtk = any(any(f.endswith('.ppprtk') for f in os.listdir(d)) for d in lst_dir if os.path.isdir(d))

    if has_flt and not has_ppprtk:
        # 直接分析原始目录中的 .flt 文件，跳过临时目录中转
        # 规则：第一个目录分析 float / pppar / ppprtk；后续目录只分析 ppprtk
        opt.dirres = dirdst
        rm_mk_dir(opt.dirres)

        true_tbl = {}
        if true_file and os.path.isfile(true_file):
            true_tbl = load_true_coords(true_file)
            print('Loaded true coordinates for %d sites from [%s]' % (len(true_tbl), true_file))
        else:
            print('No true coordinate file provided. .flt files will be skipped.')

        dict_md_files = {}

        # 第一个目录：float, pppar, ppprtk（ppprtk 使用 lst_md[0] 作为标签）
        dir0 = lst_dir[0]
        md0_ppprtk = lst_md[0] if len(lst_md) > 0 else 'ppprtk'
        for nm in os.listdir(dir0):
            if not nm.endswith('.flt') or '-augoff' in nm:
                continue
            if nm.endswith('_float.flt'):
                md = 'float'
            elif nm.endswith('_pppar.flt'):
                md = 'pppar'
            elif nm.endswith('_ppprtk.flt'):
                md = md0_ppprtk
            else:
                continue
            dict_md_files.setdefault(md, []).append((dir0, nm))

        # 后续目录：只分析 ppprtk，标签使用对应的 lst_md[i]
        for dir_i, md_i in zip(lst_dir[1:], lst_md[1:]):
            for nm in os.listdir(dir_i):
                if not nm.endswith('_ppprtk.flt') or '-augoff' in nm:
                    continue
                dict_md_files.setdefault(md_i, []).append((dir_i, nm))

        run_cvg_analyze_(dict_md_files, opt, first_conv_period, true_tbl)
    else:
        # 准备临时目录，把 flo / ar / 各目录 ppprtk 统一平铺，按 mode 重命名
        dirtemp_mixed = r'/home/lzx/code/GREAT-PVT/data/GA_2023070/compare_mixed'
        rm_mk_dir(dirtemp_mixed)

        # Flo: 第一个目录的 .ppprtk -> 0=Flo
        for nm in os.listdir(lst_dir[0]):
            if not nm.endswith('.ppprtk') or '-augoff' in nm:
                continue
            src = os.path.join(lst_dir[0], nm)
            dst = os.path.join(dirtemp_mixed, nm.split('-PPP_')[0] + '-PPP_0=Flo.ppprtk')
            shutil.copy(src, dst)

        # AR: 第一个目录的 .pppar -> 1=AR
        for nm in os.listdir(lst_dir[0]):
            if not nm.endswith('.pppar') or '-augoff' in nm:
                continue
            src = os.path.join(lst_dir[0], nm)
            dst = os.path.join(dirtemp_mixed, nm.split('-PPP_')[0] + '-PPP_1=AR.ppprtk')
            shutil.copy(src, dst)

        # 各目录的 ppprtk
        for d0, md_name in zip(lst_dir, lst_md):
            for nm in os.listdir(d0):
                if not nm.endswith('.ppprtk') or '-augoff' in nm:
                    continue
                src = os.path.join(d0, nm)
                dst = os.path.join(dirtemp_mixed, nm.split('-PPP_')[0] + '-PPP_' + md_name + '.ppprtk')
                shutil.copy(src, dst)

        # 统一分析，结果输出到 compare/
        opt.dirres = dirdst
        run_cvg_analyze_dirlst([dirtemp_mixed], opt, first_conv_period, true_file=true_file)

    print('All files processed successfully.')