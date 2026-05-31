#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制 .pos / .flt 结果文件相对于 true_crd.true_crd 的三维误差曲线（E/N/U 随时间）。

用法:
  python scripts/plot_pos_error.py INPUT_PATH [--true TRUE_FILE] [--site SITE]
                                 [--save-combined OUT.png]
                                 [--save-faceted OUT.png]
                                 [--save-individual-dir OUT_DIR]
                                 [--conv-thresh 0.1] [--conv-win 5] [--conv-mode 3d|hz]

- INPUT_PATH: 单个 .pos/.flt 文件，或包含 .pos/.flt 文件的目录（递归搜索）
- --true: 真实坐标文件路径（默认: ./true_crd.true_crd）
- --site: 指定台站代码（否则尝试从文件名或头部推断）
- --save-combined/--save-faceted/--save-individual-dir: 指定保存路径/目录

默认行为：
- 若未显式指定保存参数，脚本会在输入路径下创建 plot/ 目录，并同时输出：
  - plot/enu_combined.png（合并图：所有文件的 E/N/U 三条曲线分别叠加）
  - plot/enu_faceted.png（分面图：每个文件一张子图，内含 E/N/U 三条曲线）
  - plot/<文件名>_enu.png（每个文件单独一张图，包含 E/N/U 三条曲线）

收敛时间统计：
- 默认统计“连续 5 个历元三维误差 < 0.1 m”的首次出现，打印到终端（可用参数覆盖）。

依赖: Python 3.8+, numpy, matplotlib, pandas(可选)
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple


try:
    import numpy as np
except Exception as e:
    print("ERROR: numpy is required. Please install it: pip install numpy", file=sys.stderr)
    raise

# pandas and matplotlib are optional until plotting
try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)

@dataclass
class TrueCoord:
    site: str
    x: float
    y: float
    z: float
    lat_deg: float
    lon_deg: float
    hgt_m: float


def dms_to_deg(s: str) -> float:
    # Legacy helper, currently unused; kept for backward compatibility
    return float(s)


def _snx_load_true_coords(true_file: str) -> Dict[str, TrueCoord]:
    """Parse a SINEX file (e.g. igs21P2185.snx) and build TrueCoord table.

    Strategy: use SOLUTION/ESTIMATE block and collect STAX/STAY/STAZ as ECEF XYZ.
    """
    coords_tmp: Dict[str, Dict[str, float]] = {}
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
            # *INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__
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

    true_tbl: Dict[str, TrueCoord] = {}
    for site, rec in coords_tmp.items():
        if not all(k in rec for k in ('x', 'y', 'z')):
            continue
        x = rec['x']
        y = rec['y']
        z = rec['z']
        lat_deg, lon_deg, hgt_m = ecef_to_geodetic(x, y, z)
        true_tbl[site] = TrueCoord(site, x, y, z, lat_deg, lon_deg, hgt_m)

    if not true_tbl:
        raise ValueError(f"No station ECEF coordinates parsed from SINEX file: {true_file}")
    return true_tbl


def geodetic_to_ecef(lat_deg: float, lon_deg: float, h: float) -> Tuple[float, float, float]:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)
    N = WGS84_A / np.sqrt(1 - WGS84_E2 * sin_lat * sin_lat)
    x = (N + h) * cos_lat * cos_lon
    y = (N + h) * cos_lat * sin_lon
    z = (N * (1 - WGS84_E2) + h) * sin_lat
    return x, y, z


def ecef_to_enu_matrix(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    sl, cl = np.sin(lat), np.cos(lat)
    sb, cb = np.sin(lon), np.cos(lon)
    # ENU from ECEF rotation matrix
    R = np.array([[-sb,          cb,          0],
                  [-sl*cb,      -sl*sb,       cl],
                  [ cl*cb,       cl*sb,       sl]])
    return R


def ecef_to_geodetic(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Convert ECEF (meters) to geodetic coordinates (lat_deg, lon_deg, h).
    Uses an iterative method (Bowring) for latitude convergence.
    """
    # constants
    a = WGS84_A
    e2 = WGS84_E2
    lon = float(np.arctan2(y, x))
    p = float(np.sqrt(x * x + y * y))
    # initial guess of latitude
    lat = float(np.arctan2(z, p * (1 - e2)))
    for _ in range(10):
        sin_lat = np.sin(lat)
        N = a / np.sqrt(1 - e2 * sin_lat * sin_lat)
        h = p / np.cos(lat) - N
        lat_new = float(np.arctan2(z, p * (1 - e2 * (N / (N + h)))))
        if abs(lat_new - lat) < 1e-12:
            lat = lat_new
            break
        lat = lat_new
    sin_lat = np.sin(lat)
    N = a / np.sqrt(1 - e2 * sin_lat * sin_lat)
    h = p / np.cos(lat) - N
    return float(np.rad2deg(lat)), float(np.rad2deg(lon)), float(h)


_DATA_EXTENSIONS = ('.pos', '.flt')


def find_data_files(path: str) -> List[str]:
    if os.path.isdir(path):
        got: List[str] = []
        for root, _, files in os.walk(path):
            for f in files:
                if f.lower().endswith(_DATA_EXTENSIONS):
                    got.append(os.path.join(root, f))
        return sorted(got)
    else:
        return [path] if path.lower().endswith(_DATA_EXTENSIONS) else []


def load_true_coords(true_file: str) -> Dict[str, TrueCoord]:
    tbl: Dict[str, TrueCoord] = {}
    if not os.path.isfile(true_file):
        raise FileNotFoundError(f"True coordinate file not found: {true_file}")
    # Auto-detect SINEX format by file extension
    if true_file.lower().endswith('.snx'):
        return _snx_load_true_coords(true_file)
    with open(true_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line.startswith('%'):
                continue
            # Expect: SITE X Y Z lat lon h ...
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
            tbl[site] = TrueCoord(site, x, y, z, lat, lon, hgt)
    return tbl


def infer_site_from_filename(path: str, true_tbl: Dict[str, TrueCoord]) -> Optional[str]:
    base = os.path.basename(path)
    name = os.path.splitext(base)[0]
    # Try tokens split by non-alnum, pick any token that matches a site
    tokens = re.split(r"[^A-Za-z0-9]+", name)
    tokens = [t.upper() for t in tokens if t]
    # Prefer exact length-4 tokens common in GNSS site codes, but keep all
    for t in tokens:
        if t in true_tbl:
            return t
    # Fallback: if unique site appears as substring
    for site in true_tbl.keys():
        if site in name.upper():
            return site
    return None


def infer_site_from_header(path: str, true_tbl: Dict[str, TrueCoord]) -> Optional[str]:
    # Search the first several hundred lines for any occurrence of a known site code.
    # Some .pos files put the site in headers, others embed it differently; be permissive.
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if i > 500:
                    break
                u = line.upper()
                for site in true_tbl.keys():
                    if site in u:
                        return site
    except Exception:
        pass
    return None


@dataclass
class PosSeries:
    time: List[Any]  # pandas.Timestamp or float index fallback
    ecef: np.ndarray  # shape (N,3)
    time_label: str = 'Time'


def parse_time_safe(cols: List[str]) -> Optional[object]:
    if pd is None:
        return None
    s = ' '.join(cols[:2]) if len(cols) >= 2 else cols[0]
    try:
        return pd.to_datetime(s, errors='coerce')
    except Exception:
        return None


def detect_format_values(vals: List[float]) -> str:
    # Return 'LLH' or 'ECEF' based on magnitudes: ECEF ~ 6e6, LLH lat within +/-90, lon within +/-360
    if len(vals) < 3:
        return 'UNKNOWN'
    v0, v1, v2 = vals[0], vals[1], vals[2]
    # If deg plausible
    if -90.0 <= v0 <= 90.0 and -360.0 <= v1 <= 360.0 and -2000.0 <= v2 <= 10000.0:
        return 'LLH'
    # If ECEF magnitude
    if max(abs(v0), abs(v1), abs(v2)) > 1e5:
        return 'ECEF'
    return 'UNKNOWN'


def parse_pos_file(path: str) -> PosSeries:
    times: List[Any] = []
    ecef_list: List[Tuple[float, float, float]] = []
    # Try to detect columns by reading few lines
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line[0] in ('%', '#'):
                continue
            cols = line.strip().split()
            # Special handling for GAMP-style POS format:
            # YYYY MM DD hh mm ss week sow X Y Z E N U 3D
            if len(cols) >= 11:
                try:
                    y = int(cols[0]); mo = int(cols[1]); d = int(cols[2])
                    hh = int(cols[3]); mi = int(cols[4]); ss = float(cols[5])
                    if 1900 <= y <= 2100 and 1 <= mo <= 12 and 1 <= d <= 31 \
                       and 0 <= hh < 24 and 0 <= mi < 60 and 0.0 <= ss < 61.0:
                        x = float(cols[8]); y_e = float(cols[9]); z_e = float(cols[10])
                        ecef_list.append((x, y_e, z_e))
                        if pd is not None:
                            try:
                                # Build a timestamp like '2017-09-01 00:00:30'
                                tstr = f"{y:04d}-{mo:02d}-{d:02d} {hh:02d}:{mi:02d}:{int(ss):02d}"
                                tstamp = pd.to_datetime(tstr, errors='coerce')
                            except Exception:
                                tstamp = None
                            times.append(tstamp if (tstamp is not None and pd.notna(tstamp)) else pd.NaT)
                        else:
                            times.append(float(len(ecef_list)))
                        continue
                except Exception:
                    # fall back to generic parser below
                    pass
            # Try parse three floats from either LLH or ECEF in common positions
            # Heuristics: if first tokens look like date/time, skip them for values
            start_idx = 0
            # Try parse time with pandas if available
            tstamp: Optional[object] = None
            if pd is not None:
                # Try 2-col time
                if len(cols) >= 2:
                    ttry = parse_time_safe(cols[:2])
                    if ttry is not None and pd.notna(ttry):
                        tstamp = ttry
                        start_idx = 2
                if tstamp is None and len(cols) >= 1:
                    ttry = parse_time_safe(cols[:1])
                    if ttry is not None and pd.notna(ttry):
                        tstamp = ttry
                        start_idx = 1
            else:
                # pandas not available: apply simple heuristics to detect date/time tokens
                # Common formats: 'YYYY/MM/DD' 'YYYY-MM-DD' followed by 'HH:MM:SS(.sss)'
                if len(cols) >= 2 and (('/' in cols[0] or '-' in cols[0]) and (':' in cols[1])):
                    start_idx = 2
                # Single-token datetime like '2024/05/11T00:00:01' or '2024-05-11T00:00:01.000'
                elif len(cols) >= 1 and (('/' in cols[0] or '-' in cols[0]) and (':' in cols[0] or 'T' in cols[0])):
                    start_idx = 1
            # Extract numeric values from the remaining
            vals: List[float] = []
            for tok in cols[start_idx:]:
                try:
                    vals.append(float(tok))
                except Exception:
                    # Stop at first non-float after numeric series
                    break
            if len(vals) < 3:
                continue
            kind = detect_format_values(vals[:3])
            if kind == 'LLH':
                lat, lon, h = vals[0], vals[1], vals[2]
                x, y, z = geodetic_to_ecef(lat, lon, h)
            elif kind == 'ECEF':
                x, y, z = vals[0], vals[1], vals[2]
            else:
                # Can't detect, skip line
                continue
            ecef_list.append((x, y, z))
            if pd is not None:
                times.append(tstamp if tstamp is not None else pd.NaT)
            else:
                times.append(float(len(ecef_list)))
    if not ecef_list:
        raise ValueError(f"No position rows parsed from {path}")
    ecef = np.asarray(ecef_list, dtype=float)
    return PosSeries(time=times, ecef=ecef)


def parse_flt_file(path: str) -> PosSeries:
    """Parse GREAT-PVT .flt filter output file.

    Expected columns (space-separated):
      sow(s)  X-ECEF(m)  Y-ECEF(m)  Z-ECEF(m)  Vx  Vy  Vz  X-RMS  Y-RMS  Z-RMS ...
    Comment lines start with '#'.
    Time is stored as GPS Seconds of Week (float) for the x-axis.
    """
    times: List[float] = []
    ecef_list: List[Tuple[float, float, float]] = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line.strip().startswith('#'):
                continue
            cols = line.strip().split()
            if len(cols) < 4:
                continue
            try:
                sow = float(cols[0])
                x = float(cols[1])
                y = float(cols[2])
                z = float(cols[3])
            except Exception:
                continue
            times.append(sow)
            ecef_list.append((x, y, z))
    if not ecef_list:
        raise ValueError(f"No position rows parsed from {path}")
    ecef = np.asarray(ecef_list, dtype=float)
    return PosSeries(time=times, ecef=ecef, time_label='GPS Seconds of Week (s)')


def compute_enu_error(series: PosSeries, true_ecef: Tuple[float, float, float], true_llh: Tuple[float, float, float]) -> Tuple[np.ndarray, Dict[str, float]]:
    e_ref = np.array(true_ecef, dtype=float).reshape(3)
    lat_deg, lon_deg, hgt = true_llh
    R = ecef_to_enu_matrix(lat_deg, lon_deg)
    diff = series.ecef - e_ref
    enu = (R @ diff.T).T  # shape (N,3)
    # Ensure a plain numpy array is returned to avoid pandas extension-array
    enu = np.asarray(enu, dtype=float)
    # stats
    E, N, U = enu[:, 0], enu[:, 1], enu[:, 2]
    def rms(a: np.ndarray) -> float:
        return float(np.sqrt(np.mean(a**2)))
    stats = {
        'rmsE': rms(E), 'rmsN': rms(N), 'rmsU': rms(U),
        'meanE': float(np.mean(E)), 'meanN': float(np.mean(N)), 'meanU': float(np.mean(U)),
        'maxAbsE': float(np.max(np.abs(E))), 'maxAbsN': float(np.max(np.abs(N))), 'maxAbsU': float(np.max(np.abs(U))),
    }
    return enu, stats


def compute_convergence_time(enu: np.ndarray, series: PosSeries, thresh: float = 0.1, win: int = 5, mode: str = '3d') -> Tuple[Optional[int], Optional[object]]:
    """统计首次收敛时间：在长度为 win 的滑动窗口内，误差始终小于阈值 thresh。
    mode: '3d' -> sqrt(E^2+N^2+U^2), 'hz' -> sqrt(E^2+N^2)
    返回 (index, time)，index 为第一个满足窗口的起始索引；
    time 为从第一个历元开始的相对秒数（对 .flt 的 sow 取差值）。
    """
    if enu.size == 0:
        return None, None
    E, N, U = enu[:, 0], enu[:, 1], enu[:, 2]
    if mode.lower() == 'hz':
        err = np.sqrt(E*E + N*N)
    else:
        err = np.sqrt(E*E + N*N + U*U)
    if len(err) < win:
        return None, None
    below = err < thresh
    consec = 0
    idx_found: Optional[int] = None
    for i, ok in enumerate(below):
        consec = consec + 1 if ok else 0
        if consec >= win:
            idx_found = i - win + 1  # 窗口起始索引
            break
    if idx_found is None:
        return None, None
    t = None
    if isinstance(series.time, list) and idx_found < len(series.time):
        tv = series.time[idx_found]
        # For plain numeric times (e.g. GPS sow), return relative time from first epoch
        if isinstance(tv, (int, float)):
            t0 = series.time[0]
            if isinstance(t0, (int, float)):
                t = float(tv) - float(t0)
            else:
                t = float(tv)
        elif pd is not None:
            try:
                t = pd.to_datetime(tv)
            except Exception:
                t = None
    return idx_found, t


def _ylims_robust(values: np.ndarray, min_span: float = 0.05, n_sigma: float = 2.0) -> Tuple[float, float]:
    """计算以 0 为中心的稳健纵轴范围，用 MAD 替代标准差以抑制异常跳变影响。"""
    if values.size == 0:
        return -0.1, 0.1
    # 中位数绝对偏差 (MAD) -> 稳健标准差
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    robust_sigma = 1.4826 * mad  # 正态分布下 MAD * 1.4826 ≈ std
    # 同时考虑最大绝对值，但用 percentile 避免极端异常值
    p99 = float(np.percentile(np.abs(values), 99.0))
    span = max(n_sigma * robust_sigma, p99 * 1.1, min_span)
    # ENU 误差图以 0 为中心
    return -span, span


def plot_enu_series(enu_list: List[Tuple[str, PosSeries, np.ndarray, Dict[str, float]]], save_path: Optional[str] = None, ylim_en: Optional[float] = None, ylim_u: Optional[float] = None) -> None:
    if plt is None:
        print("ERROR: matplotlib is required to plot. pip install matplotlib", file=sys.stderr)
        return
    # Build time axis: prefer pandas timestamps if available
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    labels = ['East (m)', 'North (m)', 'Up (m)']
    colors = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0','C1','C2','C3','C4'])

    # Precompute 2-sigma y-limits across all series per component
    if len(enu_list) > 0:
        # Ensure each enu is a numpy array before concatenation to avoid
        # pandas extension/sequence objects that disallow multidim indexing.
        E_all = np.concatenate([np.asarray(enu, dtype=float)[:, 0] for _, __, enu, ___ in enu_list])
        N_all = np.concatenate([np.asarray(enu, dtype=float)[:, 1] for _, __, enu, ___ in enu_list])
        U_all = np.concatenate([np.asarray(enu, dtype=float)[:, 2] for _, __, enu, ___ in enu_list])
        e_lim = _ylims_robust(E_all)
        n_lim = _ylims_robust(N_all)
        u_lim = _ylims_robust(U_all)
    else:
        e_lim = n_lim = u_lim = (-0.1, 0.1)

    # Override with user-specified limits if provided
    if ylim_en is not None:
        e_lim = n_lim = (-ylim_en, ylim_en)
    if ylim_u is not None:
        u_lim = (-ylim_u, ylim_u)

    for idx, (name, series, enu, stats) in enumerate(enu_list):
        t = series.time
        if isinstance(t, list) and len(t) > 0:
            # 纯数字时间（如 GPS sow）直接用数值数组
            if all(isinstance(ts, (int, float)) for ts in t):
                x = np.array(t, dtype=float)
            elif pd is not None:
                # 尝试 pandas datetime
                if any(getattr(ts, 'value', None) is None for ts in t):
                    x = np.arange(len(t))
                else:
                    x = pd.to_datetime(t)
            else:
                x = np.arange(len(t))
        else:
            x = np.arange(len(enu))
        color = colors[idx % len(colors)]
        # Coerce enu to numpy in case it's a list/sequence or pandas-backed array
        enu_arr = np.asarray(enu, dtype=float)
        # 统一图例标签（不含方向 RMS，避免重复）
        label = f"{name}"
        _plot_with_gaps(axes[0], x, enu_arr[:, 0], gap_seconds=120.0, label=label, color=color)
        _plot_with_gaps(axes[1], x, enu_arr[:, 1], gap_seconds=120.0, label=label, color=color)
        _plot_with_gaps(axes[2], x, enu_arr[:, 2], gap_seconds=120.0, label=label, color=color)

    for ax, lab in zip(axes, labels):
        ax.set_ylabel(lab)
        ax.grid(True, linestyle='--', alpha=0.4)
        # 加粗零刻度参考线
        ax.axhline(0.0, color='k', linewidth=1.2, alpha=0.7, zorder=1)
        # 不在子图内显示图例
    axes[0].set_ylim(*e_lim)
    axes[1].set_ylim(*n_lim)
    axes[2].set_ylim(*u_lim)
    axes[-1].set_xlabel(enu_list[0][1].time_label if enu_list else 'Time')
    fig.suptitle('Positioning Error (ENU) relative to true_crd.true_crd')

    # 在 figure 级别统一放置图例（右上角，去重）
    from collections import OrderedDict
    handles_all, labels_all = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        handles_all.extend(h)
        labels_all.extend(l)
    unique = OrderedDict()
    for h, l in zip(handles_all, labels_all):
        if l not in unique:
            unique[l] = h
    if unique:
        fig.legend(unique.values(), unique.keys(), loc='upper right',
                   fontsize=7, frameon=True, ncol=1,
                   bbox_to_anchor=(0.99, 0.97))

    fig.tight_layout(rect=[0, 0.03, 1, 0.97])

    if save_path:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")
    else:
        plt.show()


def _time_to_x(series: PosSeries):
    if not isinstance(series.time, list) or len(series.time) == 0:
        return np.arange(len(series.ecef))
    # If all elements are plain numbers (e.g. GPS sow), use them directly
    if all(isinstance(t, (int, float)) for t in series.time):
        return np.array(series.time, dtype=float)
    if pd is not None:
        # Check for pandas Timestamp objects
        has_ts = any(hasattr(ts, 'value') for ts in series.time)
        if has_ts:
            if any(getattr(ts, 'value', None) is None for ts in series.time):
                return np.arange(len(series.time))
            return pd.to_datetime(series.time)
    return np.arange(len(series.ecef))


def _plot_with_gaps(ax, x, y, gap_seconds: float = 120.0, label: Optional[str] = None, **plot_kwargs):
    """Plot x vs y on ax but break the line where consecutive time gaps exceed gap_seconds.

    - If pandas is available and x can be converted to datetimes, splits by gaps and plots
      each continuous segment separately. Only the first segment receives the legend label
      to avoid duplicate legend entries.
    - Falls back to a single continuous plot when datetimes are not available/parsable.
    """
    plotted_any = False
    # Ensure y is a numpy array to avoid pandas' recent restriction on multidimensional
    # indexing (e.g. obj[:, None]). This makes subsequent slicing safe.
    try:
        y = np.asarray(y)
    except Exception:
        pass
    # Heuristic: if x is a plain numeric array with small positive values,
    # treat it as GPS sow or epoch index and skip pd.to_datetime to avoid
    # misinterpreting values as nanosecond timestamps.
    skip_datetime = False
    try:
        x_arr = np.asarray(x)
        if x_arr.dtype.kind in 'iuf' and x_arr.min() >= 0 and x_arr.max() < 1e7:
            skip_datetime = True
    except Exception:
        pass
    # Try to treat x as datetimes when pandas is available
    if not skip_datetime and pd is not None:
        try:
            xt = pd.to_datetime(x)
            # Ensure a Series for diff computation
            xs = xt.to_series(index=range(len(xt))) if not isinstance(xt, (pd.Series, pd.DatetimeIndex)) else pd.Series(xt)
            dif = xs.diff().dt.total_seconds().fillna(0).values
            seg_start = 0
            for i in range(1, len(dif)):
                if dif[i] > gap_seconds:
                    seg_x = xt[seg_start:i]
                    seg_y = y[seg_start:i]
                    # convert to numpy to avoid pandas extension arrays causing
                    # multidimensional indexing errors inside matplotlib
                    try:
                        seg_x = np.asarray(seg_x)
                    except Exception:
                        pass
                    try:
                        seg_y = np.asarray(seg_y)
                    except Exception:
                        pass
                    if len(seg_x) > 0:
                        kw = dict(plot_kwargs)
                        if not plotted_any and label:
                            kw['label'] = label
                        else:
                            kw['label'] = '_nolegend_'
                        ax.plot(seg_x, seg_y, **kw)
                        plotted_any = True
                    seg_start = i
            # final segment
            seg_x = xt[seg_start:len(xt)]
            seg_y = y[seg_start:len(xt)]
            try:
                seg_x = np.asarray(seg_x)
            except Exception:
                pass
            try:
                seg_y = np.asarray(seg_y)
            except Exception:
                pass
            if len(seg_x) > 0:
                kw = dict(plot_kwargs)
                if not plotted_any and label:
                    kw['label'] = label
                else:
                    kw['label'] = '_nolegend_'
                ax.plot(seg_x, seg_y, **kw)
                plotted_any = True
            return
        except Exception:
            # fall back to continuous plotting below
            pass
    # Fallback: continuous plot when no datetime handling available
    kw = dict(plot_kwargs)
    if label:
        kw['label'] = label
    # Ensure arrays to avoid pandas' new restriction on multidimensional indexing
    try:
        x = np.asarray(x)
    except Exception:
        pass
    try:
        y = np.asarray(y)
    except Exception:
        pass
    ax.plot(x, y, **kw)


def plot_enu_faceted(enu_list: List[Tuple[str, PosSeries, np.ndarray, Dict[str, float]]], save_path: Optional[str] = None, max_cols: int = 3, ylim_en: Optional[float] = None, ylim_u: Optional[float] = None) -> None:
    if plt is None:
        print("ERROR: matplotlib is required to plot. pip install matplotlib", file=sys.stderr)
        return
    import math
    n = len(enu_list)
    if n == 0:
        print("Nothing to plot.", file=sys.stderr)
        return
    cols = max(1, min(max_cols, n))
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 3*rows), squeeze=False)
    colors = {
        'E': 'C0',
        'N': 'C1',
        'U': 'C2',
    }
    for idx, (name, series, enu, stats) in enumerate(enu_list):
        r = idx // cols
        c = idx % cols
        ax = axes[r][c]
        x = _time_to_x(series)
        enu_arr = np.asarray(enu, dtype=float)
        _plot_with_gaps(ax, x, enu_arr[:, 0], gap_seconds=120.0, label=f"E (RMS={stats['rmsE']:.3f})", color=colors['E'])
        _plot_with_gaps(ax, x, enu_arr[:, 1], gap_seconds=120.0, label=f"N (RMS={stats['rmsN']:.3f})", color=colors['N'])
        _plot_with_gaps(ax, x, enu_arr[:, 2], gap_seconds=120.0, label=f"U (RMS={stats['rmsU']:.3f})", color=colors['U'])
        ax.set_title(name, fontsize=10)
        ax.set_ylabel('Error (m)')
        ax.set_xlabel(series.time_label)
        ax.grid(True, linestyle='--', alpha=0.4)
        # 加粗零刻度参考线
        ax.axhline(0.0, color='k', linewidth=1.2, alpha=0.7, zorder=1)
        ax.legend(loc='upper right', fontsize=8, ncol=3)
        # 2-sigma limits per component (per file)
        e_lim = _ylims_robust(enu_arr[:, 0])
        n_lim = _ylims_robust(enu_arr[:, 1])
        u_lim = _ylims_robust(enu_arr[:, 2])
        # Apply lims by axis by replotting? Here one axes has all three; choose a common symmetric lim across E/N/U
        common_vals = np.concatenate([enu_arr[:, 0], enu_arr[:, 1], enu_arr[:, 2]])
        y_lim = _ylims_robust(common_vals)
        # Override with user-specified limits if provided
        if ylim_en is not None:
            y_lim = (-ylim_en, ylim_en)
        ax.set_ylim(*y_lim)
    # Hide unused axes
    total_axes = rows * cols
    for j in range(n, total_axes):
        r = j // cols
        c = j % cols
        axes[r][c].axis('off')
    fig.suptitle('Positioning Error (ENU) - faceted by file')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_path:
        os.makedirs(os.path.dirname(save_path) or '.', exist_ok=True)
        plt.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")
    else:
        plt.show()


def plot_enu_individual(enu_list: List[Tuple[str, PosSeries, np.ndarray, Dict[str, float]]], out_dir: str, ylim_en: Optional[float] = None, ylim_u: Optional[float] = None) -> None:
    if plt is None:
        print("ERROR: matplotlib is required to plot. pip install matplotlib", file=sys.stderr)
        return
    os.makedirs(out_dir, exist_ok=True)
    for name, series, enu, stats in enu_list:
        fig, axes = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
        labels = ['East (m)', 'North (m)', 'Up (m)']
        x = _time_to_x(series)
        enu_arr = np.asarray(enu, dtype=float)
        _plot_with_gaps(axes[0], x, enu_arr[:, 0], gap_seconds=120.0, color='C0', label=f"E (RMS={stats['rmsE']:.3f})")
        _plot_with_gaps(axes[1], x, enu_arr[:, 1], gap_seconds=120.0, color='C1', label=f"N (RMS={stats['rmsN']:.3f})")
        _plot_with_gaps(axes[2], x, enu_arr[:, 2], gap_seconds=120.0, color='C2', label=f"U (RMS={stats['rmsU']:.3f})")
        for ax, lab in zip(axes, labels):
            ax.set_ylabel(lab)
            ax.grid(True, linestyle='--', alpha=0.4)
            # 加粗零刻度参考线
            ax.axhline(0.0, color='k', linewidth=1.2, alpha=0.7, zorder=1)
            ax.legend(loc='upper right', fontsize=9)
        # Set 2-sigma limits per axis
        e_lim = (-ylim_en, ylim_en) if ylim_en is not None else _ylims_robust(enu_arr[:, 0])
        n_lim = (-ylim_en, ylim_en) if ylim_en is not None else _ylims_robust(enu_arr[:, 1])
        u_lim = (-ylim_u, ylim_u) if ylim_u is not None else _ylims_robust(enu_arr[:, 2])
        axes[0].set_ylim(*e_lim)
        axes[1].set_ylim(*n_lim)
        axes[2].set_ylim(*u_lim)
    # 如需在单文件图上标注收敛，可在此处恢复标注逻辑（当前仅做统计输出，不绘制标注）
        axes[-1].set_xlabel(series.time_label)
        fig.suptitle(f'Positioning Error (ENU) - {name}')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        safe_name = os.path.splitext(name)[0]
        out_path = os.path.join(out_dir, f"{safe_name}_enu.png")
        plt.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved figure to {out_path}")
   

def main():
    ap = argparse.ArgumentParser(description='Plot ENU error curves from .pos/.flt files relative to true coordinates')
    ap.add_argument('input', help='Input .pos/.flt file or a directory containing .pos/.flt files')
    ap.add_argument('--true', default='true_crd.true_crd', help='Path to true coordinates file (default: %(default)s)')
    ap.add_argument('--site', default=None, help='Station/site code (override detection). E.g. HERS')
    # Output options
    ap.add_argument('--save', default=None, help='[Deprecated alias of --save-combined] Save combined E/N/U figure to PNG')
    ap.add_argument('--save-combined', default=None, help='Save combined E/N/U figure to PNG')
    ap.add_argument('--save-faceted', default=None, help='Save faceted figure (each file in its own subplot) to PNG')
    ap.add_argument('--save-individual-dir', default=None, help='Directory to save individual figures for each file')
    ap.add_argument('--max-cols', type=int, default=3, help='Max columns for faceted layout (default: %(default)s)')
    # 收敛/尺度相关选项（默认 0.1 m、5 历元、3D）
    ap.add_argument('--conv-thresh', type=float, default=0.1, help='Convergence threshold in meters (default: %(default)s)')
    ap.add_argument('--conv-win', type=int, default=10, help='Consecutive epochs required for convergence (default: %(default)s)')
    ap.add_argument('--conv-mode', choices=['3d', 'hz'], default='3d', help='Convergence error mode: 3d or hz (default: %(default)s)')
    ap.add_argument('--ylim', nargs=2, type=float, default=None, metavar=('EN_LIMIT', 'U_LIMIT'), help='Custom symmetric y-axis limits: EN_LIMIT for East/North, U_LIMIT for Up (default: auto 2-sigma)')
    args = ap.parse_args()

    true_file = args.true
    true_tbl = load_true_coords(true_file)
    if not true_tbl:
        print(f"No valid true coordinates parsed from {true_file}", file=sys.stderr)
        if not os.path.exists(true_file):
            print(f"-> True file does not exist: {true_file}", file=sys.stderr)
        else:
            try:
                with open(true_file, 'r', encoding='utf-8', errors='ignore') as _f:
                    head = ''.join([next(_f) for _ in range(5)])
                print("-> First lines of true file:\n" + head, file=sys.stderr)
            except Exception:
                print("-> Failed to read preview of true file.", file=sys.stderr)
        print("Please check the --true path or file format.", file=sys.stderr)
        sys.exit(2)

    data_paths = find_data_files(args.input)
    if not data_paths:
        print(f"No .pos/.flt files found under {args.input}", file=sys.stderr)
        if not os.path.exists(args.input):
            print(f"-> Input path does not exist: {args.input}", file=sys.stderr)
        elif os.path.isdir(args.input):
            # show a short listing to help debug
            try:
                files = os.listdir(args.input)
                sample = '\n'.join(files[:10])
                print(f"-> Directory exists. Sample entries:\n{sample}", file=sys.stderr)
            except Exception:
                print("-> Cannot list directory contents.", file=sys.stderr)
        else:
            print(f"-> Input exists but is not a .pos/.flt file (suffix check).", file=sys.stderr)
        print("Ensure the input path points to a .pos/.flt file or a directory containing such files.", file=sys.stderr)
        sys.exit(2)

    enu_all: List[Tuple[str, PosSeries, np.ndarray, Dict[str, float]]] = []

    for p in data_paths:
        site = args.site or infer_site_from_filename(p, true_tbl) or infer_site_from_header(p, true_tbl)
        use_std = False
        if site is None:
            print(f"WARNING: Cannot infer site for {p}. Will compute STD series instead.", file=sys.stderr)
            use_std = True
            tcoord = None
        else:
            tcoord = true_tbl.get(site)
            if tcoord is None:
                print(f"WARNING: Site {site} not found in true coordinates for {p}. Will compute STD series instead.", file=sys.stderr)
                use_std = True
        try:
            if p.lower().endswith('.flt'):
                series = parse_flt_file(p)
            else:
                series = parse_pos_file(p)
        except Exception as e:
            print(f"WARNING: Failed to parse {p}: {e}", file=sys.stderr)
            continue

        if use_std:
            # compute reference as the mean of the series and convert to geodetic for ENU
            ref_ecef = tuple(np.mean(series.ecef, axis=0).tolist())
            lat_deg, lon_deg, hgt_m = ecef_to_geodetic(ref_ecef[0], ref_ecef[1], ref_ecef[2])
            enu, stats = compute_enu_error(series, ref_ecef, (lat_deg, lon_deg, hgt_m))
            # When using STD reference (no true value), report standard deviation instead of RMS
            enu = np.asarray(enu, dtype=float)
            stats['rmsE'] = float(np.std(enu[:, 0]))
            stats['rmsN'] = float(np.std(enu[:, 1]))
            stats['rmsU'] = float(np.std(enu[:, 2]))
            idx_conv, t_conv = compute_convergence_time(enu, series, thresh=args.conv_thresh, win=args.conv_win, mode=args.conv_mode)
            name = os.path.basename(p) + ' (STD)'
        else:
            enu, stats = compute_enu_error(series, (tcoord.x, tcoord.y, tcoord.z), (tcoord.lat_deg, tcoord.lon_deg, tcoord.hgt_m))
            idx_conv, t_conv = compute_convergence_time(enu, series, thresh=args.conv_thresh, win=args.conv_win, mode=args.conv_mode)
            name = os.path.basename(p)

        enu_all.append((name, series, enu, stats))
        if t_conv is not None:
            if isinstance(t_conv, (int, float)):
                t_str = f"{t_conv:.1f}"
            elif pd is not None:
                t_str = pd.to_datetime(t_conv).strftime('%H:%M:%S')
            else:
                t_str = str(t_conv)
        else:
            t_str = 'N/A'
        metric_label = 'STD' if use_std else 'RMS'
        print(f"Parsed {name}: N={enu.shape[0]} | {metric_label}(E,N,U) = {stats['rmsE']:.3f}, {stats['rmsN']:.3f}, {stats['rmsU']:.3f} m | Conv(win={args.conv_win}, thresh={args.conv_thresh}m, mode={args.conv_mode}) idx={idx_conv if idx_conv is not None else 'N/A'}, time={t_str}")

    if not enu_all:
        print("No valid series to plot.", file=sys.stderr)
        print(f"-> Scanned {len(data_paths)} data paths, but all were skipped due to missing site mapping or parse errors.", file=sys.stderr)
        print("-> Tips: use --site to force a site code, or check that your files contain recognizable station codes/headers.", file=sys.stderr)
        sys.exit(2)

    # Determine outputs
    save_combined = args.save_combined or args.save  # backward compat
    save_faceted = args.save_faceted
    save_individual_dir = args.save_individual_dir

    # Determine base directory for plot outputs
    if os.path.isdir(args.input):
        base_dir = args.input  # Input is a directory
    else:
        base_dir = os.path.dirname(args.input) or '.'  # Input is a file

    # Create plot directory within the base directory
    plot_dir = os.path.join(base_dir, 'plot')
    os.makedirs(plot_dir, exist_ok=True)

    # Default save paths for combined, faceted, and individual plots
    if not save_combined:
        save_combined = os.path.join(plot_dir, 'enu_combined.png')
    if not save_faceted:
        save_faceted = os.path.join(plot_dir, 'enu_faceted.png')
    if not save_individual_dir:
        save_individual_dir = plot_dir  # Default individual plots to the same plot directory

    # Extract custom y-limits
    ylim_en = args.ylim[0] if args.ylim else None
    ylim_u = args.ylim[1] if args.ylim else None

    # Generate plots
    if save_combined:
        plot_enu_series(enu_all, save_path=save_combined, ylim_en=ylim_en, ylim_u=ylim_u)
    if save_faceted:
        plot_enu_faceted(enu_all, save_path=save_faceted, max_cols=args.max_cols, ylim_en=ylim_en, ylim_u=ylim_u)
    if save_individual_dir:
        plot_enu_individual(enu_all, out_dir=save_individual_dir, ylim_en=ylim_en, ylim_u=ylim_u)


if __name__ == '__main__':
    main()
