#!/usr/bin/env python3
# coding: utf-8
"""
fit2crd.py

从输入文件夹中提取每个站的 PPP-AR（或 float）结果文件，
取最后 20 个历元计算平均 XYZ 坐标，转换为 lat/lon/hgt 后
更新到 true_crd.true_crd 文件中。

同时在终端输出每个站的真值坐标（若已有）及最后 20 历元的标准差（std）。

用法:
    python fit2crd.py <result_dir> [true_crd_file]

示例:
    python fit2crd.py /home/lzx/code/GREAT-PVT/data/GA_2023070/svr_30_obs/result_EST_UC_GRE
"""

import os
import sys
import math
import re
import numpy as np

# ------------------------------------------------------------------
# WGS84 及坐标转换（与 GTest_analysis_cvg_cst.py 保持一致）
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
    return math.degrees(lat), math.degrees(lon), h


def ecef_to_enu_matrix(lat_deg, lon_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    return np.array([
        [-sin_lon, cos_lon, 0.0],
        [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
        [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat]
    ])


def compute_enu(ecef, ref_ecef, ref_llh):
    """计算 ENU 误差向量"""
    R = ecef_to_enu_matrix(ref_llh[0], ref_llh[1])
    diff = np.array(ecef, dtype=float) - np.array(ref_ecef, dtype=float)
    return np.dot(R, diff)


# ------------------------------------------------------------------
# true_crd 文件读写
# ------------------------------------------------------------------
def load_true_crd(path):
    """加载 true_crd.true_crd，返回 (站点坐标表, 原始行列表)"""
    tbl = {}
    lines = []
    if not os.path.isfile(path):
        return tbl, lines
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            lines.append(line)
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
    return tbl, lines


def write_true_crd(path, records):
    """
    重写 true_crd 文件。
    records: dict {site: (x, y, z, lat, lon, hgt, epnum)}
    """
    header = [
        "%%  The site coordinate summary file\n",
        "%%  program  :  fit2crd.py\n",
        "%%  comment  :  the first four cols are mandatory, the other cols are optional.\n",
        "%%  ---------------------------------------------------------------------------------------------------------\n",
        "%%  site             X/m             Y/m             Z/m         lat/deg         lon/deg      hgt/m    epnum\n",
        "%%  ---------------------------------------------------------------------------------------------------------\n",
        "%%\n",
    ]
    with open(path, 'w') as f:
        for h in header:
            f.write(h)
        for site in sorted(records.keys()):
            r = records[site]
            f.write(
                "    %-4s %14.3f %15.3f %15.3f %15.9f %15.9f %10.3f %8d    0.000    0.000    0.000    0.000\n"
                % (site, r[0], r[1], r[2], r[3], r[4], r[5], r[6])
            )


# ------------------------------------------------------------------
# 结果文件查找与解析
# ------------------------------------------------------------------
def extract_site_name(filename):
    """从文件名中提取 4 字符站名"""
    base = os.path.splitext(os.path.basename(filename))[0]

    # 模式1: 4个大写字母+年积日数字，如 CBRA2023070、GODN2023305
    m = re.search(r'([A-Z]{4})\d{4,7}', base)
    if m:
        return m.group(1)

    # 模式2: 开头4个字母，如 GODN-PPP_sta_DF_Fixed
    if len(base) >= 4 and base[:4].isalpha():
        return base[:4].upper()

    # 模式3: 任意位置连续4个大写字母
    m = re.search(r'[A-Z]{4}', base)
    if m:
        return m.group(0)

    return base[:4].upper()


def get_fixed_ratio(fpath):
    """读取 flt 文件，返回 Fixed 历元占比"""
    total = 0
    fixed = 0
    with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line.strip().startswith('#'):
                continue
            ss = line.split()
            if len(ss) < 17:
                continue
            total += 1
            # AmbStatus 列索引约为 16（0-based）
            if ss[16].lower() == 'fixed':
                fixed += 1
    return fixed / total if total > 0 else 0.0


def find_result_files(dir_path):
    """
    找出每个站的结果文件，优先级：pppar/fixed > float。
    返回 dict {site: filename}
    """
    files = [f for f in os.listdir(dir_path) if f.endswith('.flt')]

    pppar = {}
    float_ = {}
    unclassified = {}

    for f in files:
        site = extract_site_name(f)
        fname_lower = f.lower()
        # 明确标识为 PPP-AR / Fixed
        if '_pppar' in fname_lower or 'fixed' in fname_lower:
            pppar[site] = f
        # 明确标识为 float
        elif '_float' in fname_lower or fname_lower.endswith('_float.flt'):
            float_[site] = f
        else:
            unclassified[site] = f

    # 对未分类文件，按 Fixed 比例判断（>50% 视为 AR）
    for site, f in unclassified.items():
        if site in pppar:
            continue
        fpath = os.path.join(dir_path, f)
        ratio = get_fixed_ratio(fpath)
        if ratio > 0.5:
            pppar[site] = f
        else:
            float_[site] = f

    result = {}
    for site in sorted(set(list(pppar.keys()) + list(float_.keys()))):
        if site in pppar:
            result[site] = pppar[site]
        else:
            result[site] = float_[site]
    return result


def read_last_n_epochs(fpath, n=20):
    """读取 flt 文件最后 n 个有效历元，返回 [(sow, x, y, z), ...]"""
    epochs = []
    with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line.strip().startswith('#'):
                continue
            ss = line.split()
            if len(ss) < 19:
                # float 文件可能只有 18 列，但前 4 列（sow, x, y, z）一定有
                if len(ss) < 4:
                    continue
            try:
                sow = float(ss[0])
                x = float(ss[1])
                y = float(ss[2])
                z = float(ss[3])
                epochs.append((sow, x, y, z))
            except ValueError:
                continue
    return epochs[-n:]


def std_sample(vals):
    """样本标准差 (N-1)"""
    n = len(vals)
    if n < 2:
        return 0.0
    m = sum(vals) / n
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (n - 1))


# ------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python fit2crd.py <result_dir> [true_crd_file]")
        sys.exit(1)

    result_dir = sys.argv[1]
    if len(sys.argv) > 2:
        true_crd_path = os.path.abspath(sys.argv[2])
    else:
        # 默认与 script 同级目录下的 true_crd.true_crd
        script_dir = os.path.dirname(os.path.abspath(__file__))
        true_crd_path = os.path.join(script_dir, '..', 'true_crd.true_crd')
        true_crd_path = os.path.abspath(true_crd_path)

    if not os.path.isdir(result_dir):
        print(f"Error: directory not found: {result_dir}")
        sys.exit(1)

    # 加载已有真值坐标
    true_tbl, _ = load_true_crd(true_crd_path)

    # 查找结果文件
    site_files = find_result_files(result_dir)
    if not site_files:
        print(f"No .flt result files found in {result_dir}")
        sys.exit(1)

    # 准备输出记录（保留已有站点）
    output_records = {}
    for site, coord in true_tbl.items():
        output_records[site] = (*coord, 2880)  # 默认 epnum

    print(f"Found {len(site_files)} station(s) in [{result_dir}]")
    print("-" * 110)
    print(f"{'Site':>6}  {'Mode':>8}  {'EpNum':>6}  "
          f"{'X_true/m':>14} {'Y_true/m':>14} {'Z_true/m':>14}  "
          f"{'std_E/m':>10} {'std_N/m':>10} {'std_U/m':>10}")
    print("-" * 110)

    for site in sorted(site_files.keys()):
        fpath = os.path.join(result_dir, site_files[site])
        epochs = read_last_n_epochs(fpath, n=20)
        n_ep = len(epochs)

        if n_ep == 0:
            print(f"  {site:4}  {'--':>8}  {'0':>6}  {'--':>14} {'--':>14} {'--':>14}  {'--':>10} {'--':>10} {'--':>10}")
            continue

        xs = [e[1] for e in epochs]
        ys = [e[2] for e in epochs]
        zs = [e[3] for e in epochs]

        mean_x = sum(xs) / n_ep
        mean_y = sum(ys) / n_ep
        mean_z = sum(zs) / n_ep

        lat, lon, hgt = ecef_to_geodetic(mean_x, mean_y, mean_z)
        output_records[site] = (mean_x, mean_y, mean_z, lat, lon, hgt, n_ep)

        # 判断模式（仅用于显示）
        fname_lower = site_files[site].lower()
        mode_str = 'pppar' if ('_pppar' in fname_lower or 'fixed' in fname_lower) else 'float'

        # 标准差计算与输出
        if site in true_tbl:
            true_ecef = true_tbl[site][:3]
            true_llh = true_tbl[site][3:6]
            enus = [compute_enu((e[1], e[2], e[3]), true_ecef, true_llh) for e in epochs]
            es = [e[0] for e in enus]
            ns = [e[1] for e in enus]
            us = [e[2] for e in enus]

            std_e = std_sample(es)
            std_n = std_sample(ns)
            std_u = std_sample(us)

            print(f"  {site:4}  {mode_str:>8}  {n_ep:>6}  "
                  f"{true_ecef[0]:14.4f} {true_ecef[1]:14.4f} {true_ecef[2]:14.4f}  "
                  f"{std_e:10.4f} {std_n:10.4f} {std_u:10.4f}")
        else:
            # 无真值时，输出 XYZ std
            std_x = std_sample(xs)
            std_y = std_sample(ys)
            std_z = std_sample(zs)
            print(f"  {site:4}  {mode_str:>8}  {n_ep:>6}  "
                  f"{mean_x:14.4f} {mean_y:14.4f} {mean_z:14.4f}  "
                  f"{std_x:10.4f} {std_y:10.4f} {std_z:10.4f}  (XYZ std, no true coord)")

    print("-" * 110)

    # 写入/更新 true_crd.true_crd
    write_true_crd(true_crd_path, output_records)
    print(f"Updated [{true_crd_path}] with {len(output_records)} station(s).")


if __name__ == '__main__':
    main()
