#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
并行调用 GREAT_PVT PPP 解算程序处理文件夹中的所有观测文件。

用法:
  python run_greatpvt_parallel.py <XML配置文件> <数据文件夹> [选项]

选项:
  -j, --jobs      并行进程数（默认: CPU 核心数）
  --exe           GREAT_PVT 可执行文件路径（默认: ./build_Linux/Bin/GREAT_PVT）
  --dry-run       仅列出将要处理的文件，不实际执行
  -v, --verbose   显示详细输出

示例:
  python run_greatpvt_parallel.py data/GA_2023070/xml/GREAT_PPPFLT_static_DF_Fixed_AUG.xml data/GA_2023070/svr_30/ -j 4
"""

from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# ---------------------------------------------------------------------------
# Process lifecycle management (no pkill / no process-name lookup)
# ---------------------------------------------------------------------------

_shutdown_event = threading.Event()
_active_lock = threading.Lock()
_active_processes: list[subprocess.Popen] = []


def _start_process(cmd: list[str], capture: bool = True) -> subprocess.Popen:
    """Start a subprocess in a new process group so we can terminate it reliably."""
    kwargs: dict = {
        "stdout": subprocess.PIPE if capture else None,
        "stderr": subprocess.PIPE if capture else None,
        "text": True,
    }
    if os.name == "nt":
        kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        kwargs["preexec_fn"] = os.setsid

    proc = subprocess.Popen(cmd, **kwargs)

    with _active_lock:
        if not _shutdown_event.is_set():
            _active_processes.append(proc)
        else:
            # Shutdown was requested before we could register the process
            _kill_process(proc)
            raise RuntimeError("Shutdown requested before process started")
    return proc


def _kill_process(proc: subprocess.Popen) -> None:
    """Kill a single subprocess (and its process group)."""
    if proc.poll() is not None:
        return
    try:
        if os.name == "nt":
            proc.send_signal(signal.CTRL_BREAK_EVENT)
        else:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except Exception:
        pass


def _terminate_all() -> None:
    """Terminate every subprocess we are currently tracking."""
    _shutdown_event.set()
    with _active_lock:
        procs = list(_active_processes)

    # Phase 1: polite SIGTERM to the whole process group
    for proc in procs:
        _kill_process(proc)

    # Phase 2: wait a moment, then SIGKILL any survivors
    time.sleep(2)
    with _active_lock:
        procs = list(_active_processes)
    for proc in procs:
        if proc.poll() is None:
            try:
                proc.kill()
                proc.wait(timeout=5)
            except Exception:
                pass


def _signal_handler(signum: int, _frame) -> None:
    print(f"\n[run_greatpvt_parallel] 收到信号 {signum}，正在终止所有 GREAT_PVT 子进程...")
    _terminate_all()
    sys.exit(128 + signum)


signal.signal(signal.SIGINT, _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_basepath_from_xml(xml_file: str) -> str | None:
    """从 XML 配置文件中提取 <inputs><basepath> 的值。"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        inp = root.find("inputs")
        if inp is not None:
            bp = inp.findtext("basepath", default="")
            if bp:
                return bp.strip()
    except Exception:
        pass
    return None


def _extract_site_4ch(obs_file: str) -> str:
    """从观测文件名提取 4 字符站名（与 GREAT_PVT -o 逻辑一致）。"""
    fname = Path(obs_file).name
    site = fname[:min(4, len(fname))]
    return site.upper()


def _check_flt_result(obs_file: str, basepath: str | None, start_time: float) -> tuple[bool, str]:
    """检查与观测文件对应的 flt 结果文件是否有有效数据行（只考虑 start_time 之后修改的文件）。"""
    site = _extract_site_4ch(obs_file)

    # 确定搜索根目录
    search_roots = []
    if basepath:
        bp = Path(basepath)
        if bp.is_absolute():
            search_roots.append(bp)
        else:
            # 相对路径：尝试相对于当前工作目录
            search_roots.append(Path.cwd() / bp)
    # 同时也在当前工作目录搜索
    search_roots.append(Path.cwd())

    # 去重并保持顺序
    seen = set()
    unique_roots = []
    for r in search_roots:
        resolved = r.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_roots.append(r)

    # 搜索 result_* 子目录下的 flt 文件
    flt_files = []
    for root in unique_roots:
        if not root.exists():
            continue
        for subdir in root.iterdir():
            if subdir.is_dir() and subdir.name.startswith("result"):
                for f in subdir.iterdir():
                    if f.is_file() and f.suffix == ".flt" and f.name.startswith(site):
                        try:
                            if f.stat().st_mtime >= start_time:
                                flt_files.append(f)
                        except Exception:
                            pass
        # 也直接在根目录搜索（fallback）
        for f in root.iterdir():
            if f.is_file() and f.suffix == ".flt" and f.name.startswith(site):
                try:
                    if f.stat().st_mtime >= start_time:
                        flt_files.append(f)
                except Exception:
                    pass

    if not flt_files:
        return False, "未找到 .flt 结果文件"

    for flt_file in flt_files:
        try:
            with open(flt_file, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as exc:
            return False, f"读取 .flt 失败: {exc}"

        data_lines = 0
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines += 1

        if data_lines > 0:
            return True, f"结果正常 ({data_lines} 行数据) [{flt_file.name}]"

    return False, ".flt 文件无有效数据"


# ---------------------------------------------------------------------------
# File scanning (replicated from run_uduc_parallel.py)
# ---------------------------------------------------------------------------

def is_obs_file(fname: str) -> bool:
    """
    复现 GREAT_PVT/RINEX 观测文件判断逻辑：
      - 不是 .obs 格式（那是广播星历）
      - 或 RINEX2 观测文件: *.[0-9][0-9]o / *.[0-9][0-9]O
      - 或 RINEX3/4 观测文件: 扩展名包含 "rnx" / "RNX"（但排除含 "_MN." 的导航文件）
    """
    if len(fname) <= 4:
        return False

    ext = Path(fname).suffix
    if not ext:
        return False
    ext = ext[1:]  # 去掉点
    ns = len(ext)
    if ns < 3 or ns > 8:
        return False

    ext_lower = ext.lower()
    fname_lower = fname.lower()

    # b1: 不是 .obs
    b1 = ext_lower != "obs"

    # b2 && b3 && b4: RINEX2 观测文件, 如 21o, 21O
    b2 = ext.endswith("o") or ext.endswith("O")
    b3 = "0" <= ext[ns - 2] <= "9"
    b4 = "0" <= ext[ns - 3] <= "9"

    # b5: 扩展名包含 rnx，且文件名不含 _MN.
    b5_0 = "rnx" in ext_lower
    b5_2 = "_mn." in fname_lower
    b5 = b5_0 and (not b5_2)

    isobs = (b1 and (b2 and b3 and b4)) or b5
    return isobs


def find_obs_files(data_dir: str) -> list[str]:
    """扫描文件夹，返回所有观测文件的绝对路径。"""
    data_path = Path(data_dir).resolve()
    if not data_path.is_dir():
        raise ValueError(f"数据文件夹不存在: {data_dir}")

    obs_files = []
    for fpath in data_path.iterdir():
        if fpath.is_file() and is_obs_file(fpath.name):
            obs_files.append(str(fpath))

    # 按文件名排序，保证顺序一致
    obs_files.sort()
    return obs_files


# ---------------------------------------------------------------------------
# Run single GREAT_PVT instance
# ---------------------------------------------------------------------------

def run_single(obs_file: str, xml_file: str, exe_path: str,
               basepath: str | None, verbose: bool = False) -> tuple[str, bool, str]:
    """Run GREAT_PVT for one observation file.  Returns (obs_file, success, message)."""
    if _shutdown_event.is_set():
        return obs_file, False, "已取消"

    site = _extract_site_4ch(obs_file)
    cmd = [exe_path, "-x", xml_file, "-o", obs_file]
    start = time.time()

    try:
        proc = _start_process(cmd, capture=not verbose)
    except RuntimeError:
        return obs_file, False, "已取消"

    try:
        if verbose:
            ret = proc.wait()
            stdout = stderr = ""
        else:
            stdout, stderr = proc.communicate()
            ret = proc.returncode
    finally:
        with _active_lock:
            if proc in _active_processes:
                _active_processes.remove(proc)

    if _shutdown_event.is_set():
        return obs_file, False, "已取消"

    elapsed = time.time() - start
    if ret == 0:
        has_flt, flt_msg = _check_flt_result(obs_file, basepath, start)
        if has_flt:
            return obs_file, True, f"成功 ({elapsed:.1f}s) {flt_msg}"
        else:
            return obs_file, False, f"解算失败: {flt_msg}"
    else:
        stderr = stderr.strip()[:200] if stderr else ""
        return obs_file, False, f"失败 (code={ret}) {stderr}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="并行调用 GREAT_PVT PPP 解算程序",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("xml", help="GREAT_PVT XML 配置文件路径")
    parser.add_argument("data_dir", help="包含观测文件的文件夹路径")
    parser.add_argument(
        "-j", "--jobs", type=int, default=os.cpu_count(),
        help="并行进程数（默认: CPU 核心数 = %(default)s）"
    )
    parser.add_argument(
        "--exe", default="./build_Linux/Bin/GREAT_PVT",
        help="GREAT_PVT 可执行文件路径（默认: %(default)s）"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="仅列出将要处理的文件，不实际执行"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="显示子进程的标准输出（默认只捕获 stderr）"
    )

    args = parser.parse_args()

    xml_file = Path(args.xml).resolve()
    if not xml_file.is_file():
        print(f"错误: XML 配置文件不存在: {args.xml}")
        sys.exit(1)

    exe_path = Path(args.exe).resolve()
    if not exe_path.is_file():
        print(f"错误: 可执行文件不存在: {args.exe}")
        sys.exit(1)

    try:
        obs_files = find_obs_files(args.data_dir)
    except ValueError as e:
        print(f"错误: {e}")
        sys.exit(1)

    if not obs_files:
        print(f"警告: 在 {args.data_dir} 中未找到任何观测文件")
        sys.exit(0)

    # 从 XML 提取 basepath，用于后续结果文件定位
    basepath = _get_basepath_from_xml(str(xml_file))

    print(f"XML 配置 : {xml_file}")
    print(f"数据目录 : {Path(args.data_dir).resolve()}")
    print(f"可执行文件: {exe_path}")
    print(f"观测文件数: {len(obs_files)}")
    print(f"并行进程数: {args.jobs}")
    if basepath:
        print(f"输出根目录: {basepath}")
    print("-" * 60)

    if args.dry_run:
        print("[Dry-run] 将要处理的观测文件:")
        for f in obs_files:
            site = _extract_site_4ch(f)
            print(f"  - {Path(f).name}  (站名: {site})")
        sys.exit(0)

    # 记录总耗时
    t0 = time.time()
    success_count = 0
    fail_count = 0

    # Use ThreadPoolExecutor so that all subprocesses are children of the main
    # process and can be reliably tracked and terminated via Popen objects.
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_obs = {
            executor.submit(run_single, obs, str(xml_file), str(exe_path), basepath, args.verbose): obs
            for obs in obs_files
        }

        for future in as_completed(future_to_obs):
            try:
                obs_file, success, msg = future.result()
            except Exception as exc:
                obs_file = future_to_obs[future]
                success = False
                msg = f"异常: {exc}"
            name = Path(obs_file).name
            status = "[OK]" if success else "[FAIL]"
            print(f"{status} {name:<40s} {msg}")
            if success:
                success_count += 1
            else:
                fail_count += 1

    total_time = time.time() - t0
    print("-" * 60)
    print(f"全部完成: 成功 {success_count} / 失败 {fail_count} / 总计 {len(obs_files)}")
    print(f"总耗时  : {total_time:.1f}s")


if __name__ == "__main__":
    main()
