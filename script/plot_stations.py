#!/usr/bin/env python3
"""Parse station files for APPROX POSITION XYZ and OBSERVER, then plot locations.

Usage:
    python3 plot_stations.py /path/to/files/*.25d -o stations.png
    python3 plot_stations.py /path/to/files/*.25d --sitename
    python3 plot_stations.py --svr_usr /path/to/svr_files/*.25d /path/to/usr_files/*.25d

The script scans the provided files, extracts the OBSERVER string and the
APPROX POSITION XYZ (ECEF meters), converts to geodetic lat/lon, and plots
stations colored by OBSERVER.
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# 尝试使用共享的 plot_utils 来设置中文字体和负号显示
try:
    from script.plot_utils import setup_chinese_font
except Exception:
    try:
        from plot_utils import setup_chinese_font
    except Exception:
        setup_chinese_font = None

if setup_chinese_font:
    setup_chinese_font()

NUM_RE = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')


def extract_observer_and_xyz(path):
    observer = None
    approx = None
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            up = line.upper()
            if 'OBSERVER' in up:
                # OBSERVER / AGENCY usually at left columns
                left = line[:40].strip()
                if left:
                    observer = left.split()[0]
            if 'APPROX POSITION XYZ' in up or 'ACCURATE POSITION XYZ' in up or 'POSITION XYZ' in up:
                nums = NUM_RE.findall(line)
                if len(nums) >= 3:
                    try:
                        approx = (float(nums[0]), float(nums[1]), float(nums[2]))
                    except ValueError:
                        approx = None
            if 'END OF HEADER' in up:
                break
    if observer is None:
        observer = os.path.splitext(os.path.basename(path))[0]
    return observer, approx


def ecef_to_geodetic(x, y, z):
    # WGS84
    a = 6378137.0
    f = 1.0 / 298.257223563
    b = a * (1 - f)
    e2 = f * (2 - f)
    ep2 = (a * a - b * b) / (b * b)

    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    # Bowring's method
    theta = math.atan2(z * a, p * b)
    st = math.sin(theta)
    ct = math.cos(theta)
    lat = math.atan2(z + ep2 * b * st * st * st, p - e2 * a * ct * ct * ct)
    N = a / math.sqrt(1 - e2 * math.sin(lat) * math.sin(lat))
    alt = p / math.cos(lat) - N
    return math.degrees(lat), math.degrees(lon), alt


def gather(files):
    data = []  # list of (observer, lat, lon, sitename, path)
    for path in files:
        obs, approx = extract_observer_and_xyz(path)
        if approx is None:
            continue
        lat, lon, _ = ecef_to_geodetic(*approx)
        sitename = os.path.basename(path)[:4]
        data.append((obs, lat, lon, sitename, path))
    return data


def resolve_input_paths(path):
    if os.path.isfile(path):
        return [path]
    if os.path.isdir(path):
        resolved = []
        for root, _, files in os.walk(path):
            for filename in files:
                if filename.startswith('.'):
                    continue
                resolved.append(os.path.join(root, filename))
        return sorted(resolved)
    return sorted(glob.glob(path))


def load_station_names(path):
    names = set()
    for resolved_path in resolve_input_paths(path):
        if not os.path.isfile(resolved_path):
            continue
        with open(resolved_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                text = line.strip()
                if not text or text.startswith('#') or text.startswith('='):
                    continue
                parts = text.split()
                if not parts:
                    continue
                name = parts[0]
                if name.upper() in {'SERVICE', 'USER'} and len(parts) >= 2:
                    name = parts[1]
                names.add(name[:4])
    return names


def _median_step(values):
    if len(values) <= 1:
        return 0.0
    diffs = [b - a for a, b in zip(values[:-1], values[1:]) if b > a]
    if not diffs:
        return 0.0
    diffs.sort()
    mid = len(diffs) // 2
    if len(diffs) % 2:
        return diffs[mid]
    return (diffs[mid - 1] + diffs[mid]) / 2.0


def _to_float(text):
    if text is None:
        return None
    try:
        value = str(text).strip()
        if not value:
            return None
        return float(value)
    except Exception:
        return None


def load_grid_csv(path):
    """Load a CLAS-style grid CSV.

    The file format is expected to contain columns like:
    Compact Network ID, GRID No., Latitude, Longitude
    The first column may be blank for rows after the first row in each network.
    """
    if not path:
        return []
    if not os.path.exists(path):
        raise FileNotFoundError(f'grid csv not found: {path}')

    with open(path, 'r', encoding='utf-8-sig', newline='') as fh:
        reader = csv.reader(fh)
        rows = list(reader)

    if not rows:
        return []

    header = rows[0]
    if len(header) < 4:
        raise ValueError('grid csv must have at least 4 columns')

    grid_cells = []
    current_net_id = None
    current_points = []

    for row in rows[1:]:
        if not row:
            continue
        row = row + [''] * max(0, 4 - len(row))
        net_text = row[0].strip() if len(row) > 0 else ''
        lat = _to_float(row[2])
        lon = _to_float(row[3])
        if lat is None or lon is None:
            continue

        if net_text:
            if current_net_id is not None and current_points:
                grid_cells.append(_finalize_grid_cell(current_net_id, current_points))
            current_net_id = net_text
            current_points = []

        if current_net_id is None:
            current_net_id = 'Network_0'
        current_points.append((lat, lon))

    if current_net_id is not None and current_points:
        grid_cells.append(_finalize_grid_cell(current_net_id, current_points))

    return grid_cells


def _finalize_grid_cell(network_id, points):
    lats = sorted({float(lat) for lat, _ in points})
    lons = sorted({float(lon) for _, lon in points})
    d_lat = _median_step(lats)
    d_lon = _median_step(lons)
    lat_min = min(lat for lat, _ in points) - d_lat / 2.0
    lat_max = max(lat for lat, _ in points) + d_lat / 2.0
    lon_min = min(lon for _, lon in points) - d_lon / 2.0
    lon_max = max(lon for _, lon in points) + d_lon / 2.0
    return {
        'network_id': network_id,
        'points': points,
        'bounds': {
            'lat_range': (lat_min, lat_max),
            'lon_range': (lon_min, lon_max),
        },
    }


def _annotate_points(ax, points, transform=None, text_offset=(0, 7), text_color='black', fontsize=7, stagger=False):
    for index, (lon, lat, sitename) in enumerate(points):
        dx, dy = text_offset
        if stagger:
            dx += 2 if index % 2 == 0 else -2
            dy += (index % 3 - 1) * 2
        kwargs = {
            'xy': (lon, lat),
            'xytext': (dx, dy),
            'textcoords': 'offset points',
            'ha': 'center',
            'va': 'bottom' if dy >= 0 else 'top',
            'fontsize': fontsize,
            'color': text_color,
            'bbox': dict(boxstyle='round,pad=0.10', facecolor='white', alpha=0.7, edgecolor='none'),
            'zorder': 10,
        }
        if transform is not None:
            kwargs['transform'] = transform
        ax.annotate(sitename, **kwargs)


def _add_service_legend(ax, service_count, user_count):
    legend_handles = [
        Line2D([0], [0], marker='^', linestyle='None', markersize=4,
               markerfacecolor='#e74c3c', markeredgecolor='black', label=f'svr({service_count})'),
        Line2D([0], [0], marker='o', linestyle='None', markersize=3.5,
               markerfacecolor='#3498db', markeredgecolor='black', label=f'usr({user_count})'),
    ]
    ax.legend(handles=legend_handles, loc='upper right',
          frameon=True, framealpha=0.85, fontsize=7, borderpad=0.25,
          handletextpad=0.45, labelspacing=0.3, borderaxespad=0.4)


def _draw_grid_overlay(ax, cells, draw_points=False):
    boundary_color = '#4d4d4d'
    point_colors = ['#ff8c00', '#2ca02c', '#7f3c8d', '#8c564b', '#17a2b8', '#c49c94', '#6f42c1', '#bcbd22']
    for idx, cell in enumerate(cells):
        color = point_colors[idx % len(point_colors)]
        lat_min, lat_max = cell['bounds']['lat_range']
        lon_min, lon_max = cell['bounds']['lon_range']
        xs = [lon_min, lon_max, lon_max, lon_min, lon_min]
        ys = [lat_min, lat_min, lat_max, lat_max, lat_min]
        ax.plot(xs, ys, linestyle='--', linewidth=0.8, color=boundary_color, alpha=0.85, zorder=2)
        if draw_points:
            pts = cell.get('points', [])
            if pts:
                lats = [p[0] for p in pts]
                lons = [p[1] for p in pts]
                ax.scatter(lons, lats, s=10, color=color, edgecolors='black', linewidths=0.2, zorder=3)


def _plot_map_base(ax, service_mode, default_points, service_points, user_points,
                   grid_cells, show_grid_points, show_sitename, transform=None):
    """Common plotting helper for cartopy and plain matplotlib axes."""
    if grid_cells:
        _draw_grid_overlay(ax, grid_cells, draw_points=show_grid_points)

    if service_mode:
        if service_points:
            serv_lons = [p[0] for p in service_points]
            serv_lats = [p[1] for p in service_points]
            ax.scatter(serv_lons, serv_lats, s=40, marker='^', color='#e74c3c',
                       edgecolors='black', linewidths=0.5, zorder=5, transform=transform)
            _annotate_points(ax, service_points, transform=transform, text_offset=(4, 4), text_color='#c0392b', fontsize=4, stagger=True)
        if user_points:
            user_lons = [p[0] for p in user_points]
            user_lats = [p[1] for p in user_points]
            ax.scatter(user_lons, user_lats, s=25, marker='o', color='#3498db',
                       edgecolors='black', linewidths=0.5, zorder=6, transform=transform)
            _annotate_points(ax, user_points, transform=transform, text_offset=(-4, -5), text_color='#1f618d', fontsize=4, stagger=True)
        if default_points:
            lons = [p[0] for p in default_points]
            lats = [p[1] for p in default_points]
            ax.scatter(lons, lats, s=20, marker='o', color='#3498db',
                       edgecolors='black', linewidths=0.5, alpha=0.8, zorder=4, transform=transform)
            if show_sitename:
                _annotate_points(ax, default_points, transform=transform, text_offset=(0, 4), text_color='black', fontsize=4, stagger=False)
        _add_service_legend(ax, len(service_points), len(user_points))
    else:
        if default_points:
            lons = [p[0] for p in default_points]
            lats = [p[1] for p in default_points]
            ax.scatter(lons, lats, s=40, marker='o', color='#3498db',
                       edgecolors='black', linewidths=0.5, alpha=0.8, zorder=5, transform=transform)
            if show_sitename:
                _annotate_points(ax, default_points, transform=transform, text_offset=(0, 7), text_color='black')


def plot_by_observer(data, outfname=None, show_sitename=False, service_paths=None, user_paths=None, grid_cells=None, show_grid_points=False):
    service_mode = service_paths is not None or user_paths is not None
    service_paths = service_paths or set()
    user_paths = user_paths or set()

    default_points = []
    service_points = []
    user_points = []

    for obs, lat, lon, sitename, path in data:
        point = (lon, lat, sitename)
        if service_mode and path in service_paths:
            service_points.append(point)
        elif service_mode and path in user_paths:
            user_points.append(point)
        else:
            default_points.append(point)

    all_points = default_points + service_points + user_points
    station_lons = [p[0] for p in all_points]
    station_lats = [p[1] for p in all_points]
    grid_lons = []
    grid_lats = []
    if grid_cells:
        for cell in grid_cells:
            lat_min, lat_max = cell['bounds']['lat_range']
            lon_min, lon_max = cell['bounds']['lon_range']
            grid_lons.extend([lon_min, lon_max])
            grid_lats.extend([lat_min, lat_max])

    extent_lons = station_lons if station_lons else grid_lons
    extent_lats = station_lats if station_lats else grid_lats

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        fig = plt.figure(figsize=(10, 5), facecolor='white')
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_facecolor('white')
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, facecolor='white', edgecolor='gray', linewidth=0.5)
        ax.add_feature(cfeature.OCEAN, facecolor='white')
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.3, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        try:
            gl.xlabel_style = {'size': 8}
            gl.ylabel_style = {'size': 8}
        except Exception:
            pass

        if extent_lons and extent_lats:
            lon_range = max(extent_lons) - min(extent_lons)
            lat_range = max(extent_lats) - min(extent_lats)
            minlon = max(min(extent_lons) - lon_range * 0.15, -180)
            maxlon = min(max(extent_lons) + lon_range * 0.15, 180)
            minlat = max(min(extent_lats) - lat_range * 0.15, -90)
            maxlat = min(max(extent_lats) + lat_range * 0.15, 90)
            ax.set_extent([minlon, maxlon, minlat, maxlat], crs=ccrs.PlateCarree())

        _plot_map_base(ax, service_mode, default_points, service_points, user_points,
                       grid_cells, show_grid_points, show_sitename, transform=ccrs.PlateCarree())
    except Exception:
        fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')
        ax.set_facecolor('white')

        _plot_map_base(ax, service_mode, default_points, service_points, user_points,
                       grid_cells, show_grid_points, show_sitename, transform=None)

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.grid(True, linestyle='--', alpha=0.3, color='gray')
        if extent_lons and extent_lats:
            lon_range = max(extent_lons) - min(extent_lons)
            lat_range = max(extent_lats) - min(extent_lats)
            ax.set_xlim(max(min(extent_lons) - lon_range * 0.15, -180), min(max(extent_lons) + lon_range * 0.15, 180))
            ax.set_ylim(max(min(extent_lats) - lat_range * 0.15, -90), min(max(extent_lats) + lat_range * 0.15, 90))

    if service_mode:
        ax.set_title('Service and User Stations Distribution', fontsize=13, fontweight='bold')
    else:
        ax.set_title('Station distribution', fontsize=13, fontweight='bold')

    if outfname:
        fig.savefig(outfname, dpi=150, bbox_inches='tight')
        print(f'Plot written to {outfname}')
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Plot stations as blue circles by default, or service/user stations when --svr_usr is used')
    parser.add_argument('paths', nargs='*', help='file paths or glob patterns')
    parser.add_argument('-o', '--out', default='stations_by_observer.png', help='output image file')
    parser.add_argument('--sitename', action='store_true', help='annotate station names using the first four characters of each file name')
    parser.add_argument('--svr_usr', nargs=2, metavar=('SERVICE_PATH', 'USER_PATH'), help='plot service stations as red triangles and user stations as blue circles using two paths, each path can be a file or a directory')
    parser.add_argument('--grid', help='grid CSV file in CLAS format, e.g. grid/ga_grid_single_network.csv')
    parser.add_argument('--grid-points', action='store_true', help='draw the individual grid points from the grid CSV')
    args = parser.parse_args()

    files = []
    for p in args.paths:
        files.extend(resolve_input_paths(p))

    service_paths = user_paths = None
    if args.svr_usr:
        service_files = resolve_input_paths(args.svr_usr[0])
        user_files = resolve_input_paths(args.svr_usr[1])
        if not files:
            files = service_files + user_files
        service_paths = set(service_files)
        user_paths = set(user_files)

    grid_cells = None
    if args.grid:
        grid_cells = load_grid_csv(args.grid)
    show_grid_points = bool(grid_cells) or args.grid_points

    files = [f for f in files if os.path.isfile(f)]
    if not files:
        print('no files found for given patterns', file=sys.stderr)
        return

    data = gather(files)
    if not data:
        print('no station coordinates found', file=sys.stderr)
        return

    plot_by_observer(
        data,
        outfname=args.out,
        show_sitename=args.sitename,
        service_paths=service_paths,
        user_paths=user_paths,
        grid_cells=grid_cells,
        show_grid_points=show_grid_points,
    )


if __name__ == '__main__':
    main()
