#!/bin/bash
##
## usage: ./crx_downloaded_fromGA.sh /home/zhouforme/Documents/MySoftwares/proj_Cpp/UNIQ/GOOD/dataset_Linux 2023 1 2
##
datadir=$1  # data directory
year=$2     # 4-digit year
doy=$3      # day of year
ndays=$4    # number of consecutive days
yy=`echo ${year:2:2}`  # 2-digit year
USER=anonymous
PASSWORD=zhouforme@163.com

# ---------------------------------------------------------------------------
# Site list configuration (EDIT HERE)
# ---------------------------------------------------------------------------
# Path to the file containing station names or full CRX filenames.
# Each line: either a 4-char site code (e.g. "ADMN") or a full filename
# (e.g. "ADMN00AUS_S_20230700000_01D_30S_MO.crx.gz").
# If left empty or the file does not exist, the script falls back to:
#   ${datadir}/site_ga.list
# ---------------------------------------------------------------------------
site_list="/home/lzx/code/GREAT-PVT/data/ga_1d3_trimble_alloy_100uniform.txt"
# site_list="/path/to/your/site.list"
# site_list="${datadir}/site_ga.list"

crx2rnx_bin=./script/crx2rnx

download_one_site() {
    local year="$1"
    local cdoy="$2"
    local site_new="$3"

    if command -v lftp >/dev/null 2>&1; then
        lftp -u ${USER},${PASSWORD} sftp://52.65.50.124<<EOF  # lftp needs to be installed in advance
            set sftp:auto-confirm yes
            cd rinex/daily/${year}/${cdoy}
            get ${site_new}
            bye
EOF
    else
        # Fallback for environments where lftp is not installed.
        curl --silent --show-error --fail --user "${USER}:${PASSWORD}" \
            "sftp://52.65.50.124/rinex/daily/${year}/${cdoy}/${site_new}" \
            --output "${site_new}"
    fi
}

to_rinex_name() {
    local fpath="$1"
    local fname
    fname=$(basename "$fpath")

    if [[ "$fname" == *.crx.gz ]]; then
        echo "${fpath%.crx.gz}.rnx"
    elif [[ "$fname" == *.crx ]]; then
        echo "${fpath%.crx}.rnx"
    elif [[ "$fname" =~ \.[0-9][0-9]d$ ]]; then
        echo "${fpath%?}o"
    else
        echo "${fpath}.rnx"
    fi
}

postprocess_crx_gz() {
    local gz_path="$1"
    local crx_path="${gz_path%.gz}"
    local rnx_path
    rnx_path=$(to_rinex_name "$gz_path")

    if [[ ! -f "$gz_path" ]]; then
        echo "Missing file: $gz_path" >&2
        return 1
    fi

    if ! command -v gzip >/dev/null 2>&1; then
        echo "gzip is required but not found." >&2
        return 1
    fi

    if [[ ! -x "$crx2rnx_bin" ]]; then
        chmod +x "$crx2rnx_bin" 2>/dev/null || true
    fi
    if [[ ! -x "$crx2rnx_bin" ]]; then
        echo "crx2rnx is not executable: $crx2rnx_bin" >&2
        return 1
    fi

    if ! gzip -d "$gz_path"; then
        echo "Failed to gunzip: $gz_path" >&2
        return 1
    fi

    if "$crx2rnx_bin" "$crx_path" - -f > "$rnx_path"; then
        rm -f "$crx_path"
        return 0
    fi

    rm -f "$rnx_path"
    echo "Failed to convert CRX to RINEX: $crx_path" >&2
    return 1
}

while [[ ${ndays} -gt 0 ]]; do  # day-by-day
    # Resolve site_list: use explicit setting if valid, otherwise fall back to datadir
    if [[ -z "${site_list}" || ! -f "${site_list}" ]]; then
        site_list="${datadir}/site_ga.list"
    fi

    cdoy=`echo ${doy} | awk '{printf("%3.3d\n", $1)}'`
    obsdir=${datadir}/obs/${year}/${cdoy}  # the subdirectory for observations
    if [[ ! -d "${obsdir}" ]]; then
        mkdir -p ${obsdir}
    fi
    # 如果用户提供的 site_list 只包含短站码（如 ADDE），则尝试从模板生成完整文件名列表
    firstline=$(grep -m1 -E "\S" "${site_list}" 2>/dev/null | tr -d '\r' || true)
    if [[ -n "$firstline" && "$firstline" =~ ^[A-Za-z]{3,6}$ ]]; then
        # 优先使用全名模板（GAMPII 提供的完整文件名列表），其次才是 datadir 下的列表
        template_candidates=("/home/zxli/GAMPII-GOOD-master/dataset_Linux/site_ga.list" "${datadir}/site_ga.list")
        template=""
        for t in "${template_candidates[@]}"; do
            if [[ -f "$t" ]]; then
                template="$t"
                break
            fi
        done

        tmp_list="${datadir}/site_ga.full.list"
        : > "$tmp_list"
        while IFS= read -r code || [[ -n "$code" ]]; do
            code=$(echo "$code" | tr -d '\r' | tr '[:lower:]' '[:upper:]' | awk '{print $1}')
            if [[ -z "$code" ]]; then
                continue
            fi
            if [[ -n "$template" ]]; then
                # 在模板中查找以站码开头的完整文件名
                match=$(grep -m1 -E "^${code}" "$template" || true)
                if [[ -n "$match" ]]; then
                    # 确认 match 看起来像完整文件名（包含 .crx 或 .crx.gz 或包含年份标记）
                    if echo "$match" | grep -Eq '\.crx(\.gz)?$|2023[0-9]{3}|202[0-9]{3}'; then
                        echo "$match" >> "$tmp_list"
                        continue
                    else
                        match=""
                    fi
                fi
            fi
            # 没有模板或未找到匹配项，使用默认命名模式
            echo "${code}00AUS_S_20230010000_01D_30S_MO.crx.gz" >> "$tmp_list"
        done < "${site_list}"

        # 使用生成的完整文件名列表
        site_list="$tmp_list"
        echo "使用生成的完整文件名列表: $site_list" >&2
    fi
    while IFS= read -r site_old || [[ -n "$site_old" ]]; do  # site-by-site
        site_old=$(echo "$site_old" | tr -d '\r')
        [[ -z "$site_old" ]] && continue

        site_new=`echo ${site_old/2023001/${year}${cdoy}}`
        rnx_new=`echo ${site_new%.crx.gz}.rnx`
        echo ${site_new}
        if [[ -f ${obsdir}/${rnx_new} ]]; then
            continue  # final rinex exists
        fi
        if [[ -f ${obsdir}/${site_new} ]]; then
            continue  # file exists
        fi
        if download_one_site "${year}" "${cdoy}" "${site_new}"; then
            mv ${site_new} ${obsdir}
            if ! postprocess_crx_gz "${obsdir}/${site_new}"; then
                echo "Post-process failed for ${obsdir}/${site_new}" >&2
            fi
        else
            echo "Failed to download ${site_new}" >&2
        fi
    done < "${site_list}"
    let "doy+=1"
    let "ndays-=1"
done