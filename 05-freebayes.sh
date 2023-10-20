#!/usr/bin/env bash

fb() {
     REGION="${1}"
     ASSMBL="${2}"
     BAM_LIST="${3}"
     FREEBAYES_BIN="${4}"
     LOG_F="${5}"
     # RSYNC_CMD="${6}"
     THETA="${7}"
     # VCF="${LOG_F%/*}/${REGION%:*}"

     SECONDS=0
     # ----------------------------------------------------------------------------
     "${FREEBAYES_BIN}" \
          --fasta-reference "${ASSMBL}" \
          --bam-list "${BAM_LIST}" \
          --region "${REGION}" \
          --ploidy 2 \
          --theta "${THETA}" \
          --report-monomorphic
     # ----------------------------------------------------------------------------
     # printf -v SECONDS_FMT %6.6s "${SECONDS}"
     SECONDS_FMT="${SECONDS}"

     range="${REGION##*:}"
     start="${range%-*}"
     end="${range##*-}"
     len="$((end - start))"

     # echo -e "${SECONDS_FMT}\\t${len}\\t${REGION}\\t$(ulimit -n)" >> "${LOG_F}"
     echo -e "${SECONDS_FMT}\\t${len}\\t${REGION}" >> "${LOG_F}"

     # eval "${RSYNC_CMD}"
}

export -f fb

NCPU=14
# THETA=0.00001
THETA=0.001

FREEBAYES_BIN="/home/karolis/.local/bin/freebayes-1.3.6-linux-amd64-static"

OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"
ASSMBL="${OUT_PATH}/nrv_reference_transcriptome.fa"
ASSMBL_INDEX="${ASSMBL}.fai"
BAM_PREFIX="${OUT_PATH}/15-sorted_dupes_marked_bam"
BAM_FILES=("$(ls -Sr ${BAM_PREFIX}/*.rgsort.sorted.dupes_marked.bam)")
OUT_PREFIX="${OUT_PATH}/16-freebayes-theta-${THETA}"
BAM_LIST="${OUT_PREFIX}/bam-list.txt"
VCF="${OUT_PREFIX}/nrv-pitcarn.vcf"
LOG_F="${VCF}.log"

mkdir -p "${OUT_PREFIX}"

# TMP_DIR="${OUT_PATH}/temp"
# mkdir -p "${TMP_DIR}"

echo "${BAM_FILES[@]}" > "${BAM_LIST}"

FB_SCRIPTS="freebayes-scripts"
mapfile -t regions < <("${FB_SCRIPTS}/fasta_generate_regions.py" "${ASSMBL_INDEX}" "100000000000")

# ----------------------------------------------------------------------------
# echo
# regions_filtered=()
# for region in "${regions[@]}"; do
#   # echo -en "${region}\\t\\t"
#   range="${region##*:}"
#   start="${range%-*}"
#   end="${range##*-}"
#   # echo "${start}"
#   # echo "${end}"
#   len="$((end - start))"
#   # echo "${len}"

#   if [[ len -lt 600 ]]; then
#     regions_filtered+=("${region}")
#     printf -v region %-45.45s "$region"
#     echo "${region}${len}"
#   fi

#   if [[ ${#regions_filtered[@]} -eq 200 ]]; then
#      break
#   fi

# done
# echo
# regions=("${regions_filtered[@]}")
# ----------------------------------------------------------------------------

RSYNC_CMD="rsync -a ${OUT_PREFIX}/ /home/karolis/SyncThing/temp/${OUT_PREFIX##*/}/ &> /dev/null"

echo "#B: $(date)" > "${LOG_F}"
# NOTE -----------------------------------------------------------------------
#   parallel: Warning: No more file handles.
#   parallel: Warning: Raising ulimit -n or /etc/security/limits.conf may help.
# ----------------------------------------------------------------------------
ulimit -n 48000 && \
(parallel -k -j "${NCPU}" fb \
 ::: "${regions[@]}" \
 ::: "${ASSMBL}" \
 ::: "${BAM_LIST}" \
 ::: "${FREEBAYES_BIN}" \
 ::: "${LOG_F}" \
 ::: "${RSYNC_CMD}" \
 ::: "${THETA}") | vcffirstheader | vcfstreamsort -w 10000 | vcfuniq > "${VCF}"
echo "#E: $(date)" >> "${LOG_F}"

# cat -p "${LOG_F}" | grep -v "#" | awk '{print $1/60/60,$2,$3}' | sort --numeric-sort | wc -l
# vcf2tsv "${VCF}" > "${VCF}.tsv"
eval "${RSYNC_CMD}"
