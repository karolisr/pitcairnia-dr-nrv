#!/usr/bin/env bash

run() {
    THETA="${1}"
    OUT_PATH="${2}"
    OUT_PREFIX="${OUT_PATH}/16-freebayes-theta-${THETA}"

    VCF_IN="${OUT_PREFIX}/nrv-pitcarn.vcf"
    VCF_OUT="${VCF_IN}.filt.01.vcf"

    bcftools view "${VCF_IN}" \
        | bcftools view -i 'TYPE="ref" | TYPE="snp"' \
        | bcftools view -i 'NUMALT<=1' \
        | bcftools view -e '(TYPE="snp") & (QUAL<20)' \
        | bcftools view -i '(MIN(FMT/DP)>=10) & (NS>=1)' \
        | bcftools view -e '(COUNT(GT="het")>0) & ((AB<0.25) | (AB>0.75))' \
        | bcftools view -e '(COUNT(GT="het")=0) & (AB>0.05)' \
        | bcftools view -e 'AC>=1 & MAF<0.05' \
        | bcftools plugin fill-tags \
    > "${VCF_OUT}"

    # ------------------------------------------------------------------------
    # run processes and store pids in array
    i=0
    pids=()
    for n in {10..19}; do
        i=$((i + 1))
        VCF_OUT_N="${VCF_IN}.filt.${n}.vcf"
        echo "${VCF_OUT_N}"
        bcftools view -i "NS>=${n}" "${VCF_OUT}" > "${VCF_OUT_N}" &
        pids[${i}]=$!
    done

    # wait for all pids
    for pid in "${pids[@]}"; do
        wait ${pid}
    done
    # ------------------------------------------------------------------------

    RSYNC_CMD="rsync -a ${OUT_PREFIX}/ /home/karolis/SyncThing/temp/${OUT_PREFIX##*/}/ &> /dev/null"
    eval "${RSYNC_CMD}"
}

export -f run

NCPU=2
THETAS=(0.00001 0.001)
OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"

echo "#B: $(date)"
# NOTE -----------------------------------------------------------------------
#   parallel: Warning: No more file handles.
#   parallel: Warning: Raising ulimit -n or /etc/security/limits.conf may help.
# ----------------------------------------------------------------------------
ulimit -n 48000 && \
(parallel -k -j "${NCPU}" run \
 ::: "${THETAS[@]}" \
 ::: "${OUT_PATH}")
echo "#E: $(date)"
