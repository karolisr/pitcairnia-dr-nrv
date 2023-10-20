#!/usr/bin/env bash

OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"
IN_PREFIX="${OUT_PATH}/13-star_pass_2_rgsort"
OUT_PREFIX="${OUT_PATH}/14-sam_to_sorted_bam"
TMP_DIR="${OUT_PATH}/temp"

mkdir -p "${OUT_PREFIX}"
mkdir -p "${TMP_DIR}"

sam2bam () {

	echo -ne "\\t$(date)"

	NTHREADS=2
	SAMTOOLS_BIN="/usr/bin/samtools"

	SAM="${1}"
	OUT_PREFIX="${2}"
	TMP_DIR="${3}"

	BAM="${SAM/.sam/.bam}"
	BAM="${BAM##*/}"
	BAM="${TMP_DIR}/${BAM}"

	BAM_SORTED="${2}/${BAM/.bam/.sorted.bam}"
	BAM_SORTED="${BAM_SORTED##*/}"
	LOG_F="${OUT_PREFIX}/${BAM_SORTED}.log"
	BAM_SORTED="${OUT_PREFIX}/${BAM_SORTED}"

	echo -ne " | ${SAM##*/} --> ${BAM_SORTED##*/} | "

	"${SAMTOOLS_BIN}" view  -bh --verbosity 0 -@ ${NTHREADS} -o "${BAM}"        "${SAM}"        &>  "${LOG_F}"
	"${SAMTOOLS_BIN}" sort      --verbosity 0 -@ ${NTHREADS} -o "${BAM_SORTED}" "${BAM}"        &>> "${LOG_F}"
	"${SAMTOOLS_BIN}" index -b                -@ ${NTHREADS}                    "${BAM_SORTED}" &>> "${LOG_F}"

	rm "${BAM}"

	echo -ne "$(date)"
	echo
}

export -f sam2bam

echo -ne "\\nB:\\t$(date)"
echo
IN_FILES=("$(ls -S ${IN_PREFIX}/*.rgsort.sam)")
# echo "${IN_FILES[@]}"
parallel -j 15 sam2bam ::: "${IN_FILES[@]}" ::: "${OUT_PREFIX}" ::: "${TMP_DIR}"
echo -ne "E:\\t$(date)\\n"
