#!/usr/bin/env bash

NCPU=5

OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"
IN_PREFIX="${OUT_PATH}/14-sam_to_sorted_bam"
OUT_PREFIX="${OUT_PATH}/15-sorted_dupes_marked_bam"
TMP_DIR="${OUT_PATH}/temp"

mkdir -p "${OUT_PREFIX}"
mkdir -p "${TMP_DIR}"

markdups() {
	echo -ne "\\t$(date)"

	PICARD_BIN="PicardCommandLine"

	BAM_SORTED="${1}"
	OUT_PREFIX="${2}"
	TMP_DIR="${3}"

	BAM_DUPES="${BAM_SORTED/.bam/.dupes_marked.bam}"
	BAM_DUPES="${BAM_DUPES##*/}"
	LOG_F="${OUT_PREFIX}/${BAM_DUPES}.log"
	TXT_METRICS="${OUT_PREFIX}/${BAM_DUPES}.txt"
	BAM_DUPES="${OUT_PREFIX}/${BAM_DUPES}"

	echo -ne " | ${BAM_SORTED##*/} --> ${BAM_DUPES##*/} | "

	"${PICARD_BIN}" MarkDuplicates \
		TMP_DIR="${TMP_DIR}" \
		I="${BAM_SORTED}" \
		O="${BAM_DUPES}" \
		M="${TXT_METRICS}" \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT &> "${LOG_F}"

	echo -ne "$(date)"
	echo
}

export -f markdups

echo -ne "\\nB:\\t$(date)"
echo
IN_FILES=("$(ls -S ${IN_PREFIX}/*.rgsort.sorted.bam)")
# echo "${IN_FILES[@]}"
parallel -j ${NCPU} markdups ::: "${IN_FILES[@]}" ::: "${OUT_PREFIX}" ::: "${TMP_DIR}"
echo -ne "E:\\t$(date)\\n"
