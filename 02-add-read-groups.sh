#!/usr/bin/env bash

NTHREADS=6

PICARD_BIN="PicardCommandLine"

FQ_PATH="/home/karolis/scratch/nrv-pitcairnia-transcriptomes"
OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"
TMP_DIR="${OUT_PATH}/temp"
ASSMBL="${OUT_PATH}/nrv_reference_transcriptome.fa"
IN_PREFIX="${OUT_PATH}/12-star_pass_2"
OUT_PREFIX="${OUT_PATH}/13-star_pass_2_rgsort"
SAM_FILE_SUFFIX="_Genomes_Aligned.out.sam"

mkdir -p "${TMP_DIR}"
mkdir -p "${OUT_PREFIX}"

for INPUT_SAM_FILE in $(ls "${IN_PREFIX}" | grep "${SAM_FILE_SUFFIX}"); do

    ((i=i%NTHREADS)); ((i++==0)) && wait

    SAMPLE_NAME="${INPUT_SAM_FILE/_Genomes_Aligned.out.sam/}"

    if [[ ! "${SAMPLE_NAME}" =~ PI[AS][LT] ]]; then

        IF="${IN_PREFIX}/${INPUT_SAM_FILE}"
        OF="${OUT_PREFIX}/${SAMPLE_NAME}.rgsort.sam"
        LF="${OUT_PREFIX}/${SAMPLE_NAME}.rgsort.log"

        echo "${SAMPLE_NAME}" :: "${IF##*/}" --\> "${OF##*/}"

        "${PICARD_BIN}" AddOrReplaceReadGroups \
            TMP_DIR="${TMP_DIR}" \
            I="${IF}" \
            O="${OF}" \
            SO=coordinate \
            RGID="${SAMPLE_NAME}" \
            RGLB="LIB1" \
            RGPL="ILLUMINA" \
            RGPU="HCVWHDSX3.3" \
            RGSM="${SAMPLE_NAME}" &> "${LF}" &
    fi
done

# SO=SortOrder Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.
# ID=String  Read-Group ID Default value: 1. This option can be set to 'null' to clear the default value.
# LB=String  Read-Group library. [Required]
# PL=String  Read-Group platform (e.g. illumina, solid). [Required]
# PU=String  Read-Group platform unit (eg. run barcode). [Required]
# SM=String  Read-Group sample name. [Required]
