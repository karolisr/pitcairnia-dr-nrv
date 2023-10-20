#!/usr/bin/env bash

NTHREADS=12

STAR_BIN="STAR"
# STAR_BIN="${HOME}/.local/bin/STAR"
# STAR_BIN="/opt/homebrew/bin/STAR"

FQ_PATH="/home/karolis/scratch/nrv-pitcairnia-transcriptomes"
OUT_PATH="/home/karolis/scratch/nrv-pop-gen-analyses"
ASSMBL="${OUT_PATH}/nrv_reference_transcriptome.fa"

# Note: STAR REFERENCE INDEX ("GENOME")
# "${STAR_BIN}" \
# --runThreadN ${NTHREADS} \
# --runMode genomeGenerate \
# --genomeDir "${ASSMBL}.star" \
# --genomeFastaFiles "${ASSMBL}" \
# --genomeSAindexNbases 10

# Note: STAR PASS 1
# for transname in $(ls "${FQ_PATH}"); do
#     transpath="${FQ_PATH}/${transname}/Nuclear"
#     echo "${transpath}"
#         "${STAR_BIN}" \
#             --runThreadN ${NTHREADS} \
#             --genomeDir "${ASSMBL}.star" \
#             --readFilesIn \
#                 "${transpath}/nuclear_1.fq.gz" \
#                 "${transpath}/nuclear_2.fq.gz" \
#             --readFilesCommand gunzip \
#             -c \
#             --outFileNamePrefix "${OUT_PATH}/11-star_pass_1/${transname}_"
# done

# Note: GENERATE A UNIQUE JUNCTIONS FILE
# cat "${OUT_PATH}/11-star_pass_1/"*_SJ.out.tab | sort | uniq > "${OUT_PATH}/11-star_pass_1/unique_junctions.tab"

# Note: STAR PASS 2
for transname in $(ls "${FQ_PATH}"); do
    transpath="${FQ_PATH}/${transname}/Nuclear"
    echo "${transpath}"
        "${STAR_BIN}" \
            --runThreadN ${NTHREADS} \
            --genomeDir "${ASSMBL}.star" \
            --readFilesIn \
                "${transpath}/nuclear_1.fq.gz" \
                "${transpath}/nuclear_2.fq.gz" \
            --readFilesCommand gunzip \
            -c \
            --sjdbFileChrStartEnd "${OUT_PATH}/11-star_pass_1/unique_junctions.tab" \
            --outFileNamePrefix "${OUT_PATH}/12-star_pass_2/${transname}_"
done
