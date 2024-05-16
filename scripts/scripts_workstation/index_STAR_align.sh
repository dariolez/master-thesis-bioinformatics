#!/bin/bash
#
# DESCRIPTION: Script to index BAM files from STAR
# AUTHOR: Darío González

wd=/home/dario/TFM/2019.12.31_Meyer/data
sample_folders=$1

while read -r folder; do
    # State the file being processed
    run=$(basename $folder)
    echo -e "\n[ $(date +'%Y/%m/%d %T.%3N') ] Processing files in $run\n"

    # Prepare working directory and output directory
    wd_run=${wd}/${run}

    # Transfer files from nuptse to folia
    rsync -avR --progress -hh --include="*.bam" --include="*.bai" --exclude="/*/*/*" \
        dario@nuptse:$(dirname $folder)/./${run}/STAR_out/ ${wd}

    # Index BAM file
    bam_dir=${wd_run}/STAR_out

    if [[ ! -e ${bam_dir}/Aligned.sortedByCoord.out.bam.bai ]]; then
        samtools index -b -@ 4 ${bam_dir}/Aligned.sortedByCoord.out.bam
    fi

    # Transfer result files back to nuptse
    rsync -av --progress -hh ${bam_dir}/*.bai dario@nuptse:${folder}/STAR_out/

    # Clean up
    #rm -rf $wd_run
    rm -f ${bam_dir}/*.bam
    rm -f ${bam_dir}/*.bai

    # State end of processing
    echo -e "\n[ $(date +'%Y/%m/%d %T.%3N') ] Finished processing $run"

done < $sample_folders
