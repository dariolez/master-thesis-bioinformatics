#!/bin/bash
#
# AUTHOR: Darío González
# DESCRIPTION: Shell script to run STAR in the lab server folia.

# Assign arguments to variables
input_file=$1

# Change to the working directory in folia
mkdir -p /dario-work/TFM/2019.12.31_Meyer/STAR_align/
cd /dario-work/TFM/2019.12.31_Meyer/STAR_align

# Read names of files and process them
while read -r folder
do
    # State the file being processed
    run=$(basename ${folder})
    echo -e "\n[ $(date +'%Y/%m/%d %T.%3N') ] Processing files in ${run}\n"

    # Transfer files from nuptse to folia
    rsync -avR --progress -hh dario@nuptse:/mnt/bazar/dario_TFM/2019.12.31_Meyer/data+analysis/esophagus/singlecell/data/./${run}/*.fastq.gz .

    # Create a folder for the results
    cd ${run}
    #rm -rf STAR_out
    mkdir -p STAR_out

    # Process files
    ~/bin/STAR_2.7.11a/Linux_x86_64_static/STAR \
        --runMode alignReads \
        --genomeDir /home/dario/TFM/reference_genomes/STAR_index/ \
        --outFileNamePrefix STAR_out/ \
        --readFilesCommand zcat \
        --readFilesIn *_R2_*.fastq.gz *_R1_*.fastq.gz \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist /home/dario/TFM/10x_barcodes/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloStrand Forward \
        --soloFeatures Gene GeneFull Velocyto \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped unmapped.fa \
        --runThreadN 10

    # Transfer result files back to nuptse
    rsync -av --progress -hh STAR_out dario@nuptse:/mnt/bazar/dario_TFM/2019.12.31_Meyer/data+analysis/esophagus/singlecell/data/${run}/

    # Clean up
    cd ..
    rm -rf ${run}
    echo -e "\n[ $(date +'%Y/%m/%d %T.%3N') ] Finished processing ${run}"
done < ${input_file}

# Remove working directory
rmdir /dario-work/TFM/2019.12.31_Meyer/STAR_align/
