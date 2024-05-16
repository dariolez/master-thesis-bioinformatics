#!/bin/bash
#
# AUTHOR: Darío González
# DESCRIPTION: Calculate loss of sequencing reads in each step of the scRNA-seq analysis pipeline

# Create function to print statistics for each run
function reads_stats() {
    
    FOLDER=$1

    # get read counts for each step of the pipeline
    original=$(unzip -p fastqc/*_R2_*_fastqc.zip *_R2_001_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print $3}')
    corrected=$(grep -E "yes*" ${FOLDER}/Solo.out/Barcodes.stats | sed -E 's/^\s+//g; s/\s+/\t/g' | awk 'BEGIN{total=0} {total += $2} END{print total}')
    unique_aligned=$(grep "Uniquely mapped reads number" ${FOLDER}/Log.final.out | sed -E "s/^\s+//g; s/\s{2,}/\t/g" | awk 'BEGIN{FS="\t"} {print $2}')
    spliced=$(tail -n +4 ${FOLDER}/Solo.out/Velocyto/raw/spliced.mtx | cut -d ' ' -f 3 | awk 'BEGIN{total=0} {total += $1} END{print total}')
    unspliced=$(tail -n +4 ${FOLDER}/Solo.out/Velocyto/raw/unspliced.mtx | cut -d ' ' -f 3 | awk 'BEGIN{total=0} {total += $1} END{print total}')
    ambiguous=$(tail -n +4 ${FOLDER}/Solo.out/Velocyto/raw/ambiguous.mtx | cut -d ' ' -f 3 | awk 'BEGIN{total=0} {total += $1} END{print total}')
    spliced_unspliced=$(($spliced + $unspliced))
    total=$(($spliced + $unspliced + $ambiguous))
    

    # Calculate percentages
    original_pc=100
    corrected_pc=$(echo "(${corrected}/${original}) * 100" | bc -l)
    unique_aligned_pc=$(echo "(${unique_aligned}/${original}) * 100" | bc -l)
    spliced_pc=$(echo "(${spliced}/${original}) * 100" | bc -l)
    unspliced_pc=$(echo "(${unspliced}/${original}) * 100" | bc -l)
    ambiguous_pc=$(echo "(${ambiguous}/${original}) * 100" | bc -l)
    spliced_unspliced_pc=$(echo "(${spliced_unspliced}/${original}) * 100" | bc -l)
    total_pc=$(echo "(${total}/${original}) * 100" | bc -l)

    # Output read counts
    LC_NUMERIC=en_US.UTF-8  # to have the decimal separator set as a point

    printf "step\tnumber_reads\tpercentage\n"
    printf "original\t%s\t%.2f\n" $original $original_pc
    printf "corrected\t%s\t%.2f\n" $corrected $corrected_pc
    printf "unique_aligned\t%s\t%.2f\n" $unique_aligned $unique_aligned_pc
    printf "spliced\t%s\t%.2f\n" $spliced $spliced_pc
    printf "unspliced\t%s\t%.2f\n" $unspliced $unspliced_pc
    printf "ambiguous\t%s\t%.2f\n" $ambiguous $ambiguous_pc
    printf "spliced+unspliced\t%s\t%.2f\n" $spliced_unspliced $spliced_unspliced_pc
    printf "spliced+unspliced+ambiguous\t%s\t%.2f\n" $total $total_pc
}

# Get folder names
run_folders=$1

# Set variables
FOLDER=STAR_out

# Run stats function for all runs
while read -r run_folder; do
    cd $run_folder
    run=$(basename ${run_folder})
    
    mkdir -p ${run_folder}/${FOLDER}/read_stats
    reads_stats $FOLDER > ${run_folder}/${FOLDER}/read_stats/read_stats_STAR.tsv

    cd ..
done < ${run_folders}
