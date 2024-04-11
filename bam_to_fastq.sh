#!/bin/bash

# Find folders ending with "_all" in the current directory
for folder in *_all/; do
    # Check if the folder exists and is a directory
    if [[ -d "$folder" ]]; then
        # Get the folder name without the trailing slash
        folder_name="${folder%/}"
        # Extract the prefix from the folder name
        prefix="${folder_name%%_*}"
        
        # Find genomic and chimeric BAM files in the folder
        genomic_file=$(find "$folder" -type f -name "*.genomic*.bam")
        chimeric_file=$(find "$folder" -type f -name "*.chimeric*.bam")
        
        # Check if both genomic and chimeric BAM files exist
        if [[ -n "$genomic_file" && -n "$chimeric_file" ]]; then
            # Extract the file names without extensions
            genomic_filename=$(basename "$genomic_file" .bam)
            chimeric_filename=$(basename "$chimeric_file" .bam)
                        
            # Perform samtools merge and sort
            samtools merge -@ 10 - "$genomic_file" "$chimeric_file" | samtools sort -@ 10 -n -o "$folder/${chimeric_filename}_merged_sorted.bam" -
            
            # Delete the genomic and chimeric BAM files
            rm "$genomic_file" "$chimeric_file"
            
            # Generate paired-end FASTQ files from the merged and sorted BAM file
            samtools fastq -@ 10 "$folder/${chimeric_filename}_merged_sorted.bam" -1 "$folder/${prefix}_1.fastq.gz" -2 "$folder/${prefix}_2.fastq.gz" -0 /dev/null -s /dev/null -n
                       
            # Delete the merged and sorted BAM file
            rm "$folder/${chimeric_filename}_merged_sorted.bam"
        fi
    fi
done

