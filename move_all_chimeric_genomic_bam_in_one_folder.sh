#This script is made for files that are arriving from TCGA
#all files downloaded from there are separated in different folders 
# since i need to combine chimeric and genomic files in one i need them in one folder
#this script is in directory: ~/projects/alkFus_bench/downloads/TCGA/LUAD


#!/bin/bash

# Read the file and process each line
while IFS=$'\t' read -r id filename _; do
    # Check if the filename has "chimeric" in its name
    if [[ $filename == *".chimeric."* ]]; then
        # Extract the prefix of the filename
        prefix="${filename%%.*}"

        # Determine the destination folder based on the ID
        destination_folder="${id}_all"

        # Create the destination folder if it doesn't exist
        if [[ ! -d $destination_folder ]]; then
            mkdir -p "$destination_folder"
        fi

        # Find all files with the same prefix in subdirectories and move them to the destination folder
        find . -type f -name "${prefix}.*" -exec mv {} "$destination_folder/" \;
    fi

done < "download_gdc_manifest.2023-06-04.txt"

