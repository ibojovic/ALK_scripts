#!/bin/bash


# Output directory
output_dir=~/projects/alkFus_bench/temp/data/LUAD_out

# Iterate through folders in the input directory
for folder in ~/projects/alkFus_bench/downloads/TCGA/LUAD/*/; do
  # Extract sample name from folder path
  sample=$(basename "$folder")

  # Create output directories for each fusion caller command
  arriba_dir="$output_dir/arriba/$sample"
  starfusion_dir="$output_dir/starfusion/$sample"
  infusion_dir="$output_dir/infusion/$sample"
  fusioncatcher_dir="$output_dir/fusioncatcher/$sample"
  fcirc_dir="$output_dir/fcirc/$sample"
  jaffa_dir="$output_dir/jaffa/$sample"

  # Create the directories if they don't exist
  mkdir -p "$arriba_dir"
  mkdir -p "$starfusion_dir"
  mkdir -p "$infusion_dir"
  mkdir -p "$fusioncatcher_dir"
  mkdir -p "$fcirc_dir"
  mkdir -p "$jaffa_dir"

  # Find input files with full path
  input1=$(find "$folder" -name "*_1_fastq.gz" -type f)
  input2=$(find "$folder" -name "*_2_fastq.gz" -type f)

  # Perform fusion calling using arriba and define current aviable cores
  echo "Running arriba for sample: $sample"
  source activate arriba
  cd "$arriba_dir"
  cpu_count=$(echo "scale=0; ($(top -bn1 | awk '/^%Cpu/{print 100 - $8}') * $(nproc)) / 100" | bc)
  run_arriba.sh /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/STAR_index_GRCh38viral_ENSEMBL104/ /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/ENSEMBL104.gtf /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/GRCh38viral.fa /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3 "$cpu_count" "$input1" "$input2"
  find . ! -name "fusions.tsv" ! -name "fusions.discarded.tsv" -type f -delete
    
  
  # Perform fusion calling using star-fusion
  echo "Running star-fusion for sample: $sample"
  source activate star-fusion
  cd "$starfusion_dir"
  cpu_count=$(echo "scale=0; ($(top -bn1 | awk '/^%Cpu/{print 100 - $8}') * $(nproc)) / 100" | bc)
  STAR-Fusion --genome_lib_dir /home/labgroups/ccgg/anaconda3/envs/star-fusion/bin/CTAT_resource_lib/ctat_genome_lib_build_dir --left_fq "$input1" --right_fq "$input2" --CPU "$cpu_count"
  find . ! -name "star-fusion.fusion_predictions.tsv"  -type f -delete

  # Perform fusion calling using infusion
  echo "Running infusion for sample: $sample"
  source activate py2
  cd "$infusion_dir"
  python ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/infusion -1 "$input1" -2 "$input2" ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/index_folder_infusion/infusion.cfg
  find . ! -name "fusions.detailed.full.txt"  -type f -delete

  # Perform fusion calling using fusioncatcher
  echo "Running fusioncatcher for sample: $sample"
  source activate fusioncatcher
  cd "$fusioncatcher_dir"
  /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/bin/fusioncatcher -d /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/data/current -i "$sample" -o .
  find . ! -name "final-list_candidate-fusion-genes.txt" ! -name "summary_candidate_fusions.txt"  -type f -delete

  # Perform fusion calling using fcirc
  echo "Running fcirc for sample: $sample"
  source activate py3
  cd "$fcirc_dir"
  cpu_count=$(echo "scale=0; ($(top -bn1 | awk '/^%Cpu/{print 100 - $8}') * $(nproc)) / 100" | bc)
  ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fcirc.py -t "$cpu_count" -o . -x ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index/fusiongenes_ref_U -f ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index -1 "$input1" -2 "$input2"
  find . ! -name "fcircRNA_results.tsv" ! -name "fusion_results.tsv"  -type f -delete

  # Perform fusion calling using jaffa assembly
  echo "Running jaffa assembly for sample: $sample"
  source activate jaffa
  cd "$jaffa_dir"
  cpu_count=$(echo "scale=0; ($(top -bn1 | awk '/^%Cpu/{print 100 - $8}') * $(nproc)) / 100" | bc)
  ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/tools/bin/bpipe run -n "$cpu_count" ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/JAFFA_assembly.groovy  \"\$input1\" \"\$input2\"
  find . ! -name "jaffa_results.csv"  -type f -delete


done
