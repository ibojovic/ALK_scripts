#!/bin/bash

export PATH="/home/labgroups/ccgg/anaconda3/bin:$PATH"

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
  input1=$(find "$folder" -name "*_1.fastq.gz" -type f)
  input2=$(find "$folder" -name "*_2.fastq.gz" -type f)

  # Define the commands to be executed in parallel
  commands=(
    "echo 'Running arriba for sample: $sample'; conda activate arriba && cd '$arriba_dir' && /home/labgroups/ccgg/anaconda3/envs/arriba/bin/run_arriba.sh /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/STAR_index_GRCh38viral_ENSEMBL104/ /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/ENSEMBL104.gtf /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/GRCh38viral.fa /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3 \"$cpu_count\" \"$input1\" \"$input2\" && find . ! -name 'fusions.tsv' ! -name 'fusions.discarded.tsv' -type f -delete && conda deactivate"
    "echo 'Running star-fusion for sample: $sample'; conda activate star-fusion && cd '$starfusion_dir' && STAR-Fusion --genome_lib_dir /home/labgroups/ccgg/anaconda3/envs/star-fusion/bin/CTAT_resource_lib/ctat_genome_lib_build_dir --left_fq \"$input1\" --right_fq \"$input2\" --CPU \"$cpu_count\" && find . ! -name 'star-fusion.fusion_predictions.tsv'  -type f -delete && conda deactivate"
    "echo 'Running infusion for sample: $sample'; conda activate py2 && cd '$infusion_dir' && sed -i \"s/^num_threads = .*/num_threads = $cpu_count/\" ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/index_folder_infusion/infusion.cfg && python ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/infusion -1 \"$input1\" -2 \"$input2\" ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/index_folder_infusion/infusion.cfg && find . ! -name 'fusions.detailed.full.txt'  -type f -delete && conda deactivate"
    "echo 'Running fusioncatcher for sample: $sample'; conda activate fusioncatcher && cd '$fusioncatcher_dir' && sed -i \"s/^threads = .*/threads = $cpu_count/\" /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/etc/configuration.cfg && /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/bin/fusioncatcher -d /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/data/current -i \"$sample\" -o . && find . ! -name 'final-list_candidate-fusion-genes.txt' ! -name 'summary_candidate_fusions.txt'  -type f -delete && conda deactivate"
    "echo 'Running fcirc for sample: $sample'; conda activate py3 && cd '$fcirc_dir' && ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fcirc.py -t \"$cpu_count\" -o . -x ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index/fusiongenes_ref_U -f ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index -1 \"$input1\" -2 \"$input2\" && find . ! -name 'fcircRNA_results.tsv' ! -name 'fusion_results.tsv'  -type f -delete && conda deactivate"
    "echo 'Running jaffa assembly for sample: $sample'; conda activate jaffa && cd '$jaffa_dir' && ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/tools/bin/bpipe run -n \"$cpu_count\" ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/JAFFA_assembly.groovy \"$input1\" \"$input2\" && find . ! -name 'jaffa_results.csv'  -type f -delete && conda deactivate"
  )

  # Execute commands in parallel
  parallel --halt-on-error now,fail=1 ::: "${commands[@]}"
done
