#!/bin/bash

export PATH="/home/labgroups/ccgg/anaconda3/bin:$PATH"

# Output directory
output_dir=~/projects/alkFus_bench/temp/data/LUAD_out

# Define a function that returns the number of cores to use for each command
cores_per_command() {
  # Get the number of online processors
  local nproc=$(getconf _NPROCESSORS_ONLN)
  # Get the load average for the last minute
  local load=$(cut -d " " -f 1 /proc/loadavg)
  # Calculate the number of free processors
  local free=$(echo "$nproc - $load" | bc)
  # If there are no free processors, return 1
  if [ "$free" -le 0 ]; then
    echo 1
  else
    # Otherwise, return the number of free processors divided by the number of commands
    echo $(echo "$free / ${#commands[@]}" | bc)
  fi
}


# Define the commands to be executed in parallel
commands=(
    "echo 'Running arriba for sample: {#}'; conda activate arriba && cd '$output_dir/arriba/{#}' && /home/labgroups/ccgg/anaconda3/envs/arriba/bin/run_arriba.sh /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/STAR_index_GRCh38viral_ENSEMBL104/ /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/ENSEMBL104.gtf /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/GRCh38viral.fa /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz /home/labgroups/ccgg/anaconda3/envs/arriba/var/lib/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3 "{= $_=cores_per_command =}\" {1} {2} && find . ! -name 'fusions.tsv' ! -name 'fusions.discarded.tsv' -type f -delete && conda deactivate"

    "echo 'Running star-fusion for sample: {#}'; conda activate star-fusion && cd '$output_dir/starfusion/{#}' && STAR-Fusion --genome_lib_dir /home/labgroups/ccgg/anaconda3/envs/star-fusion/bin/CTAT_resource_lib/ctat_genome_lib_build_dir --left_fq {1} --right_fq {2} --CPU "{= $_=cores_per_command =}\" && find . ! -name 'star-fusion.fusion_predictions.tsv'  -type f -delete && conda deactivate"
    
    "echo 'Running infusion for sample: {#}'; conda activate py2 && cd '$output_dir/infusion/{#}' && sed -i \"s/^num_threads = .*/num_threads = "{= $_=cores_per_command =}\" ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/index_folder_infusion/infusion.cfg && python ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/infusion -1 {1} -2 {2} ~/projects/alkFus_bench/downloads/fusion_callers/InFusion-0.8/index_folder_infusion/infusion.cfg && find . ! -name 'fusions.detailed.full.txt'  -type f -delete && conda deactivate"
    
    "echo 'Running fusioncatcher for sample: {#}'; conda activate fusioncatcher && cd '$output_dir/fusioncatcher/{#}' && sed -i \"s/^threads = .*/threads = "{= $_=cores_per_command =}\" /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/etc/configuration.cfg && /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/bin/fusioncatcher -d /home/labgroups/ccgg/anaconda3/envs/fusioncatcher/data/current -i \"{#}\" -o . && find . ! -name 'final-list_candidate-fusion-genes.txt' ! -name 'summary_candidate_fusions.txt'  -type f -delete && conda deactivate"
    
    "echo 'Running fcirc for sample: {#}'; conda activate py3 && cd '$output_dir/fcirc/{#}' && ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fcirc.py -t "{= $_=cores_per_command =}\" -o . -x ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index/fusiongenes_ref_U -f ~/projects/alkFus_bench/downloads/fusion_callers/fcirc/fusion_total_index -1 {1} -2 {2} && find . ! -name 'fcircRNA_results.tsv' ! -name 'fusion_results.tsv'  -type f -delete && conda deactivate"
    
    "echo 'Running jaffa assembly for sample: {#}'; conda activate jaffa && cd '$output_dir/jaffa/{#}' && ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/tools/bin/bpipe run -n "{= $_=cores_per_command =}\" ~/projects/alkFus_bench/downloads/fusion_callers/JAFFA-version-2.3/JAFFA_assembly.groovy {1} {2} && find . ! -name 'jaffa_results.csv'  -type f -delete && conda deactivate"
)    
	
# Execute commands in parallel
  parallel --progress --jobs 100% --limit '~/resource_limit.sh' --xapply --results "$output_dir/{#}" --shuf --delay 7 --load 75% --nice 10 ::: "${commands[@]}" ::: $(find ~/projects/alkFus_bench/downloads/TCGA/LUAD/*/ -name "*_1.fastq.gz" -type f) ::: $(find ~/projects/alkFus_bench/downloads/TCGA/LUAD/*/ -name "*_2.fastq.gz" -type f)
