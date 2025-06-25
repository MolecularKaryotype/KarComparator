## download and unzip simulation data (Bioanno Solve output)
## Google Drive: https://drive.google.com/file/d/1v59HRpO9dofnMQh5whVLEPdbm5QCGywX/view?usp=sharing

## Modify the absolute paths
sim_data=/Users/zhaoyangjia/PyCharm_Repos/KarCheck_refactor/manuscript_bionano_data/
omkar_exe_dir=/Users/zhaoyangjia/PyCharm_Repos/KarCheck_refactor/OMKar/

#################################################################################
build_version=${1:-default_build}
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pipeline_output=${SCRIPT_PATH}/simulation_karcheck/${build_version}/
version_register=${SCRIPT_PATH}/simulation_karcheck/omkar_versions.txt

mkdir -p ${pipeline_output}
current_dir=$(pwd)
cd "$omkar_exe_dir"
commit_code=$(git log -1 --format=%h)
cd "$current_dir"
echo "${build_version}    ${commit_code}" >> ${version_register}


## Run OMKar Batch Mode (activate corresponding virtual/conda env if needed)
## HTML report is not used in the comparison
echo "---------------RUNNING OMKAR-------------------"
python3 ${omkar_exe_dir}/main.py -dir ${sim_data} -o ${pipeline_output}/omkar_outs/

# Extract OMKar output for comparison
echo "---------------EXTRACTING OMKAR-------------------"
data_dir=${pipeline_output}/omkar_outs/omkar_output/
mkdir -p ${pipeline_output}/omkar_paths/
for folder in ${data_dir}/*
do
  folder_name="$(basename "$folder")"
	cp ${folder}/${folder_name}_raw.txt ${pipeline_output}/omkar_paths/${folder_name}.txt
done

## Clustering each genome into independent chromosomal clusters
echo "---------------CLUSTERING-------------------"
cmd="python ${SCRIPT_PATH}/omkar_analyses_pipeline/batch_clustering.py \
  --omkar_output ${pipeline_output}/omkar_outs/omkar_output/ \
  --omkar_paths_output ${pipeline_output}/omkar_paths/ \
  --output_dir ${pipeline_output}/cluster_files/ \
  --karsim_dir ${SCRIPT_PATH}/manuscript_karsim_data/karsim/"
echo "$cmd"
eval "$cmd"

## KarCheck
echo "---------------COMPARING-------------------"
cmd="python ${SCRIPT_PATH}/KarCheck.py \
  --matching_allowance 200000 \
  --approx_allowance 50000 \
  --karsim_kt_dir ${SCRIPT_PATH}/manuscript_karsim_data/karsim/ \
  --karsim_edgelabel_dir ${SCRIPT_PATH}/manuscript_karsim_data/karsim_edge_labels/ \
  --metadata_dir ${SCRIPT_PATH}/KarUtils/Metadata/ \
  --output_dir ${pipeline_output}/karcheck/ \
  --clusters_dir ${pipeline_output}/cluster_files/ \
  --omkar_output_dir ${pipeline_output}/omkar_outs/omkar_output \
  --simulation_bionano_output ${sim_data}"
echo "$cmd"
eval "$cmd"

