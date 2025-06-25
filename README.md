# KarComparator

KarCheck is designed to be a symmetric comparator between two unphased karyotypes, by comparing their SV and CNV calls. In addition, to accommodate molecular
techniques that do not have a nucleotide-level resolution, KarCheck has an adjustable tolerance for
small breakpoint mis-matching for an SV that is present in both karyotypes.

Currently, pipeline for comparing Simulation and OMKar's reconstruction are readily available. Pipeline for direct comparison between
two Molecular Karyotype is available upon request.

Download via

`git clone --recurse-submodules https://github.com/MolecularKaryotype/KarComparator`

## Pipeline: Running OMKar on simulated Bionano OGM, and then run KarCheck for comparison scores
This pipeline was used for generating the simulation results in *OMKar: optical map based automated karyotyping of genomes to identify constitutional abnormalities* (https://www.medrxiv.org/content/10.1101/2025.02.13.25322211v1)
### Running the Pipeline
`bash omkar_simulation_karcheck.sh [build_version]`

For example, `bash omkar_simulation_karcheck.sh b14`

Follow these steps to run the pipeline:
1) Download Simulated Data (Bionano Solve output) from Google Drive: https://drive.google.com/file/d/1v59HRpO9dofnMQh5whVLEPdbm5QCGywX/view?usp=sharing
2) Download latest OMKar from https://github.com/siavashre/OMKar or by running `git clone --recurse-submodules https://github.com/siavashre/OMKar.git`
3) Fill in the absolute path to your unzipped simulated data DIR and OMKar DIR in `omkar_simulation_karcheck.sh`
4) Run `bash omkar_simulation_karcheck.sh`

`omkar_simulation_karcheck.sh` accepts an optional parameter [build_version] to specify the output directory. All of its output will be in `simulation_karcheck/[build_version]/`.

KarCheck results can be found in `simulation_karcheck/[build_version]/karcheck/summary_stats.txt`.
Tabular results used to generate the summary statistics is in `simulation_karcheck/[build_version]/karcheck/whole_analysis_table_with_cn_005.csv`.

## Tool: Generate Molecular Karyotype file from partial cytogenetic annotation
To allow faster work flow in generating potential phased Molecular Karyotype
file, we developed a tool that assumes all un-annotated paths are WT and diploid.

The real data results in *OMKar: optical map based automated karyotyping of genomes to identify constitutional abnormalities* (https://www.medrxiv.org/content/10.1101/2025.02.13.25322211v1) were assisted by this function.
The real human samples' Molecular Karyotype file were generated using this function and expert cytogeneticists' annotation as input.

Please follow [Example Partial Annotation](sample_input/create_MK_file.txt) as input file format.
Here are some requirements:
- All cases in the batch run should be written in the same partial annotation file
- Each case is required to have
  1) a case name line begin with ">"
  2) #Sex: documenting all sex chromosomes, regardless of whether they are mentioned later
  3) #Segments: segments used later in the next section; use pter and qter to save time (the tool will automatically fill it out for you)
  4) #Paths: for a single chromosome that is affected, all of its homologous chromosomes need to be mentioned here (eg. if there is a deletion you mentioned on one path, you must mention the other homologous chromosomes, otherwise, we assume the other chromosomes is deleted; similarly, if there is a trisomy, mention all three chromosomes). 
  If a path is WT, you can omit the segment list in the path

To run it
`python Molecular_Karyotype.py partial_annotation [input_file] [output_dir]`


| Argument     | Type | Description                           |
|--------------|------|---------------------------------------|
| `input_file` | FILE | file path to the input file           |
| `output_dir` | DIR  | dir path to store the output MK files |



## Pipeline: Compare cases between two runs
**(Pending Update)**
This pipeline compares the cluster-by-cluster accuracy between two OMKar builds to help debug newly arised errors. 
For each cluster, it checks the new errors and resolved errors between builds.

### Running the Pipeline
`python compare_two_runs.py [version1] [version2]`

For example, `python compare_two_runs.py b12 b14`

### Outputs
Each run's output will be in `KarComparator/omkar_analyses_pipeline/comparisons`

1. CSV of cases where "n_different_path" changed between the versions: `<version1>vs<version2>_path_diffs.csv`
   1. "n_diff_path" = #OMKar_Path - #KarSim_Path
   2. For example, if n_path_diff_df1 = 1 and n_path_diff_df2 = 2, 
   3. this means we used to report 1 additional path in version1, 
   4. but we are now reporting 2 additional paths in version2
2. CSV of cases where number of SV missed changed between the versions: `<version1>vs<version2>_path_diffs.csv`
   1. "SV_missed_delta" = #Version2_SV_missed - #Version1_SV_missed

## Tool: Visualization of preILP/postILP graph for OMKar
Visualize the preILP and postILP graphs of OMKar runtime intermediate output for debugging.

Use under `KarComparator/` as working directory

`python debug_omkar.py [data_dir] [case_name] [output_dir] [chr_of_interest]`

| Argument     | Type | Description                                                                                           |
|--------------|------|-------------------------------------------------------------------------------------------------------|
| `data_dir`   | DIR  | OMKar output DIR                                                                                      |
| `case_name`  | STR  | basename of the case                                                                                  |
| `output_dir` | DIR  | DIR to store the output images                                                                        |
| `chr_of_interest`             | STR  | chromosomes to be graphed (must be all the chromosomes, and from the same cluster). e.g. "['3', '5']" |

sample script
`python debug_omkar.py 0717_output/ 1301 debug/ ['3','5']`

## Implementation Details of KarCheck
INPUT: cluster files
OUTPUT: cluster-unit statistics with 1) edge-based SV, 2) bin-based CN, and 3) path diffs


- `prep_df`
  - prepares the cluster file paths
  - check the origin chromosomes of each cluster
  - note the number of paths in each cluster
  - note whether a cluster is an event cluster, based on the karsim history log
  - document event counts for that cluster
- `KarCheck`
  - Two identical graphs are form according to the discription in the paper. One for SV use and one for CNV use.
  - Graph's prefix/suffix forbidden regions are ignored
  - Source and Sink nodes are added and edge land in the prefix/suffix forbidden regions are transfered to the start/end node, with an edge from the source/sink
  - SV
    - Event's unique IDs (relative to each cluster) are formed from KarSim history log; terminal events' IDs are noted 
    - Graph is pruned, approximated, and matched
      - matching distances are returned for downstream analyses (using Event Status Object)
      - after approximation, all residual KarSim edges are significant, and they are the edges we try to capture; 
      we filter KarSim's significant edges from the event-ids, such that only the nT-events' edges are considered for capturing
    - Final statistics (jaccard, recall, precision) computed based on nT-edges in KarSim, and OMKar's edges with proximity to the terminal are not counted as FP
  - Event Matching Distances
    - This is only analyzed for fully captured event (that all edges associated with that event-ids are matched)
    - distances are collected in terms of the endpoint distance for each edge (majority of edges will contribute to two entries in the distance histogram)
    - each edge is filtered with the SMAP such that we try to see if both endpoint of that edge's SV types appear in the SMAP
    - the distances histogram plot distances with edge present in SMAP/not present in SMAP; the ones not present in SMAP, OMKar was not provided the precise coordinate, so we are not responsible
  - CNV
    - genome's whole region is binned into 50kbp bins (excluding prefix/suffix forbidden regions)
    - each chromosome in the cluster is determined for a WT count = rounded avg CN from KarSim
    - A bin is considered "with variation" if it deviates from the expected count by more than 0.05
    - Jaccard score is computed among the bins (KarSim bins vs. OMKar bins), where a boolean value is assigned to each "with variation" bin, to see if any variation was included

