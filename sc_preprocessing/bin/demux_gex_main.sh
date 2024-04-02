#!/bin/bash

#SBATCH --cpus-per-task 50
#SBATCH --mem 1800G

source gex_functions.sh
module load bcl2fastq

# Initialize variables that will be set by command-line arguments
SEQ_NAME=""
CORES=${SLURM_CPUS_PER_TASK}
SEQ_TYPE=""
TENx_TYPE=""

# Parse command-line options
while getopts 's:c:t:x:' flag; do
  case "${flag}" in
    s) SEQ_NAME="${OPTARG}" ;;
    t) SEQ_TYPE="${OPTARG}" ;;
    x) TENx_TYPE="${OPTARG}" ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

# Check for mandatory options
if [ -z "$SEQ_NAME" ] || [ -z "$SEQ_TYPE" ] || [ -z "$TENx_TYPE" ]; then
  echo "Usage: $0 -s SEQ_NAME -t SEQ_TYPE (mouse, human, tapseq) -x TENx_TYPE (3.0 v 3.1)"
  exit 1
fi

# Set paths based on SEQ_NAME
FASTQ_DIR="/projects/perslab/people/lhv464/data/sc-10x/data-mkfastq/200123_SCOP-37/${SEQ_NAME}"
OUTPATH="/projects/perslab/people/lhv464/glp1r_leptin/sc_preprocessing"
BCL_DIR="/projects/perslab/people/lhv464/data/bcl/200123_SCOP-37"
FLOWCELL_DIR=$(ls -d ${BCL_DIR}/* | grep "$SEQ_NAME")
SAMPLE_SHEET_PATH="${OUTPATH}/sample_sheets/${SEQ_NAME}_samplesheet.csv"

# Define the path to the specific version of Cell Ranger
CELLRANGER_PATH="/projects/perslab/people/lhv464/glp1r_leptin/sc_preprocessing/cellranger-7.0.1"
# Define the Conda environment name
CONDA_ENV_NAME="kbpython_update"

# Prepend the Cell Ranger path to the existing PATH variable
# This ensures that the specific version of Cell Ranger is found first
export PATH="$CELLRANGER_PATH:$PATH"


# Check if BCL conversion is needed and run it
if [ ! -d "${FASTQ_DIR}/Reports/" ]; then
  run_bcl2fastq "$FLOWCELL_DIR" "$FASTQ_DIR" "$SAMPLE_SHEET_PATH" "$TENx_TYPE"
fi

# Iterate through each subdirectory in FASTQ_DIR
for dir in "${FASTQ_DIR}"/P*; do
  subdir=$(basename "$dir")
  seq_lane=$(echo "$subdir" | sed 's/^.//; s/[E-].*//')
  feature_file=("${OUTPATH}/hto_feature_files/*P${seq_lane}*")

  # Check if Cell Ranger count is needed and run it
  if [ ! -f "${OUTPATH}/gex_counts/${subdir}/outs/raw_feature_bc_matrix.h5" ]; then
    run_cellranger "${subdir}" "${FASTQ_DIR}/${subdir}" "${CORES}" "${SEQ_TYPE}" "${OUTPATH}"
  fi

  # Check if Kite count is needed and run it
  hto_dir="${OUTPATH}/hto_counts/${subdir}"
  if [ ! -f "${hto_dir}/kb_info.json" ]; then
    # Initialize Conda for bash shell
    eval "$(conda shell.bash hook)"
    #Activate the Conda environment
    conda activate "$CONDA_ENV_NAME"
    # Check if Conda environment was successfully activated
    if [ $? -eq 0 ]; then
      echo "Conda environment '$CONDA_ENV_NAME' activated."
      run_kite "${FASTQ_DIR}/${subdir}" "${hto_dir}" "${feature_file}" "${CORES}"
    else
      echo "Failed to activate Conda environment '$CONDA_ENV_NAME'."
      exit 1
    fi
  fi

done
