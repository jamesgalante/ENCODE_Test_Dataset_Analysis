#!/bin/bash
#SBATCH --job-name=dnase_cre_processing
#SBATCH --output=dnase_processing_%j.log
#SBATCH --error=dnase_processing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --partition=normal

# Activate the conda environment
echo "Initializing conda..."
# source /path/to/conda.sh
echo "Activating conda environment..."
conda activate dnase_processing

# Make sure we're in the right directory
cd $SLURM_SUBMIT_DIR

# Run the processing script
echo "Starting DNase-seq processing pipeline..."
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

# Execute the processing script
bash creating_candidate_cres_bed_file.sh

echo "Pipeline completed!"
echo "End time: $(date)"

# Deactivate conda environment
source deactivate