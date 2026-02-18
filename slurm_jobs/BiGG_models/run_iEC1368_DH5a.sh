#!/bin/bash
#SBATCH --account=def-mahadeva
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=iEC1368_DH5a
#SBATCH --output=logs/BiGG_models/iEC1368_DH5a_%j.out
#SBATCH --error=logs/BiGG_models/iEC1368_DH5a_%j.err

# Load modules
module load StdEnv/2023
module load python/3.11

# Navigate to project directory
cd $SLURM_SUBMIT_DIR

# Activate virtual environment
source venv/bin/activate

# Create logs directory
mkdir -p logs/BiGG_models

# Run pipeline
echo "Starting kinGEMs pipeline for iEC1368_DH5a"
echo "Config: configs/BiGG_models/iEC1368_DH5a.json"
echo "Time: $(date)"
echo ""

python scripts/run_pipeline.py configs/BiGG_models/iEC1368_DH5a.json

echo ""
echo "Completed at: $(date)"
