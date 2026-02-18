#!/bin/bash
#SBATCH --account=def-mahadeva
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=iND750
#SBATCH --output=logs/BiGG_models/iND750_%j.out
#SBATCH --error=logs/BiGG_models/iND750_%j.err

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
echo "Starting kinGEMs pipeline for iND750"
echo "Config: configs/BiGG_models/iND750.json"
echo "Time: $(date)"
echo ""

python scripts/run_pipeline.py configs/BiGG_models/iND750.json

echo ""
echo "Completed at: $(date)"
