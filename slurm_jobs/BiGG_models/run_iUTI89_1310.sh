#!/bin/bash
#SBATCH --account=def-mahadeva
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=iUTI89_1310
#SBATCH --output=logs/BiGG_models/iUTI89_1310_%j.out
#SBATCH --error=logs/BiGG_models/iUTI89_1310_%j.err

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
echo "Starting kinGEMs pipeline for iUTI89_1310"
echo "Config: configs/BiGG_models/iUTI89_1310.json"
echo "Time: $(date)"
echo ""

python scripts/run_pipeline.py configs/BiGG_models/iUTI89_1310.json

echo ""
echo "Completed at: $(date)"
