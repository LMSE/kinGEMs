#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for Post-tuning kinGEMs validation (enzyme-constrained)
# ---------------------------------------------------------------------
#SBATCH --job-name=posttuning_kinGEMs
#SBATCH --account=def-mahadeva
#SBATCH --cpus-per-task=4
#SBATCH --time=0-15:00:00
#SBATCH --mem=60G
#SBATCH --output=logs/posttuning_%j.out
#SBATCH --error=logs/posttuning_%j.err
#SBATCH --mail-user=ranamoneim@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
# ---------------------------------------------------------------------
echo "========================================="
echo "Job: Post-tuning kinGEMs Validation"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "========================================="
# ---------------------------------------------------------------------

# Load required modules
module load python/3.11

# Optional: Load CPLEX for faster optimization (10-20x speedup)
# Uncomment if CPLEX is available on your cluster:
# module load cplex/22.1.1

# Activate virtual environment
source venv/bin/activate

# Create logs directory if it doesn't exist
mkdir -p logs results/validation_parallel

# Run post-tuning validation only
echo "Running Post-tuning kinGEMs validation..."
python scripts/run_validation_parallel.py \
    --mode posttuning \
    --config configs/validation_iML1515.json \
    --output results/validation_parallel

exitcode=$?

# Monitor memory usage
echo ""
echo "Memory Usage Statistics:"
if command -v sacct &> /dev/null; then
    sacct -j $SLURM_JOB_ID --format=JobID,MaxRSS,AveRSS,MaxVMSize --noheader
    PEAK_MEM=$(sacct -j $SLURM_JOB_ID --format=MaxRSS --noheader | sort -n | tail -1)
    echo "Peak memory usage: $PEAK_MEM"
else
    echo "  sacct not available - memory stats unavailable"
fi

# ---------------------------------------------------------------------
echo "========================================="
echo "Finished at: `date`"
echo "Exit code: $exitcode"
echo "========================================="

exit $exitcode
