#!/usr/bin/env python3
"""
Generate SLURM job scripts for all BiGG models.
"""

import os
from pathlib import Path

def create_slurm_job(model_name, config_path):
    """Create a SLURM job script for a specific model."""

    slurm_script = f'''#!/bin/bash
#SBATCH --account=def-mahadeva
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name={model_name}
#SBATCH --output=logs/BiGG_models/{model_name}_%j.out
#SBATCH --error=logs/BiGG_models/{model_name}_%j.err

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
echo "Starting kinGEMs pipeline for {model_name}"
echo "Config: {config_path}"
echo "Time: $(date)"
echo ""

python scripts/run_pipeline.py {config_path}

echo ""
echo "Completed at: $(date)"
'''

    return slurm_script


def main():
    """Generate SLURM job scripts for all BiGG models."""

    project_root = Path(__file__).parent.parent
    configs_dir = project_root / "configs" / "BiGG_models"
    slurm_jobs_dir = project_root / "slurm_jobs" / "BiGG_models"

    # Create slurm jobs directory
    slurm_jobs_dir.mkdir(parents=True, exist_ok=True)

    # Get all config files
    config_files = sorted(configs_dir.glob("*.json"))

    print(f"Generating SLURM job scripts for {len(config_files)} BiGG models")
    print(f"Output directory: {slurm_jobs_dir}\n")

    generated = 0

    for config_file in config_files:
        model_name = config_file.stem
        config_path = f"configs/BiGG_models/{config_file.name}"

        # Create SLURM script
        slurm_script = create_slurm_job(model_name, config_path)

        # Save to file
        output_path = slurm_jobs_dir / f"run_{model_name}.sh"
        with open(output_path, 'w') as f:
            f.write(slurm_script)

        # Make executable
        os.chmod(output_path, 0o755)

        generated += 1
        if generated <= 5 or generated % 20 == 0:
            print(f"  Generated: {output_path.name}")

    print(f"\n✓ Generated {generated} SLURM job scripts in {slurm_jobs_dir}")
    print(f"\nTo submit all jobs, run:")
    print(f"  bash slurm_jobs/submit_all_bigg_models.sh")


if __name__ == '__main__':
    main()
