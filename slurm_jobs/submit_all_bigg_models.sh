#!/bin/bash
# =====================================================================
# Master submission script for all BiGG models
# =====================================================================
# This script submits SLURM jobs for all 108 BiGG models
#
# Usage:
#   ./slurm_jobs/submit_all_bigg_models.sh          # Submit all models
#   ./slurm_jobs/submit_all_bigg_models.sh --dry-run  # Preview without submitting
#
# Options:
#   --dry-run    : Show what would be submitted without actually submitting
#   --batch N    : Submit only N jobs at a time (default: submit all)
#   --organism X : Submit only models containing organism name X (e.g., "coli")
# =====================================================================

DRY_RUN=false
BATCH_SIZE=0
ORGANISM_FILTER=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --batch)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --organism)
            ORGANISM_FILTER="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--batch N] [--organism NAME]"
            exit 1
            ;;
    esac
done

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIGG_DIR="$SCRIPT_DIR/BiGG_models"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "========================================"
echo "BiGG Models SLURM Job Submission"
echo "========================================"
echo "Date: $(date)"
echo "Script directory: $SCRIPT_DIR"
echo "BiGG jobs directory: $BIGG_DIR"
echo ""

# Count total job scripts
TOTAL_JOBS=$(ls -1 "$BIGG_DIR"/run_*.sh 2>/dev/null | wc -l)

if [ $TOTAL_JOBS -eq 0 ]; then
    echo "Error: No job scripts found in $BIGG_DIR"
    exit 1
fi

echo "Found $TOTAL_JOBS BiGG model job scripts"

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN MODE - No jobs will be submitted"
fi

if [ $BATCH_SIZE -gt 0 ]; then
    echo "Batch mode: Submitting $BATCH_SIZE jobs at a time"
fi

if [ -n "$ORGANISM_FILTER" ]; then
    echo "Organism filter: '$ORGANISM_FILTER'"
fi

echo ""

# Submit jobs
SUBMITTED=0
SKIPPED=0
FAILED=0

for script in "$BIGG_DIR"/run_*.sh; do
    MODEL_NAME=$(basename "$script" .sh | sed 's/run_//')

    # Apply organism filter if specified
    if [ -n "$ORGANISM_FILTER" ]; then
        if ! echo "$MODEL_NAME" | grep -qi "$ORGANISM_FILTER"; then
            ((SKIPPED++))
            continue
        fi
    fi

    # Apply batch size limit if specified
    if [ $BATCH_SIZE -gt 0 ] && [ $SUBMITTED -ge $BATCH_SIZE ]; then
        echo "Batch limit reached ($BATCH_SIZE jobs). Stopping."
        break
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would submit: $MODEL_NAME"
        ((SUBMITTED++))
    else
        echo -n "Submitting job for model: $MODEL_NAME ... "
        JOB_ID=$(sbatch "$script" 2>&1)

        if [ $? -eq 0 ]; then
            # Extract job ID from output
            JOB_NUM=$(echo "$JOB_ID" | grep -oE '[0-9]+' | head -1)
            echo "OK (Job ID: $JOB_NUM)"
            ((SUBMITTED++))

            # Small delay to avoid overwhelming scheduler
            sleep 0.2
        else
            echo "FAILED"
            echo "  Error: $JOB_ID"
            ((FAILED++))
        fi
    fi
done

echo ""
echo "========================================"
echo "Submission Summary"
echo "========================================"
echo "Total job scripts: $TOTAL_JOBS"
echo "Jobs submitted: $SUBMITTED"
if [ $SKIPPED -gt 0 ]; then
    echo "Jobs skipped (filter): $SKIPPED"
fi
if [ $FAILED -gt 0 ]; then
    echo "Jobs failed: $FAILED"
fi
echo ""

if [ "$DRY_RUN" = false ] && [ $SUBMITTED -gt 0 ]; then
    echo "Monitor jobs with: squeue -u \$USER"
    echo "Check logs in: logs/bigg_*"
    echo ""
fi

echo "Done!"
