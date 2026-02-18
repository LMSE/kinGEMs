#!/bin/bash
#
# Run kinGEMs pipeline for all BiGG models
# Usage: bash scripts/run_all_bigg_models.sh
#

# Activate virtual environment
source venv/bin/activate

# Directory containing config files
CONFIG_DIR="configs/BiGG_models"

# Log directory
LOG_DIR="logs/BiGG_models_runs"
mkdir -p "$LOG_DIR"

# Get total number of config files
TOTAL=$(ls $CONFIG_DIR/*.json 2>/dev/null | wc -l)

if [ $TOTAL -eq 0 ]; then
    echo "No config files found in $CONFIG_DIR"
    exit 1
fi

echo "======================================================================"
echo "Running kinGEMs pipeline for $TOTAL BiGG models"
echo "Config directory: $CONFIG_DIR"
echo "Log directory: $LOG_DIR"
echo "======================================================================"
echo ""

COMPLETED=0
FAILED=0
SKIPPED=0

# Process each config file
for config_file in $CONFIG_DIR/*.json; do
    # Extract model name
    model_name=$(basename "$config_file" .json)

    echo "[$((COMPLETED + FAILED + SKIPPED + 1))/$TOTAL] Processing $model_name..."

    # Check if model XML exists
    model_xml="data/raw/BiGG_models/${model_name}.xml"
    if [ ! -f "$model_xml" ]; then
        echo "  ⚠️  Model file not found: $model_xml - SKIPPING"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Check if predictions exist
    predictions_found=0
    for pred_file in data/interim/CPI-Pred\ predictions/*${model_name}*.csv; do
        if [ -f "$pred_file" ]; then
            predictions_found=1
            break
        fi
    done

    if [ $predictions_found -eq 0 ]; then
        echo "  ⚠️  No CPI-Pred predictions found for $model_name - SKIPPING"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Run pipeline
    log_file="$LOG_DIR/${model_name}_$(date +%Y%m%d_%H%M%S).log"

    if python scripts/run_pipeline.py "$config_file" > "$log_file" 2>&1; then
        echo "  ✓ Completed successfully (log: $log_file)"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "  ✗ Failed (check log: $log_file)"
        FAILED=$((FAILED + 1))
    fi

    echo ""
done

echo "======================================================================"
echo "Batch processing complete!"
echo "  Completed: $COMPLETED"
echo "  Failed: $FAILED"
echo "  Skipped: $SKIPPED"
echo "  Total: $TOTAL"
echo "======================================================================"
