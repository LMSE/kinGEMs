#!/bin/bash
#
# Batch Pipeline Runner for kinGEMs
#
# This script provides convenient batch processing of multiple GEMs.
#
# Usage:
#   bash scripts/batch_run.sh all              # Run all models
#   bash scripts/batch_run.sh ecoli            # Run all E. coli models
#   bash scripts/batch_run.sh genome           # Run all ModelSEED models
#   bash scripts/batch_run.sh standard         # Run all standard (non-genome) models
#   bash scripts/batch_run.sh <model1> <model2> ...  # Run specific models
#
# Options:
#   --force, -f   Force regeneration of intermediate files
#
# Examples:
#   bash scripts/batch_run.sh all --force
#   bash scripts/batch_run.sh iML1515_GEM yeast-GEM9
#   bash scripts/batch_run.sh ecoli

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIGS_DIR="$PROJECT_ROOT/configs"
PIPELINE_SCRIPT="$SCRIPT_DIR/run_pipeline.py"

# Parse force flag
FORCE_FLAG=""
if [[ " $@ " =~ " --force " ]] || [[ " $@ " =~ " -f " ]]; then
    FORCE_FLAG="--force"
    # Remove force flag from arguments
    set -- "${@/--force/}"
    set -- "${@/-f/}"
fi

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print header
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  kinGEMs Batch Pipeline Runner${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Model groups
ECOLI_MODELS=(
    "iML1515_GEM"
    "e_coli_core"
    "382_genome_cpd03198"
    "376_genome_cpd03198"
    "378_genome_cpd03198"
    "380_genome_cpd01262"
)

GENOME_MODELS=(
    "382_genome_cpd03198"
    "376_genome_cpd03198"
    "378_genome_cpd03198"
    "380_genome_cpd01262"
)

STANDARD_MODELS=(
    "iML1515_GEM"
    "e_coli_core"
    "yeast-GEM9"
    "Kmarxianus_GEM"
    "Lmajor_GEM"
    "Pputida_iJN1463"
    "Pputida_iJN746"
    "Selongatus_iJB785"
)

ALL_MODELS=(
    "${STANDARD_MODELS[@]}"
    "${GENOME_MODELS[@]}"
)

# Function to run a single model
run_model() {
    local model_name=$1
    local config_file="$CONFIGS_DIR/${model_name}.json"
    
    if [[ ! -f "$config_file" ]]; then
        echo -e "${RED}✗ Config not found: $config_file${NC}"
        return 1
    fi
    
    echo -e "${YELLOW}→ Processing: $model_name${NC}"
    echo -e "  Config: $config_file"
    if [[ -n "$FORCE_FLAG" ]]; then
        echo -e "  Mode: Force regeneration"
    else
        echo -e "  Mode: Use cached files"
    fi
    echo ""
    
    if python "$PIPELINE_SCRIPT" "$config_file" $FORCE_FLAG; then
        echo -e "${GREEN}✓ Success: $model_name${NC}"
        echo ""
        return 0
    else
        echo -e "${RED}✗ Failed: $model_name${NC}"
        echo ""
        return 1
    fi
}

# Parse arguments
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <preset|model_names...> [--force]"
    echo ""
    echo "Presets:"
    echo "  all        - Run all models"
    echo "  ecoli      - Run all E. coli models"
    echo "  genome     - Run all ModelSEED (genome) models"
    echo "  standard   - Run all standard (non-genome) models"
    echo ""
    echo "Options:"
    echo "  --force, -f  - Force regeneration of intermediate files"
    echo ""
    echo "Examples:"
    echo "  $0 all"
    echo "  $0 ecoli --force"
    echo "  $0 iML1515_GEM yeast-GEM9"
    exit 1
fi

# Determine which models to run
MODELS_TO_RUN=()

case "$1" in
    all)
        MODELS_TO_RUN=("${ALL_MODELS[@]}")
        echo -e "${BLUE}Running ALL models (${#MODELS_TO_RUN[@]} total)${NC}"
        ;;
    ecoli)
        MODELS_TO_RUN=("${ECOLI_MODELS[@]}")
        echo -e "${BLUE}Running E. coli models (${#MODELS_TO_RUN[@]} total)${NC}"
        ;;
    genome)
        MODELS_TO_RUN=("${GENOME_MODELS[@]}")
        echo -e "${BLUE}Running ModelSEED genome models (${#MODELS_TO_RUN[@]} total)${NC}"
        ;;
    standard)
        MODELS_TO_RUN=("${STANDARD_MODELS[@]}")
        echo -e "${BLUE}Running standard models (${#MODELS_TO_RUN[@]} total)${NC}"
        ;;
    *)
        # Custom list of models
        MODELS_TO_RUN=("$@")
        echo -e "${BLUE}Running custom model list (${#MODELS_TO_RUN[@]} total)${NC}"
        ;;
esac

echo ""

# Track results
SUCCESS_COUNT=0
FAIL_COUNT=0
FAILED_MODELS=()

# Run all models
for model in "${MODELS_TO_RUN[@]}"; do
    if run_model "$model"; then
        ((SUCCESS_COUNT++))
    else
        ((FAIL_COUNT++))
        FAILED_MODELS+=("$model")
    fi
done

# Print summary
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Batch Processing Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Successful: $SUCCESS_COUNT${NC}"
echo -e "${RED}Failed: $FAIL_COUNT${NC}"

if [[ ${#FAILED_MODELS[@]} -gt 0 ]]; then
    echo ""
    echo -e "${RED}Failed models:${NC}"
    for model in "${FAILED_MODELS[@]}"; do
        echo -e "  ${RED}✗ $model${NC}"
    done
    exit 1
else
    echo -e "${GREEN}All models processed successfully!${NC}"
    exit 0
fi
