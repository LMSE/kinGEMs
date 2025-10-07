#!/bin/bash
#
# kinGEMs Pipeline Runner with Gurobi Module Loading
#
# This script loads the required Gurobi module and then runs the pipeline.
# Use this instead of calling run_pipeline.py directly on HPC systems.
#
# Usage:
#   bash scripts/run_with_gurobi.sh configs/iML1515_GEM.json
#   bash scripts/run_with_gurobi.sh configs/382_genome_cpd03198.json --force
#

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PIPELINE_SCRIPT="$SCRIPT_DIR/run_pipeline.py"

# Color output
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  kinGEMs Pipeline Runner (HPC)${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if config file is provided
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <config_file> [--force]"
    echo ""
    echo "Examples:"
    echo "  $0 configs/iML1515_GEM.json"
    echo "  $0 configs/382_genome_cpd03198.json --force"
    exit 1
fi

# Load Gurobi module (Compute Canada)
echo -e "${YELLOW}Loading Gurobi module...${NC}"
if module load gurobi 2>/dev/null; then
    echo -e "${GREEN}✓ Gurobi module loaded successfully${NC}"
else
    echo -e "${YELLOW}⚠ Warning: Could not load gurobi module${NC}"
    echo -e "${YELLOW}  This may be expected if not on Compute Canada cluster${NC}"
    echo -e "${YELLOW}  Or if Gurobi is already in PATH${NC}"
fi
echo ""

# Check if Gurobi is available
if command -v gurobi_cl &> /dev/null; then
    GUROBI_VERSION=$(gurobi_cl --version 2>&1 | head -n 1 || echo "unknown")
    echo -e "${GREEN}✓ Gurobi found: $GUROBI_VERSION${NC}"
else
    echo -e "${YELLOW}⚠ Warning: gurobi_cl not found in PATH${NC}"
    echo -e "${YELLOW}  The pipeline may fail if Gurobi is required${NC}"
fi
echo ""

# Activate virtual environment if it exists
VENV_PATH="$PROJECT_ROOT/venv"
if [[ -d "$VENV_PATH" ]]; then
    echo -e "${YELLOW}Activating virtual environment...${NC}"
    source "$VENV_PATH/bin/activate"
    echo -e "${GREEN}✓ Virtual environment activated${NC}"
    echo ""
fi

# Run the pipeline
echo -e "${BLUE}Running pipeline with config: $1${NC}"
if [[ " $@ " =~ " --force " ]] || [[ " $@ " =~ " -f " ]]; then
    echo -e "${YELLOW}Mode: Force regeneration${NC}"
fi
echo ""

python "$PIPELINE_SCRIPT" "$@"

exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    echo ""
    echo -e "${GREEN}✓ Pipeline completed successfully${NC}"
else
    echo ""
    echo -e "${YELLOW}⚠ Pipeline exited with code $exit_code${NC}"
fi

exit $exit_code
