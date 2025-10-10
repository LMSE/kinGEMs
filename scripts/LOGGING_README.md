# Running the kinGEMs Pipeline with Logging

## Quick Start

Instead of running the pipeline directly, use the wrapper script to automatically save logs:

```bash
# Run with logging
./scripts/run_pipeline_with_logging.sh configs/iML1515_GEM.json

# Run with force regeneration
./scripts/run_pipeline_with_logging.sh configs/iML1515_GEM.json --force
```

## What it does

The wrapper script (`run_pipeline_with_logging.sh`):

1. Creates the `results/logs/` directory if it doesn't exist
2. Runs `run_pipeline.py` with all your arguments
3. Displays output in real-time on your terminal (using `tee`)
4. Saves all output to a log file
5. Names the log file with the same run ID as the tuning results folder

## Output Location

Log files are saved as:
```
results/logs/{run_id}.out
```

For example, if your tuning results are in:
```
results/tuning_results/ecoli_iML1515_20250826_4941/
```

Then your log will be:
```
results/logs/ecoli_iML1515_20250826_4941.out
```

## Original Script

You can still run the original script directly without logging:

```bash
python scripts/run_pipeline.py configs/iML1515_GEM.json
```

## Viewing Logs Later

To view a log file:

```bash
# View entire log
cat results/logs/ecoli_iML1515_20250826_4941.out

# View last 100 lines
tail -100 results/logs/ecoli_iML1515_20250826_4941.out

# View with pagination
less results/logs/ecoli_iML1515_20250826_4941.out
```
