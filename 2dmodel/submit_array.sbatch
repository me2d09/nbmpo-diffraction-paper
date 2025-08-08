#!/bin/bash
#SBATCH --job-name=vampire_array
#SBATCH --array=0-30%10
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=23:00:00
#SBATCH --output=slurm-%A_%a.out

# source /opt/intel/oneapi/setvars.sh
# PATH=$PATH:$HOME/codes/vampire

SECONDS=0

START=0.05
STEP=0.05
OUTPUT_DIR=output-Field8-0T-Ba_high_res_steps10000x10000
INPUT_TEMPLATE="input"
SCRIPT_DIR=$(pwd)

# ---- Compute current temperature ----
TEMP=$(python3 -c "print(round(${START} + ${SLURM_ARRAY_TASK_ID} * ${STEP}, 5))")
TEMP_STR=$(printf "%.3f" $TEMP)
RUN_DIR="$OUTPUT_DIR/run_T_${TEMP_STR}"

echo "[$TEMP_STR] Setting up run in $RUN_DIR : $(date)"

mkdir -p "$RUN_DIR" "$OUTPUT_DIR"
cp "$OUTPUT_DIR/$INPUT_TEMPLATE" "$RUN_DIR"/
cp $OUTPUT_DIR/vampire*.mat "$RUN_DIR"/

sed -i -E "s/^(sim:temperature\s*=\s*).*/\1$TEMP/" "$RUN_DIR/input"

cd "$RUN_DIR"
mpirun -np 48 vampire-parallel > vampire.log
cd "$SCRIPT_DIR"

cp "$RUN_DIR/output" "$OUTPUT_DIR/output_${TEMP_STR}"
cp "$RUN_DIR/log" "$OUTPUT_DIR/"
cp "$RUN_DIR/vampire.log" "$OUTPUT_DIR/"
rm -rf "$RUN_DIR"

echo
echo "[$TEMP_STR] Done and cleaned up : $(date)"

# Format SECONDS into hh:mm:ss
hours=$((SECONDS / 3600))
minutes=$(((SECONDS % 3600) / 60))
seconds=$((SECONDS % 60))
printf -v formatted "%02d:%02d:%02d" "$hours" "$minutes" "$seconds"
echo "Elapsed time: $formatted"
