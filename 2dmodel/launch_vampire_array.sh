#!/bin/bash

# ----------- User Configuration -------------
START=0.05
STOP=1.55
STEP=0.05
OUTPUT_DIR="output-Field8-0T-Ba_high_res_steps10000x10000"
MAX_PARALLEL=10              # Max concurrent jobs (e.g. 10 nodes)
SBATCH_TIME="23:00:00"       # Time per job
SBATCH_JOBNAME="vampire_array"
SBATCH_CORES=48
SBATCH_NODES=1
SBATCH_SCRIPT_NAME="submit_array.sbatch"
# --------------------------------------------

# ----------- Compute Array Size -------------
TOTAL_STEPS=$(python3 -c "print(int(round((${STOP} - ${START}) / ${STEP})))")
ARRAY_MAX=$((TOTAL_STEPS))
echo "[+] Detected ${TOTAL_STEPS} temperature steps (array: 0 to $ARRAY_MAX)"
# --------------------------------------------

PAIR_FILE=$(grep "material:unit-cell-file" input | sed -E 's/.*=\s*"?([^"]+)"?/\1/')
mkdir -p $OUTPUT_DIR
cp input vampire*.mat $PAIR_FILE $OUTPUT_DIR

# ----------- Generate SLURM Script ----------
cat > "$SBATCH_SCRIPT_NAME" <<EOF
#!/bin/bash
#SBATCH --job-name=${SBATCH_JOBNAME}
#SBATCH --array=0-${ARRAY_MAX}%${MAX_PARALLEL}
#SBATCH --nodes=${SBATCH_NODES}
#SBATCH --ntasks=${SBATCH_CORES}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --output=slurm-%A_%a.out

# source /opt/intel/oneapi/setvars.sh
# PATH=\$PATH:\$HOME/codes/vampire

SECONDS=0

START=${START}
STEP=${STEP}
OUTPUT_DIR=${OUTPUT_DIR}
INPUT_TEMPLATE="input"
SCRIPT_DIR=\$(pwd)

# ---- Compute current temperature ----
TEMP=\$(python3 -c "print(round(\${START} + \${SLURM_ARRAY_TASK_ID} * \${STEP}, 5))")
TEMP_STR=\$(printf "%.3f" \$TEMP)
RUN_DIR="\$OUTPUT_DIR/run_T_\${TEMP_STR}"

echo "[\$TEMP_STR] Setting up run in \$RUN_DIR : \$(date)"

mkdir -p "\$RUN_DIR" "\$OUTPUT_DIR"
cp "\$OUTPUT_DIR/\$INPUT_TEMPLATE" "\$RUN_DIR"/
cp \$OUTPUT_DIR/vampire*.mat "\$RUN_DIR"/

sed -i -E "s/^(sim:temperature\s*=\s*).*/\1\$TEMP/" "\$RUN_DIR/input"

cd "\$RUN_DIR"
mpirun -np ${SBATCH_CORES} vampire-parallel > vampire.log
cd "\$SCRIPT_DIR"

cp "\$RUN_DIR/output" "\$OUTPUT_DIR/output_\${TEMP_STR}"
cp "\$RUN_DIR/log" "\$OUTPUT_DIR/"
cp "\$RUN_DIR/vampire.log" "\$OUTPUT_DIR/"
rm -rf "\$RUN_DIR"

echo
echo "[\$TEMP_STR] Done and cleaned up : \$(date)"

# Format SECONDS into hh:mm:ss
hours=\$((SECONDS / 3600))
minutes=\$(((SECONDS % 3600) / 60))
seconds=\$((SECONDS % 60))
printf -v formatted "%02d:%02d:%02d" "\$hours" "\$minutes" "\$seconds"
echo "Elapsed time: \$formatted"
EOF
# --------------------------------------------

# ----------- Submit Job ---------------------
echo "[+] Submitting job array via sbatch"
sbatch "$SBATCH_SCRIPT_NAME"
