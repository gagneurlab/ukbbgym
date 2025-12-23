#!/bin/bash

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="logs/regenie_prep_${TIMESTAMP}.log"

nohup python -u -m prep_files_for_REGENIE2.py > "$LOGFILE" 2>&1 &

echo "Job started. Log file: $LOGFILE"