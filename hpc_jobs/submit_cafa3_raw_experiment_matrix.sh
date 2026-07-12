#!/bin/bash
# Submit four new raw-GOA CAFA3 policy variants. The completed 7061922 baseline
# is intentionally not resubmitted; SGE decides when each independent job runs.
set -euo pipefail

WRAPPER="hpc_jobs/active/hpc_cafa3_historical_validation.sh"
[ -f "$WRAPPER" ] || { echo "Run this script from the repository root" >&2; exit 1; }

COMMON="HISTORICAL_TRAINING_SNAPSHOT=september-2016,TARGET_UNIVERSE_POLICY=official-cafa3-targets,HISTORICAL_TEST_SOURCE=raw-goa"

qsub -N c3_snap_bf_feb -v "${COMMON},HISTORICAL_T1_ENDPOINT_POLICY=snapshot-membership,HISTORICAL_BACKFILL_POLICY=exclude-pre-t0,HISTORICAL_BENCHMARK_ONTOLOGY=february-go-basic" "$WRAPPER"
qsub -N c3_snap_nobf_feb -v "${COMMON},HISTORICAL_T1_ENDPOINT_POLICY=snapshot-membership,HISTORICAL_BACKFILL_POLICY=allow,HISTORICAL_BENCHMARK_ONTOLOGY=february-go-basic" "$WRAPPER"
qsub -N c3_snap_bf_dgp -v "${COMMON},HISTORICAL_T1_ENDPOINT_POLICY=snapshot-membership,HISTORICAL_BACKFILL_POLICY=exclude-pre-t0,HISTORICAL_BENCHMARK_ONTOLOGY=deepgoplus-packaged" "$WRAPPER"
qsub -N c3_date_bf_dgp -v "${COMMON},HISTORICAL_T1_ENDPOINT_POLICY=assigned-date-proxy,HISTORICAL_BACKFILL_POLICY=exclude-pre-t0,HISTORICAL_BENCHMARK_ONTOLOGY=deepgoplus-packaged" "$WRAPPER"
