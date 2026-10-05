#!/usr/bin/env bash
# Launcher for the Comparative Gene & Promoter Pipeline.
# Uses the currently active Conda environment and preserves the user's
# current working directory.

set -euo pipefail
shopt -s nullglob

echo "=============================================================================="
echo "                   COMPARATIVE GENE & PROMOTER PIPELINE"
echo "=============================================================================="

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "[ERROR] No conda environment is active."
    echo "        Activate the environment where comparative-gene-promoter is installed."
    echo "        Example: conda activate comparative-gene-promoter-test"
    exit 1
fi

echo "[INFO] Conda environment : ${CONDA_PREFIX}"
echo "[INFO] Working directory : $(pwd)"
echo "       (input genome folders are read from here; outputs/ is written here)"

export COMPARATIVE_PIPELINE_WORKDIR="$(pwd)"

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -n "${COMPARATIVE_PIPELINE_SCRIPTS:-}" ] \
    && [ -d "${COMPARATIVE_PIPELINE_SCRIPTS}" ]; then
    SCRIPTS_DIR="${COMPARATIVE_PIPELINE_SCRIPTS}"
elif [ -n "$(echo "${SELF_DIR}"/master_script_*.py)" ]; then
    SCRIPTS_DIR="${SELF_DIR}"
elif [ -n "$(echo "${SELF_DIR}"/scripts/master_script_*.py)" ]; then
    SCRIPTS_DIR="${SELF_DIR}/scripts"
else
    echo "[ERROR] Could not find the pipeline scripts."
    echo "        Searched: ${SELF_DIR}"
    echo "                  ${SELF_DIR}/scripts"
    exit 1
fi

export COMPARATIVE_PIPELINE_SCRIPTS="${SCRIPTS_DIR}"

masters=("${SCRIPTS_DIR}"/master_script_*.py)

if [ "${#masters[@]}" -ne 1 ]; then
    echo "[ERROR] Expected exactly one master_script_*.py"
    echo "        Directory: ${SCRIPTS_DIR}"
    echo "        Found: ${#masters[@]}"
    if [ "${#masters[@]}" -gt 0 ]; then
        printf '        %s\n' "${masters[@]}"
    fi
    exit 1
fi

MASTER="${masters[0]}"
SETUP_DB="${SCRIPTS_DIR}/00_setup_motif_databases.py"

echo "[INFO] Pipeline scripts  : ${SCRIPTS_DIR}"
echo "[INFO] Master script     : $(basename "${MASTER}")"

echo
echo "[INFO] Verifying required tools..."

if command -v comparative-gene-promoter-check >/dev/null 2>&1; then
    if ! check_output="$(comparative-gene-promoter-check 2>&1)"; then
        echo "${check_output}"
        echo "[ERROR] The environment is incomplete."
        exit 1
    fi
    echo "[OK] Environment check passed."
else
    required_tools=(
        datasets mafft meme tomtom ceqlogo
        orthofinder diamond mcl fasttree
    )

    for tool in "${required_tools[@]}"; do
        if command -v "${tool}" >/dev/null 2>&1; then
            echo "[OK] ${tool}: $(command -v "${tool}")"
        else
            echo "[ERROR] '${tool}' is not on PATH."
            exit 1
        fi
    done
fi

echo
echo "=============================================================================="
echo "                   MOTIF REFERENCE DATABASES"
echo "                         JASPAR + UniProbe"
echo "=============================================================================="

if [ -f "${SETUP_DB}" ]; then
    python -u "${SETUP_DB}"
else
    echo "[ERROR] Motif database setup script was not found:"
    echo "        ${SETUP_DB}"
    exit 1
fi

echo
echo "=============================================================================="
echo "                   STARTING PIPELINE"
echo "                   STAGE 01 IS INTERACTIVE"
echo "=============================================================================="
echo

# Do not use conda run here: Stage 01 requires inherited terminal stdin.
exec python -u "${MASTER}" "$@"
