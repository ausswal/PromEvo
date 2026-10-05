#!/usr/bin/env bash

set -euo pipefail

echo "=============================================================================="
echo "             BUILDING COMPARATIVE GENE & PROMOTER PIPELINE"
echo "=============================================================================="

echo "[INFO] Conda prefix : ${PREFIX}"
echo "[INFO] Source       : ${SRC_DIR}"

PACKAGE_DIR="${PREFIX}/share/comparative-gene-promoter"
SCRIPTS_DIR="${PACKAGE_DIR}/scripts"

mkdir -p "${PREFIX}/bin"
mkdir -p "${PACKAGE_DIR}"
mkdir -p "${SCRIPTS_DIR}"

# =============================================================================
# INSTALL PYTHON PIPELINE SCRIPTS
# =============================================================================

echo
echo "[INFO] Installing Python pipeline scripts..."

cp "${SRC_DIR}"/scripts/*.py "${SCRIPTS_DIR}/"

echo "[OK] Python scripts installed."
echo "     ${SCRIPTS_DIR}"

# =============================================================================
# INSTALL PIPELINE LAUNCHER
# =============================================================================

if [ -f "${SRC_DIR}/scripts/run_pipeline.sh" ]; then

    cp "${SRC_DIR}/scripts/run_pipeline.sh" \
       "${PACKAGE_DIR}/run_pipeline.sh"

    chmod +x "${PACKAGE_DIR}/run_pipeline.sh"

    echo "[OK] Pipeline launcher installed."

else

    echo "[ERROR] scripts/run_pipeline.sh was not found."
    exit 1

fi

# =============================================================================
# INSTALL ENVIRONMENT CHECKER
# =============================================================================

if [ -f "${SRC_DIR}/scripts/environment_check.py" ]; then

    cp "${SRC_DIR}/scripts/environment_check.py" \
       "${PACKAGE_DIR}/environment_check.py"

    echo "[OK] Environment checker installed."

else

    echo "[WARN] environment_check.py was not found."

fi

# =============================================================================
# INSTALL DOCUMENTATION
# =============================================================================

if [ -f "${SRC_DIR}/README.md" ]; then

    cp "${SRC_DIR}/README.md" \
       "${PACKAGE_DIR}/README.md"

    echo "[OK] README installed."

fi

if [ -f "${SRC_DIR}/LICENSE.txt" ]; then

    cp "${SRC_DIR}/LICENSE.txt" \
       "${PACKAGE_DIR}/LICENSE.txt"

    echo "[OK] LICENSE installed."

fi

if [ -f "${SRC_DIR}/environment.yml" ]; then

    cp "${SRC_DIR}/environment.yml" \
       "${PACKAGE_DIR}/environment.yml"

    echo "[OK] environment.yml installed."

fi

# =============================================================================
# IDENTIFY MASTER SCRIPT
# =============================================================================

masters=("${SCRIPTS_DIR}"/master_script_*.py)

if [ "${#masters[@]}" -ne 1 ]; then

    echo
    echo "[ERROR] Expected exactly one master_script_*.py."
    echo "        Found ${#masters[@]}."

    if [ "${#masters[@]}" -gt 0 ]; then
        printf '        %s\n' "${masters[@]}"
    fi

    exit 1

fi

MASTER_SCRIPT="${masters[0]}"
MASTER_BASENAME="$(basename "${MASTER_SCRIPT}")"

echo "[OK] Master script detected:"
echo "     ${MASTER_BASENAME}"

# =============================================================================
# MAIN PIPELINE COMMAND
# =============================================================================

cat > "${PREFIX}/bin/comparative-gene-promoter" <<'EOF'
#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="${CONDA_PREFIX}/share/comparative-gene-promoter"

if [ ! -d "${PIPELINE_DIR}" ]; then

    echo "[ERROR] Comparative Gene Promoter installation directory not found:"
    echo "        ${PIPELINE_DIR}"

    exit 1

fi

if [ ! -f "${PIPELINE_DIR}/run_pipeline.sh" ]; then

    echo "[ERROR] Pipeline launcher not found:"
    echo "        ${PIPELINE_DIR}/run_pipeline.sh"

    exit 1

fi

exec bash "${PIPELINE_DIR}/run_pipeline.sh" "$@"
EOF

chmod +x "${PREFIX}/bin/comparative-gene-promoter"

# =============================================================================
# MASTER SCRIPT COMMAND
# =============================================================================

cat > "${PREFIX}/bin/comparative-gene-promoter-master" <<'EOF'
#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="${CONDA_PREFIX}/share/comparative-gene-promoter"
MASTER="${PIPELINE_DIR}/scripts/master_script_10.py"

if [ ! -f "${MASTER}" ]; then

    echo "[ERROR] Master script not found:"
    echo "        ${MASTER}"

    exit 1

fi

exec python "${MASTER}" "$@"
EOF

chmod +x "${PREFIX}/bin/comparative-gene-promoter-master"

# =============================================================================
# ENVIRONMENT CHECK COMMAND
# =============================================================================

cat > "${PREFIX}/bin/comparative-gene-promoter-check" <<'EOF'
#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="${CONDA_PREFIX}/share/comparative-gene-promoter"
CHECKER="${PIPELINE_DIR}/environment_check.py"

if [ ! -f "${CHECKER}" ]; then

    echo "[ERROR] Environment checker not found:"
    echo "        ${CHECKER}"

    exit 1

fi

exec python "${CHECKER}" "$@"
EOF

chmod +x "${PREFIX}/bin/comparative-gene-promoter-check"

# =============================================================================
# FINAL BUILD INFORMATION
# =============================================================================

echo
echo "=============================================================================="
echo "                     BUILD INSTALLATION COMPLETE"
echo "=============================================================================="

echo "[INFO] Package directory : ${PACKAGE_DIR}"
echo "[INFO] Scripts directory : ${SCRIPTS_DIR}"
echo "[INFO] Master script      : ${MASTER_BASENAME}"

echo
echo "[INFO] Installed commands:"
echo "       comparative-gene-promoter"
echo "       comparative-gene-promoter-master"
echo "       comparative-gene-promoter-check"

echo
echo "[INFO] Motif databases are NOT packaged."
echo "[INFO] JASPAR and UniProbe databases will be downloaded by Stage 00."
echo
echo "=============================================================================="