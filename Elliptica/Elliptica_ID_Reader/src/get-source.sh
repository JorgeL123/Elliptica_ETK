#!/bin/bash
# get-source.sh
# Descarga (si hace falta) y compila la librería externa
# Elliptica_ID_Reader (https://github.com/rashti-alireza/Elliptica_ID_Reader)
# Se invoca desde src/make.code.defn en cada configuración/compilación.

set -e

THORN_SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_DIR="${THORN_SRC_DIR}/Elliptica_ID_Reader_src"
REPO_URL="https://github.com/rashti-alireza/Elliptica_ID_Reader.git"

if [ ! -d "${LIB_DIR}" ]; then
  echo "Elliptica_ID_Reader: clonando libreria externa en ${LIB_DIR}"
  git clone --depth=1 "${REPO_URL}" "${LIB_DIR}"
fi

echo "Elliptica_ID_Reader: compilando libreria externa"
( cd "${LIB_DIR}" && make -f GNUmakefile )