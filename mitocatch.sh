#!/bin/bash

MODO="getorganelle"
ORGANELLE_TYPE=""
CPUS=6

while [[ $# -gt 0 ]]; do
  case $1 in
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --modo) MODO="$2"; shift 2 ;;
    --organelle_type) ORGANELLE_TYPE="$2"; shift 2 ;;
    *) shift ;;
  esac
done

if [ -z "$R1" ] || [ -z "$R2" ]; then
    echo "Erro: Deves fornecer as duas reads (--r1 e --r2)."
    exit 1
fi

R1_BASE=$(basename "$R1")
R2_BASE=$(basename "$R2")

if [ "$MODO" == "getorganelle" ]; then
    ASSEMBLERS="getorganelle"
elif [ "$MODO" == "trinity" ]; then
    ASSEMBLERS="trinity"
else
    echo "Erro: Modo desconhecido. Escolha 'getorganelle' ou 'trinity'."
    exit 1
fi

WORKDIR="run_${MODO}"
mkdir -p "$WORKDIR"

mkdir -p "$WORKDIR/data/raw"
ln -snf "../workflow" "$WORKDIR/workflow"

ln -snf "$(realpath "$R1")" "$WORKDIR/data/raw/$R1_BASE"
ln -snf "$(realpath "$R2")" "$WORKDIR/data/raw/$R2_BASE"

snakemake --nolock --use-conda --snakefile workflow/Snakefile \
          --cores $CPUS \
          --directory "$WORKDIR" \
          --config r1="data/raw/$R1_BASE" r2="data/raw/$R2_BASE" assemblers="$ASSEMBLERS" outdir="results" organelle_type="$ORGANELLE_TYPE"

