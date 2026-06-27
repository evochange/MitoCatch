#!/bin/bash

# 1. Valores padrão
MODO="getorganelle"
ORGANELLE_TYPE=""
CPUS=6

# 2. Ler a linha de comando
while [[ $# -gt 0 ]]; do
  case $1 in
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --modo) MODO="$2"; shift 2 ;;
    --organelle_type) ORGANELLE_TYPE="$2"; shift 2 ;;
    *) shift ;;
  esac
done

# Verificar se os inputs obrigatórios existem
if [ -z "$R1" ] || [ -z "$R2" ]; then
    echo "Erro: Deves fornecer as duas reads (--r1 e --r2)."
    exit 1
fi

R1_BASE=$(basename "$R1")
R2_BASE=$(basename "$R2")

# 3. Traduzir o modo para os assemblers do Snakemake
if [ "$MODO" == "getorganelle" ]; then
    ASSEMBLERS="getorganelle"
elif [ "$MODO" == "trinity" ]; then
    ASSEMBLERS="trinity"
else
    echo "Erro: Modo desconhecido. Escolha 'getorganelle' ou 'trinity'."
    exit 1
fi

# 4. AUTOMAÇÃO: Definir uma pasta de trabalho única para esta corrida
WORKDIR="run_${MODO}"
mkdir -p "$WORKDIR"

# Criar a estrutura interna necessária para as reads serem encontradas pelo Snakefile
mkdir -p "$WORKDIR/data/raw"
ln -snf "../workflow" "$WORKDIR/workflow"

# Criar links simbólicos para manter o isolamento
ln -snf "$(realpath "$R1")" "$WORKDIR/data/raw/$R1_BASE"
ln -snf "$(realpath "$R2")" "$WORKDIR/data/raw/$R2_BASE"


# 5. Chamar o Snakemake apontando o outdir de forma relativa interna
snakemake --nolock --use-conda --snakefile workflow/Snakefile \
          --cores $CPUS \
          --directory "$WORKDIR" \
          --config r1="data/raw/$R1_BASE" r2="data/raw/$R2_BASE" assemblers="$ASSEMBLERS" outdir="results" organelle_type="$ORGANELLE_TYPE"

