#!/bin/bash
# ==========================================================
# rodar_cpc.sh — rotina diaria completa do CPC
#
# Portado em 2026-11-16 de
#   <usuario>@<servidor>:~/PRODUCAO/CPC/OPERA/rodar_cpc.sh
# que roda no cron do servidor todo dia as 09:00:
#   <agendamento> .../OPERA/rodar_cpc.sh >> .../OPERA/cron.log 2>&1
#
# Este arquivo nao existia em nenhuma das copias locais — era a
# peca que amarrava o sistema inteiro. Copia verbatim do original
# em _legado/servidor/rodar_cpc.sh.servidor
#
# Diferencas em relacao ao original: os scripts agora moram em
# subpastas (OPERACIONAL/, SAE/, SIN/, NFE/, scripts/) e sao
# chamados por caminho absoluto, mas o diretorio de trabalho
# continua plano (OPERA/), porque as rotinas SIMOP leem e gravam
# as planilhas acumuladoras por nome simples no diretorio corrente.
#
# Janela: ultimos 90 dias ate hoje.
# ==========================================================

set -euo pipefail

export TZ=America/Sao_Paulo

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="$RAIZ/OPERA"
LOG="$WORKDIR/log_cpc.txt"
LOCK="/tmp/cpc_cron.lock"

PY="${PY:-python3}"

# Diretorios de shapes das rotinas SIN e NFE (sobrepoem o cpc.config).
# Padroes do servidor de producao; ajuste se rodar em outra maquina.
SHAPEDIR_SIN="${SHAPEDIR_SIN:-/home/<usuario>/STUFF/SHAPES/SIN}"
SHAPEDIR_NFE="${SHAPEDIR_NFE:-/home/<usuario>/STUFF/SHAPES/NFE}"

HOJE=$(date +%d/%m/%Y)
INICIO=$(date -d "90 days ago" +%d/%m/%Y)

mkdir -p "$WORKDIR"

# Evita duas execucoes simultaneas
exec 9>"$LOCK"
flock -n 9 || { echo "[$(date '+%F %T')] Execucao ja em andamento; saindo." >> "$LOG"; exit 0; }

registra() { echo "[$(date '+%F %T')] $*" >> "$LOG"; }

registra "Inicio: janela $INICIO -> $HOJE"
cd "$WORKDIR"

# --- 0) atualiza o repositorio de dados CPC ---
$PY "$RAIZ/OPERACIONAL/cpc_download.py" -i "$INICIO" -f "$HOJE"

# --- 1) SAE (bacias padrao do cpc.config) ---
$PY "$RAIZ/OPERACIONAL/cpc_calcbacia.py" "$INICIO" "$HOJE" --saida chuva_diaria_CPC.xlsx
$PY "$RAIZ/SAE/cpc_SIMOP_SAE.py"

# --- 2) SIN ---
$PY "$RAIZ/OPERACIONAL/cpc_calcbacia.py" "$INICIO" "$HOJE" \
    --saida chuva_SIN.xlsx --shapefile "$SHAPEDIR_SIN" --padrao "SIN*.shp"
$PY "$RAIZ/SIN/cpc_SIMOP_SIN.py"
$PY "$RAIZ/SIN/cpc_SIMOP_SIN_GRAFICOS.py"

# --- 3) NFE ---
$PY "$RAIZ/OPERACIONAL/cpc_calcbacia.py" "$INICIO" "$HOJE" \
    --saida chuva_NFE.xlsx --shapefile "$SHAPEDIR_NFE" --padrao "NFE*.shp"
$PY "$RAIZ/NFE/cpc_SIMOP_NFE.py"
$PY "$RAIZ/NFE/cpc_SIMOP_NFE_GRAFICOS.py"

# --- 4) acumulado mensal SAE + mapas ---
bash "$RAIZ/scripts/roda_SAE_opera.sh"

registra "Fim OK: janela $INICIO -> $HOJE"
