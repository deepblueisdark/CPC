#!/bin/bash
# ==========================================================
# roda_SAE_opera.sh
#
# Gera o acumulado de chuva do mes corrente (dia 1 ate hoje-2)
# para as bacias SAE e plota os mapas.
#
# Atualizado em 2026-11-16 na consolidacao do sistema:
#   - os scripts mudaram para PLOTAGEM/, os caminhos agora sao
#     resolvidos a partir da raiz do projeto (nao mais do CWD);
#   - as saidas vao para SAIDAS/CLIMA/ em vez da pasta do codigo;
#   - a regiao [SAE] e a chave escala_cores_mensal passaram a
#     existir no cpc.config (antes o script abortava com
#     ValueError nos dois casos).
#
# Requer GNU date (Linux). Em macOS, instale coreutils e use gdate.
# ==========================================================

set -euo pipefail

# ----- raiz do projeto (pasta acima de scripts/) -----
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ----- shapefile das bacias SAE -----
# Ajuste SHAPEFILE se o arquivo estiver em outro lugar. O valor
# padrao segue o shapefile_dir do cpc.config.
# ATENCAO: este .shp NAO esta dentro do repositorio.
SHAPEFILE="${SHAPEFILE:-${RAIZ}/SHAPES/SAE/SAE_TODAS_AS_BACIAS_HICON_POLY.shp}"

SAIDA_DIR="${RAIZ}/SAIDAS/CLIMA"
NETCDF="${SAIDA_DIR}/SAE_opera.nc"

if [[ ! -f "${SHAPEFILE}" ]]; then
    echo "[ERRO] Shapefile nao encontrado: ${SHAPEFILE}" >&2
    echo "       Defina a variavel SHAPEFILE antes de rodar, ex:" >&2
    echo "       SHAPEFILE=/caminho/SAE_TODAS_AS_BACIAS_HICON_POLY.shp $0" >&2
    exit 1
fi

mkdir -p "${SAIDA_DIR}"

# ----- monta as datas -----
# O CPC so tem dados consolidados ate "hoje - 2 dias". A data-base
# com defasagem de 2 dias tambem define o mes, para funcionar quando
# a virada cai no comeco do mes (ex: hoje e dia 1 ou 2).
DIA_FIM=$(date -d "-2 days" +%-d)
MES=$(date -d "-2 days" +%-m)
ANO=$(date -d "-2 days" +%Y)

DATA_INICIO="1/${MES}/${ANO}"
DATA_FIM="${DIA_FIM}/${MES}/${ANO}"

echo "Periodo: ${DATA_INICIO} ate ${DATA_FIM}"
echo "Saidas em: ${SAIDA_DIR}"

# ----- 1. gera o NetCDF acumulado -----
python "${RAIZ}/PLOTAGEM/cpc_gera_plota_chuva_acumulada.py" "${DATA_INICIO}" "${DATA_FIM}" \
    --saida "${NETCDF}" \
    --resolucao 0.05 \
    --interpolacao \
    --shapefile "${SHAPEFILE}" \
    --figura "${SAIDA_DIR}/SAE_OPERA" \
    --titulo "CHUVA TOTAL " \
    --regiao SAE \
    --cores escala_cores_mensal \
    --politica off \
    --cor-shapefile "black" \
    --unidade mm/mes

echo "NetCDF gerado: ${NETCDF}"

# ----- 2. gera as imagens a partir do NetCDF -----
echo "Gerando imagens..."
python "${RAIZ}/PLOTAGEM/cpc_plota_clima_unificado.py" "${NETCDF}" \
    --shapefile "${SHAPEFILE}" \
    --titulo "CLIMATOLOGIA DE PRECIPITACAO" \
    --regiao SAE \
    --cores escala_cores_mensal \
    --unidade mm/mes \
    --politica off \
    --cor-shapefile "black" \
    --tipo ambos \
    --saida-dir "${SAIDA_DIR}"

echo "Concluido. Resultados em ${SAIDA_DIR}"
