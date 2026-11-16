#!/bin/bash
# ==========================================================
# coleta_servidor.sh
#
# Coleta o sistema CPC do servidor operacional para comparacao
# com esta copia de desenvolvimento.
#
# RODE ISTO NO WSL (Ubuntu) da sua maquina, nao no servidor:
#
#     bash <raiz-do-repositorio>/scripts/coleta_servidor.sh
#
# Gera em _backup/:
#     INVENTARIO_SERVIDOR.txt   arvore, hashes, crontab, pacotes
#     CPC_SERVIDOR.tar.gz       codigo e config do servidor
# ==========================================================

set -uo pipefail

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# ----------------------------------------------------------
# Identidade do ambiente — vem de cpc.local, que fica fora do git
# ----------------------------------------------------------
# Este repositorio e publico e nao carrega onde a producao roda:
# servidor, usuario e caminhos ficam em cpc.local, barrado pelo
# .gitignore. Modelo em cpc.local.example.
CPC_LOCAL="${CPC_LOCAL:-$RAIZ/cpc.local}"
if [ -f "$CPC_LOCAL" ]; then
    . "$CPC_LOCAL"
fi

if [ -z "${SERVIDOR:-}" ]; then
    echo "ERRO: SERVIDOR nao definido." >&2
    echo "      Crie o cpc.local a partir do modelo:" >&2
    echo "        cp '$RAIZ/cpc.local.example' '$RAIZ/cpc.local'" >&2
    echo "      e preencha. Ou passe pelo ambiente:" >&2
    echo "        SERVIDOR=usuario@host bash \$0" >&2
    exit 1
fi

DEST="${DEST:-$RAIZ/_backup}"

mkdir -p "$DEST"
INV="$DEST/INVENTARIO_SERVIDOR.txt"
TGZ="$DEST/CPC_SERVIDOR.tar.gz"

echo ">> Coletando inventario de ${SERVIDOR} ..."

# ----- etapa 1: descobrir onde mora o CPC e inventariar -----
ssh -o ConnectTimeout=15 "$SERVIDOR" 'bash -s' > "$INV" <<'REMOTO'
echo "host=$(hostname)  usuario=$(whoami)  data=$(date -Is)"
echo

# acha a raiz do CPC: o diretorio que contem cpc.config ou cpc_download.py
CAND=$(find "$HOME" /opt /srv /usr/local -maxdepth 5 -type f \
         \( -iname "cpc.config" -o -iname "cpc_download.py" \) 2>/dev/null \
       | xargs -r -n1 dirname | sort -u)

echo "### candidatos a raiz do CPC"
echo "${CAND:-(nenhum encontrado)}"
echo

RAIZ=$(echo "$CAND" | head -1)
echo "CPC_DIR=${RAIZ}"
echo

if [ -n "$RAIZ" ] && [ -d "$RAIZ" ]; then
  echo "### arvore de ${RAIZ}"
  find "$RAIZ" -not -path "*/.git/*" -not -name "*.nc" \
       -printf "%y %10s  %TY-%Tm-%Td %TH:%TM  %p\n" 2>/dev/null | sort -k5
  echo
  echo "### hashes (md5 | bytes | mtime | caminho) — codigo e config"
  find "$RAIZ" -type f \( -name "*.py" -o -name "*.sh" -o -name "*.config" -o -name "*.cfg" -o -name "*.txt" -o -name "*.md" \) \
       -not -path "*/.git/*" 2>/dev/null \
    | while IFS= read -r f; do
        printf "%s | %s | %s | %s\n" \
          "$(md5sum "$f" | cut -d' ' -f1)" \
          "$(stat -c%s "$f")" \
          "$(stat -c%y "$f" | cut -d. -f1)" \
          "$f"
      done
  echo
  echo "### git (se houver)"
  if [ -d "$RAIZ/.git" ]; then
    git -C "$RAIZ" log --oneline -10 2>&1
    echo "-- status --"; git -C "$RAIZ" status --short 2>&1 | head -40
    echo "-- remote --"; git -C "$RAIZ" remote -v 2>&1
  else
    echo "(nao e repositorio git)"
  fi
  echo
  echo "### conteudo do cpc.config"
  cat "$RAIZ/cpc.config" 2>/dev/null || echo "(sem cpc.config)"
fi

echo
echo "### crontab do usuario"
crontab -l 2>/dev/null || echo "(sem crontab)"
echo
echo "### entradas de cron do sistema mencionando cpc"
grep -ril "cpc" /etc/cron* 2>/dev/null | head -20
echo
echo "### servicos systemd mencionando cpc"
systemctl list-units --all 2>/dev/null | grep -i cpc || echo "(nenhum)"
echo
echo "### python e pacotes relevantes"
python3 -V 2>&1
pip3 list 2>/dev/null | grep -iE "cartopy|geopandas|netCDF4|shapely|openpyxl|scipy|pandas|numpy|matplotlib" || echo "(pip3 nao listou)"
REMOTO

rc=$?
if [ $rc -ne 0 ]; then
    echo "!! Falha ao conectar em ${SERVIDOR} (codigo $rc)." >&2
    echo "   Verifique a chave SSH ou rode: ssh ${SERVIDOR}" >&2
    exit $rc
fi

echo ">> Inventario salvo em: $INV"

# ----- etapa 2: empacotar o codigo do servidor -----
RAIZ=$(grep -m1 '^CPC_DIR=' "$INV" | cut -d= -f2-)

if [ -z "$RAIZ" ]; then
    echo "!! Nao localizei a raiz do CPC no servidor." >&2
    echo "   Abra $INV, veja a secao 'candidatos', e rode de novo assim:" >&2
    echo "   CPC_DIR=/caminho/correto bash \$0" >&2
    exit 1
fi

echo ">> Empacotando ${RAIZ} ..."
ssh -o ConnectTimeout=15 "$SERVIDOR" \
    "tar -czf - --exclude='*.nc' --exclude='*.png' --exclude='*.xlsx' \
         --exclude='.git' --exclude='__pycache__' \
         -C \"\$(dirname '$RAIZ')\" \"\$(basename '$RAIZ')\"" > "$TGZ"

if [ -s "$TGZ" ]; then
    echo ">> Pacote salvo em: $TGZ  ($(du -h "$TGZ" | cut -f1))"
    echo
    echo "Pronto. Volte na conversa e diga que a coleta terminou."
else
    echo "!! O pacote saiu vazio — confira o inventario." >&2
    exit 1
fi
