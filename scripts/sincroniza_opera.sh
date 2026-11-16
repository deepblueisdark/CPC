#!/bin/bash
# ==========================================================
# sincroniza_opera.sh
#
# Traz o OPERA/ do servidor de producao para dentro deste
# repositorio, em _producao/OPERA/, para que o que roda de fato
# fique sob historico de versao.
#
# E UMA SINCRONIZACAO DE CAPTURA, EM UM SENTIDO SO:
#   servidor  -->  repositorio
# O script NUNCA escreve no servidor. Todo acesso remoto e
# leitura (ssh + tar para stdout).
#
# RODE ISTO NO WSL (Ubuntu) da sua maquina, nao no servidor:
#
#     bash <raiz-do-repositorio>/scripts/sincroniza_opera.sh
#
# Gera/atualiza:
#     _producao/OPERA/              codigo e config do OPERA/, verbatim
#     _producao/MANIFESTO.txt       md5 do que foi capturado
#     _producao/ARVORE_OPERA.txt    arvore completa do OPERA/ no servidor,
#                                   inclusive os dados que nao sao copiados
#
# Variaveis que dao para sobrepor:
#     SERVIDOR=<usuario>@<servidor>
#     OPERA_REMOTO=/home/<usuario>/PRODUCAO/CPC/OPERA
# ==========================================================

set -uo pipefail

SERVIDOR="${SERVIDOR:-<usuario>@<servidor>}"
OPERA_REMOTO="${OPERA_REMOTO:-/home/<usuario>/PRODUCAO/CPC/OPERA}"

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$RAIZ/_producao"
SNAP="$DEST/OPERA"
MANIFESTO="$DEST/MANIFESTO.txt"
ARVORE="$DEST/ARVORE_OPERA.txt"

# Extensoes de codigo e config: o que entra no snapshot.
# Planilhas, PNG, NetCDF e logs sao dado operacional e ficam fora —
# note que log_cpc.txt casaria com *.txt, e por isso excluido por nome
# junto com *.log na etapa 2.
INCLUI=( "*.py" "*.sh" "*.config" "*.cfg" "*.txt" "*.md" "*.json" "*.yaml" "*.yml" )

if [ ! -f "$RAIZ/cpc.config" ]; then
    echo "!! Nao parece ser a raiz do repositorio CPC: $RAIZ" >&2
    exit 1
fi

echo ">> Servidor : $SERVIDOR"
echo ">> Remoto   : $OPERA_REMOTO"
echo ">> Destino  : $SNAP"
echo

# ----------------------------------------------------------
# etapa 1 — arvore completa do OPERA/, para registro
# ----------------------------------------------------------
echo ">> Lendo a arvore do OPERA/ ..."
ssh -o ConnectTimeout=15 "$SERVIDOR" "bash -s" > "$ARVORE.novo" <<REMOTO
set -u
OP="$OPERA_REMOTO"
echo "# arvore de \$OP"
echo "# host=\$(hostname)  usuario=\$(whoami)  data=\$(date -Is)"
echo
if [ ! -d "\$OP" ]; then
    echo "!! \$OP nao existe no servidor"
    exit 2
fi
find "\$OP" -not -path "*/__pycache__/*" \
     -printf "%y %10s  %TY-%Tm-%Td %TH:%TM  %P\n" 2>/dev/null | sort -k5
echo
echo "# --- hashes do codigo e dos configs ---"
find "\$OP" -maxdepth 2 -type f \
     \( -name "*.py" -o -name "*.sh" -o -name "*.config" -o -name "*.cfg" \
        -o -name "*.txt" -o -name "*.md" -o -name "*.json" -o -name "*.yaml" -o -name "*.yml" \) \
     -not -path "*/__pycache__/*" 2>/dev/null \
  | sort \
  | while IFS= read -r f; do
      printf "%s  %s\n" "\$(md5sum "\$f" | cut -d" " -f1)" "\${f#\$OP/}"
    done
REMOTO

rc=$?
if [ $rc -ne 0 ]; then
    echo "!! Falha ao ler o servidor (codigo $rc)." >&2
    echo "   Teste a conexao com: ssh $SERVIDOR" >&2
    rm -f "$ARVORE.novo"
    exit $rc
fi
mkdir -p "$DEST"
mv "$ARVORE.novo" "$ARVORE"
echo "   arvore salva em: $ARVORE"
echo

# ----------------------------------------------------------
# etapa 2 — trazer o codigo
# ----------------------------------------------------------
echo ">> Empacotando o codigo do OPERA/ ..."

ARGS_TAR=""
for padrao in "${INCLUI[@]}"; do
    ARGS_TAR="$ARGS_TAR -name '$padrao' -o"
done
ARGS_TAR="${ARGS_TAR% -o}"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ssh -o ConnectTimeout=15 "$SERVIDOR" \
    "cd '$OPERA_REMOTO' && find . -maxdepth 2 -type f \\( $ARGS_TAR \\) \
     -not -path './__pycache__/*' \
     -not -name 'log_cpc.txt' -not -name '*.log' -print0 | tar -czf - --null -T -" \
  | tar -xzf - -C "$TMP"

if [ ! "$(ls -A "$TMP" 2>/dev/null)" ]; then
    echo "!! O pacote saiu vazio. Confira $ARVORE." >&2
    exit 1
fi

# Troca o snapshot inteiro, para que arquivo apagado no servidor
# tambem desapareca daqui. O snapshot e derivado: da para refazer.
rm -rf "$SNAP"
mkdir -p "$SNAP"
cp -a "$TMP"/. "$SNAP"/

N=$(find "$SNAP" -type f | wc -l)
echo "   $N arquivos capturados em: $SNAP"
echo

# ----------------------------------------------------------
# etapa 3 — manifesto
# ----------------------------------------------------------
{
    echo "# MANIFESTO — snapshot de $OPERA_REMOTO"
    echo "# origem : $SERVIDOR"
    echo "# data   : $(date -Is)"
    echo "# gerado por scripts/sincroniza_opera.sh"
    echo
    (
      cd "$SNAP" || exit 1
      find . -type f | sed 's|^\./||' | sort \
        | while IFS= read -r f; do
            printf "%s  %s\n" "$(md5sum "$f" | cut -d' ' -f1)" "$f"
          done
    )
} > "$MANIFESTO"
echo ">> Manifesto em: $MANIFESTO"
echo

# ----------------------------------------------------------
# etapa 4 — o que mudou desde o snapshot anterior
# ----------------------------------------------------------
if git -C "$RAIZ" rev-parse --git-dir >/dev/null 2>&1; then
    echo ">> Mudancas na producao desde o ultimo snapshot:"
    if git -C "$RAIZ" diff --quiet -- _producao/OPERA 2>/dev/null \
       && [ -z "$(git -C "$RAIZ" ls-files --others --exclude-standard -- _producao/OPERA)" ]; then
        echo "   (nenhuma — a producao esta igual ao ultimo snapshot commitado)"
    else
        git -C "$RAIZ" status --short -- _producao/OPERA
        echo
        echo "   Detalhe: git -C \"$RAIZ\" diff -- _producao/OPERA"
    fi
    echo
fi

# ----------------------------------------------------------
# etapa 5 — snapshot x codigo ativo do repositorio
# ----------------------------------------------------------
# Compara por nome de arquivo, ignorando CRLF e espaco no fim da
# linha. Diferencas ESPERADAS estao em _producao/DIVERGENCIAS.md
# (bloco _raiz_projeto(), assinaturas com CONFIG_FILE, layout em
# subpastas e as latitudes [SAE]). O que aparecer fora dessa lista
# e drift de verdade e merece olhada.
echo ">> Snapshot x codigo ativo do repositorio (normalizando CRLF):"
norm() { sed -e 's/\r$//' -e 's/[[:space:]]*$//' "$1" | md5sum | cut -d' ' -f1; }

while IFS= read -r f; do
    base="$(basename "$f")"
    ativo="$(find "$RAIZ" -path "$RAIZ/_producao" -prune -o \
                          -path "$RAIZ/_legado" -prune -o \
                          -path "$RAIZ/.git" -prune -o \
                          -type f -name "$base" -print 2>/dev/null | head -1)"
    if [ -z "$ativo" ]; then
        printf "   %-34s  SO NA PRODUCAO\n" "$base"
    elif [ "$(norm "$f")" = "$(norm "$ativo")" ]; then
        printf "   %-34s  identico\n" "$base"
    else
        printf "   %-34s  difere de %s\n" "$base" "${ativo#$RAIZ/}"
    fi
done < <(find "$SNAP" -type f | sort)

# ----------------------------------------------------------
# etapa 6 — cpc.config: comparacao chave por chave
# ----------------------------------------------------------
# No cpc.config o diff de texto engana: as duas copias tem secoes
# em ordem diferente, e o espacamento entre os valores das escalas
# conta como diferenca sem ser uma. Aqui a comparacao e semantica.
# Quais diferencas sao esperadas: _producao/DIVERGENCIAS.md
if [ -f "$SNAP/cpc.config" ] && command -v python3 >/dev/null 2>&1; then
    echo ">> cpc.config — producao x repositorio, chave por chave:"
    python3 - "$SNAP/cpc.config" "$RAIZ/cpc.config" <<'FIMPY'
import configparser, sys

def carrega(caminho):
    c = configparser.ConfigParser()
    c.read(caminho, encoding="utf-8")
    d = {}
    for k, v in c.defaults().items():
        d[("DEFAULT", k)] = v.strip()
    padrao = set(c.defaults())
    for s in c.sections():
        for k in c.options(s):
            if k not in padrao:
                d[(s, k)] = c.get(s, k).strip()
    return d

def corta(v):
    if v is None:
        return "(ausente)"
    return v if len(v) <= 46 else v[:43] + "..."

prod, repo = carrega(sys.argv[1]), carrega(sys.argv[2])
n = 0
for chave in sorted(prod.keys() | repo.keys()):
    p, r = prod.get(chave), repo.get(chave)
    if p == r:
        continue
    n += 1
    secao, k = chave
    if p is not None and r is not None and p.split() == r.split():
        print("   [%s] %s: so espacamento, valores iguais" % (secao, k))
        continue
    print("   [%s] %s" % (secao, k))
    print("       producao    : %s" % corta(p))
    print("       repositorio : %s" % corta(r))
print("   (nenhuma diferenca)" if n == 0 else
      "   %d chaves diferentes — confira contra _producao/DIVERGENCIAS.md" % n)
FIMPY
    echo
fi

echo
echo "Pronto. Volte na conversa e diga que a sincronizacao terminou."
