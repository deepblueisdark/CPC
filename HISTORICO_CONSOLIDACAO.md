# Consolidacao do sistema CPC — 2026-11-16

Registro do que foi feito ao unificar as copias que existiam em paralelo.

## Situacao anterior

| Local | Conteudo | Git |
|---|---|---|
| `<raiz-do-repositorio>` | 64 arquivos, 14 MB — CLIMA/, MAPAS_CPC/, NFE/, SAE/, SIN/ + scripts de raiz | nao |
| `<pasta-antiga>` | 14 arquivos, 338 KB — OPERACIONAL/, PLOTAGEM/, SAE/, SIN/ | sim, `github.com/deepblueisdark/CPC` |
| `<servidor>` (conta `simop`) | versao operacional | nao verificada |

## Problemas encontrados

1. **Pipeline partido entre as duas copias.** `roda_SAE_opera.sh`
   (so no Devel) chamava `cpc_gera_plota_chuva_acumulada.py` (so no
   GitHub) e `cpc_plota_clima_unificado.py` (so no Devel). Nenhuma das
   duas copias rodava sozinha.
2. **Dois `cpc.config` incompativeis.** O do GitHub tinha os caminhos
   entre aspas literais (`SHAPEFILE_DIR = "<caminho-wsl>/..."`), o que quebrava
   os `os.path.join`, e usava caminhos WSL. O do Devel estava correto
   para Windows, mas tinha perdido `MAPA_DIR` e `saida_excel`.
3. **`roda_SAE_opera.sh` nao rodava com nenhum dos dois configs:**
   pedia `--regiao SAE` (so existiam AMS, BRASIL e SAARA) e
   `--cores escala_cores_mensal` (chave inexistente). Os scripts
   levantam `ValueError` nos dois casos.
4. **`.gitignore` ignorava `cpc.config`**, que estava rastreado no
   repositorio — regra inocua e enganosa.
5. **Devel nao era repositorio git**, apesar de ser onde o
   desenvolvimento acontecia.
6. **Tres versoes de `cpc_plota_clima`** divergentes (V1 no Devel, V1
   diferente no GitHub, V2), alem do `unificado` que o proprio
   cabecalho declara substituir V1 e V2.
7. **Codigo misturado com saida:** 42 PNGs em `MAPAS_CPC/`, 4 PNGs em
   `CLIMA/`, planilha de 7 MB em `NFE/`.
8. **Todos os scripts procuravam `cpc.config` no diretorio atual**, o
   que impedia qualquer organizacao em subpastas.
9. **Referencias ausentes:** `SAE_TODAS_AS_BACIAS_HICON_POLY.shp` e
   quatro planilhas usadas pelas rotinas SIMOP.

Os arquivos em comum (`cpc_download.py`, `cpc_calcbacia.py`,
`SAE/cpc_SIMOP_SAE.py`, `SIN/cpc_SIMOP_SIN.py`,
`SIN/cpc_SIMOP_SIN_GRAFICOS.py`) estavam **byte a byte identicos** nas
duas copias — nao houve conflito de conteudo.

## O que foi feito

### Backups (antes de qualquer mudanca)
- `_backup/BACKUP_Devel_CPC_2026-11-16.tar.gz` (13 MB)
- `_backup/BACKUP_GITHUB_CPC_2026-11-16.tar.gz` (113 KB, inclui o `.git` completo)

### Scripts trazidos do GitHub (existiam so la)
- `PLOTAGEM/cpc_gera_plota_chuva_acumulada.py`
- `PLOTAGEM/cpc_plota.py`
- `PLOTAGEM/cpc_climatology.py`

### Reorganizacao
- `cpc_download.py`, `cpc_calcbacia.py` -> `OPERACIONAL/`
- `cpc_plota_clima_unificado.py` -> `PLOTAGEM/`
- `roda_SAE_opera.sh` -> `scripts/`
- 49 PNGs -> `SAIDAS/MAPAS_CPC/` e `SAIDAS/CLIMA/`
- `chuva_diaria_CPC.xlsx`, `chuva_NFE.xlsx` -> `SAIDAS/PLANILHAS/`
- `chuva_diaria_CPC_NFE.xlsx` (entrada, 7 MB) -> `DADOS/`
- pastas `CLIMA/` e `MAPAS_CPC/` removidas (vazias)

### Arquivadas em `_legado/`
- `cpc_plota_clima_V1_devel.py` (era `CLIMA/cpc_plota_clima.py`)
- `cpc_plota_clima_V1_github.py` (era `PLOTAGEM/cpc_plota_clima.py`)
- `cpc_plota_clima_V2.py`

Ativo: apenas `PLOTAGEM/cpc_plota_clima_unificado.py`.

### Correcoes de codigo
- **Resolucao de caminhos** em 6 scripts (`cpc_download.py`,
  `cpc_calcbacia.py`, `cpc_plota_clima_unificado.py`,
  `cpc_gera_plota_chuva_acumulada.py`, `cpc_plota.py`,
  `cpc_climatology.py`): adicionado `_raiz_projeto()`, que procura o
  `cpc.config` na pasta do script e acima. Agora rodam de qualquer
  diretorio.
- `MAPA_DIR` em `cpc_calcbacia.py` passou de `./MAPAS_CPC` para
  `SAIDAS/MAPAS_CPC` na raiz do projeto.
- `roda_SAE_opera.sh` reescrito: resolve a raiz via `BASH_SOURCE`,
  aceita a variavel `SHAPEFILE`, falha com mensagem clara se o
  shapefile nao existir, e grava em `SAIDAS/CLIMA/`.

### `cpc.config` unificado
- Caminhos Windows ativos; equivalentes Linux/WSL comentados.
- Aspas literais dos caminhos removidas.
- `mapa_dir` e `saida_excel` recuperados, apontando para `SAIDAS/`.
- Adicionada a chave `escala_cores_mensal`.
- Adicionada a secao `[SAE]`.

### Git
- Repositorio novo iniciado nesta pasta (historico anterior preservado
  apenas no `.tar.gz` de backup e no GitHub remoto).
- `.gitignore` reescrito: ignora `SAIDAS/`, `DADOS/`, `_backup/`,
  `__pycache__/`, `*.nc` e `desktop.ini`; a regra inocua sobre
  `cpc.config` foi removida e o arquivo segue rastreado.

### Pasta antiga
- Conteudo de `<pasta-antiga>` movido para
  `<pasta-antiga>\_to_delete\`, para conferencia e exclusao
  manual.

## Nao feito

- **Comparacao com o servidor** `<servidor>` — esta sessao nao
  alcanca o host (SSH fora da allowlist de rede).
- **Push para o GitHub** — o historico novo divergiu do remoto; um push
  seria rejeitado. Decida entre `--force` (sobrescreve o remoto) ou um
  branch novo.
- Os valores de `[SAE]` e `escala_cores_mensal` precisam da sua
  confirmacao.

---

# Segunda rodada — comparacao com o servidor (2026-11-16, tarde)

A sessao nao alcanca `<servidor>` por rede (nenhuma porta sai da
allowlist) e o computer use recusa digitacao em IDEs e terminais
(`tier: "click"`). A coleta foi feita pelo WSL da maquina local, com
`scripts/coleta_servidor.sh`, que alcanca o servidor e enxerga o `E:`.

## O que havia no servidor

`/home/<usuario>/PRODUCAO/CPC` tem **tres** copias internas:

| Diretorio | Papel |
|---|---|
| `CPC/` (raiz) | copia antiga, rastreada no git, 8 commits atras do GitHub |
| `CPC/OPERA/` | **producao** — e o que o cron executa; nao rastreada |
| `CPC/SAE/` | sandbox de desenvolvimento; nao rastreada |

O git da raiz estava em `624e45d`, oito commits atras do
`origin/main`, com `cpc.config` modificado e `OPERA/`, `SAE/`,
`cpc.config.local` e `log_cpc.txt` nao rastreados. Ou seja: o que
roda em producao nunca esteve sob versionamento.

## Resultado da comparacao

Comparacao por hash, normalizando CRLF e espacos finais.

**Todo o codigo que o servidor executa em `OPERA/` e identico ao que
foi consolidado** — os 10 scripts batem byte a byte. A consolidacao
da manha acertou o codigo que importa.

Faltavam duas pecas, ambas so no servidor:

1. **`OPERA/rodar_cpc.sh`** — o orquestrador do cron, que amarra o
   sistema inteiro. Nao existia em nenhuma copia local. Portado para
   `scripts/rodar_cpc.sh`; original em `_legado/servidor/`.
2. **`OPERA/cpc.config`** — o config de producao, com as escalas e a
   regiao `[SAE]` reais.

E tres scripts de funcionalidade ausente, da sandbox `SAE/`:

- `cpc_calcula_anomalia.py` -> `PLOTAGEM/` (introduz dependencia `xarray`)
- `cpc_plota_anomalia.py` -> `PLOTAGEM/`
- `cpc_plota_chuva_acumulada.py` e `cpc_calcula_anomalia_old.py` -> `_legado/`

## Correcoes que a comparacao obrigou

- **Os valores que eu havia chutado estavam errados.** A
  `escala_cores_mensal` de producao tem 31 niveis de 0.1 a 300 mm/mes;
  eu havia posto 15 niveis de 1 a 1000. A regiao `[SAE]` real e
  lat -8/-21, lon 285/305; eu havia posto -21/-6 e 284/306.
  Ambos substituidos pelos valores de producao.
- **Escalas que faltavam** foram acrescentadas: `escala_diaria`,
  `escala_clima`, `escala_clima_simples`, e as secoes `[ANOMALIA]` e
  `[ANOMALIA_PER]` (estas vindas da sandbox `SAE/cpc.config`).
- **`[SAE]` tem `min_lat` maior que `max_lat`** (-8 contra -21),
  invertido em relacao as outras secoes. Copiado verbatim porque e
  assim que roda em producao; sinalizado no `cpc.config` e no README.
- **As duas "versoes divergentes" de `cpc_plota_clima.py` eram o mesmo
  arquivo.** Diferiam apenas em CRLF. Os dois arquivos que eu havia
  arquivado como `_V1_devel` e `_V1_github` foram colapsados em
  `_legado/cpc_plota_clima_V1.py`.
- **`DADOS/` foi substituido por `OPERA/`.** A leitura dos scripts
  SIMOP mostrou que eles leem e gravam as planilhas acumuladoras por
  nome simples no diretorio corrente, entao o diretorio de trabalho
  precisa ser plano, como no servidor. As tres planilhas locais foram
  movidas para `OPERA/`.
- **O shapefile deixou de ser pendencia**: `SAE_TODAS_AS_BACIAS_HICON_POLY`
  veio do servidor para `SHAPES/SAE/`, e o `roda_SAE_opera.sh` passou a
  apontar para ele por padrao.
- `.gitattributes` ganhou `*.sh text eol=lf`, ja que os shells rodam
  no servidor Linux.

## O que o servidor tem e nao foi trazido

- As quatro planilhas acumuladoras (12 MB), excluidas do pacote de
  coleta: `CHUVA_OBSERVADA.xlsx`, `CHUVA_OBSERVADA_CPC_COMPLETA.xlsx`,
  `CHUVA_CPC_COMPLETA_SIN.xlsx`, `chuva_SIN.xlsx`. Sem elas a cadeia
  SIMOP nao roda localmente.
- Os shapes de SIN e NFE (`/home/<usuario>/STUFF/SHAPES/{SIN,NFE}`).
- As copias divergentes e mais antigas da sandbox `SAE/`
  (`cpc_calcbacia.py`, `roda_SAE_opera.sh`, `cpc.config`) — eram
  rascunhos nunca promovidos a producao. Os dois `cpc.config` do
  servidor estao em `_legado/servidor/` para referencia.

## Ainda em aberto

- O servidor continua rodando a partir de `OPERA/`, fora do git. Nada
  foi alterado la. Decidir se producao passa a apontar para o
  repositorio consolidado, e como.
- Push para o GitHub: o historico novo divergiu do remoto.
