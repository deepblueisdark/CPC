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
