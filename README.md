# CPC — Chuva observada (CPC/NOAA) para bacias operacionais

Sistema de aquisicao, calculo por bacia e plotagem de precipitacao
observada a partir dos dados **CPC Unified Gauge-Based Analysis
(NOAA/NCEP)**, usado nas rotinas operacionais SAE, SIN e NFE.

> **Esta pasta e a unica copia de desenvolvimento valida:**
> `<raiz-do-repositorio>`
>
> Consolidada em 2026-11-16 a partir de tres fontes que estavam
> todas incompletas: duas copias locais e o servidor de producao
> `<usuario>@<servidor>:~/PRODUCAO/CPC`.
> Ver `HISTORICO_CONSOLIDACAO.md`.
>
> **Escopo: este repositorio guarda codigo, nao operacional.** Os
> scripts das tres aplicacoes (SAE, SIN, NFE) ficam aqui; as planilhas
> acumuladoras, os shapefiles, os mapas gerados e os logs sao dado
> operacional e vivem no servidor. `OPERA/` e `SAIDAS/` existem vazias
> apenas como pontos de montagem da execucao local.

## Como o sistema roda em producao

No servidor, o cron dispara uma rotina unica todo dia as 09:00:

```
<agendamento> /home/<usuario>/PRODUCAO/CPC/OPERA/rodar_cpc.sh >> .../cron.log 2>&1
```

Ela percorre, em ordem, com janela dos **ultimos 90 dias**:

```
cpc_download.py          baixa/atualiza os NetCDF do CPC/NOAA
  -> cpc_calcbacia.py --saida chuva_diaria_CPC.xlsx        (bacias SAE)
     -> cpc_SIMOP_SAE.py
  -> cpc_calcbacia.py --saida chuva_SIN.xlsx  --shapefile SIN
     -> cpc_SIMOP_SIN.py  -> cpc_SIMOP_SIN_GRAFICOS.py
  -> cpc_calcbacia.py --saida chuva_NFE.xlsx  --shapefile NFE
     -> cpc_SIMOP_NFE.py  -> cpc_SIMOP_NFE_GRAFICOS.py
  -> roda_SAE_opera.sh     acumulado mensal SAE + mapas
```

O equivalente local e `scripts/rodar_cpc.sh`.

## Estrutura

```
CPC/
├── cpc.config               Configuracao unica (caminhos, regioes, escalas)
├── OPERACIONAL/             Aquisicao e calculo por bacia
│   ├── cpc_download.py          Baixa e recorta os NetCDF do CPC/NOAA
│   └── cpc_calcbacia.py         Estatisticas por bacia + mapas diarios
├── PLOTAGEM/                Geracao de mapas
│   ├── cpc_plota.py                      Mapa simples por data
│   ├── cpc_gera_plota_chuva_acumulada.py Acumulado entre datas -> NetCDF
│   ├── cpc_plota_clima_unificado.py      Climatologia (ATIVO)
│   ├── cpc_climatology.py                Climatologia (versao enxuta)
│   ├── cpc_calcula_anomalia.py           Anomalia absoluta e percentual
│   └── cpc_plota_anomalia.py             Mapas de anomalia
├── SAE/ SIN/ NFE/           Rotinas SIMOP (planilhas e graficos por bacia)
├── scripts/
│   ├── rodar_cpc.sh             Rotina diaria completa (equiv. ao cron)
│   ├── roda_SAE_opera.sh        Acumulado mensal SAE
│   ├── coleta_servidor.sh       Inventaria o CPC do servidor (roda no WSL)
│   └── sincroniza_opera.sh      Captura o OPERA/ do servidor (roda no WSL)
├── OPERA/                   Diretorio de trabalho, vazio no git
├── SAIDAS/                  Mapas e NetCDF gerados, vazio no git
├── _producao/               Snapshot versionado do OPERA/ do servidor
└── _legado/                 Versoes antigas e copias verbatim do servidor
```

### Por que `OPERA/` e um diretorio plano

As rotinas SIMOP leem e gravam as planilhas acumuladoras por **nome
simples**, no diretorio corrente — `CHUVA_OBSERVADA.xlsx`,
`CHUVA_CPC_COMPLETA_SIN.xlsx`, `chuva_diaria_CPC_NFE.xlsx` e outras.
Elas sao a base historica: entrada e saida ao mesmo tempo. Por isso
`OPERA/` espelha o diretorio homonimo do servidor, e `rodar_cpc.sh`
entra nele antes de chamar os scripts por caminho absoluto.

## Dependencias

```
cartopy  geopandas  shapely  netCDF4  xarray  numpy  pandas  scipy
matplotlib  openpyxl
```

`xarray` e usado so por `cpc_calcula_anomalia.py`.
As rotinas SIN e NFE usam `locale pt_BR.UTF-8`.

## Uso

Os scripts localizam o `cpc.config` sozinhos, a partir da raiz do
projeto — podem ser chamados de qualquer diretorio.

```bash
# rotina diaria completa
bash scripts/rodar_cpc.sh

# so o acumulado mensal SAE + mapas
bash scripts/roda_SAE_opera.sh

# passos isolados
python OPERACIONAL/cpc_download.py -i 01/09/2026 -f 16/11/2026
python OPERACIONAL/cpc_calcbacia.py 01/09/2026 16/11/2026 --saida chuva_diaria_CPC.xlsx
```

## Configuracao (`cpc.config`)

| Chave | Para que serve |
|---|---|
| `repositorio` | Pasta dos NetCDF baixados do CPC |
| `prefixo` | Prefixo dos arquivos NetCDF |
| `shapefile_dir` / `shapefile_padrao` | Shapefiles das bacias (default SAE) |
| `mapa_dir` / `saida_excel` | Destino dos mapas e da planilha |
| `escala_cores` / `escala_diaria` | Escalas diarias (mm/dia) |
| `escala_cores_mensal` | Escala mensal/acumulada (mm/mes) |
| `escala_clima` / `escala_clima_simples` | Escalas de climatologia |

Regioes: `[AMS]`, `[BRASIL]`, `[SAARA]`, `[SAE]`.
Secoes `[ANOMALIA]` e `[ANOMALIA_PER]` alimentam os mapas de anomalia.

As escalas sao os valores que rodam em producao. A regiao `[SAE]` e a
excecao: as latitudes foram corrigidas aqui em 2026-11-16 e **diferem
da producao de proposito** — ver `_producao/DIVERGENCIAS.md`.
**Ao rodar no servidor**, descomente em `[DEFAULT]` as linhas
`repositorio` e `shapefile_dir` com os caminhos `/home/<usuario>/...`.

## O que nao esta aqui, e de onde vem

Por decisao de projeto o repositorio nao carrega dado operacional.
Para rodar a cadeia completa fora do servidor, estes itens precisam
ser providenciados:

| Item | Onde vive | Como apontar |
|---|---|---|
| Planilhas acumuladoras (`CHUVA_OBSERVADA*.xlsx`, `CHUVA_CPC_COMPLETA_SIN.xlsx`, `chuva_diaria_CPC_NFE.xlsx`, `chuva_SIN.xlsx`) | `~/PRODUCAO/CPC/OPERA/` no servidor | copiar para `OPERA/` |
| Shapefiles SAE | `/home/<usuario>/STUFF/SHAPES/SAE/OPERACIONAL/` | `shapefile_dir` no `cpc.config`, ou a variavel `SHAPEFILE` |
| Shapefiles SIN | `/home/<usuario>/STUFF/SHAPES/SIN/` | variavel `SHAPEDIR_SIN` |
| Shapefiles NFE | `/home/<usuario>/STUFF/SHAPES/NFE/` | variavel `SHAPEDIR_NFE` |
| NetCDF do CPC/NOAA | `repositorio` no `cpc.config` | baixados por `cpc_download.py` |

## Sincronizacao com a producao

O servidor continua rodando de `~/PRODUCAO/CPC/OPERA/`, fora do git, e
assim fica. Para que o codigo que roda de fato tenha historico, o
repositorio guarda um snapshot dele em `_producao/OPERA/`:

```bash
# no WSL da maquina de desenvolvimento, nao no servidor
bash scripts/sincroniza_opera.sh
```

A sincronizacao tem **um sentido so** — servidor para repositorio. O
script apenas le do servidor (`ssh` + `tar`), nunca escreve la. Ao
terminar, ele mostra o que mudou na producao desde o ultimo snapshot
commitado. As diferencas intencionais entre as duas copias estao em
`_producao/DIVERGENCIAS.md`; o que aparecer fora daquela lista e drift
de verdade.

## Pendencias abertas

- **`[SAE]`:** as latitudes foram corrigidas neste repositorio em
  2026-11-16 (`min_lat = -21`, `max_lat = -8`), mas **a producao segue
  com os valores invertidos** (`-8` / `-21`), como esta desde 2025-10.
  Antes de levar o config ao servidor, gerar uma imagem do SAE com os
  valores novos e comparar com a do dia: o recorte pode mudar.
- **Push para o GitHub:** o historico local divergiu de `origin/main`
  (`ec25759`); um push seria rejeitado. Decidir entre `--force` e um
  branch novo. O historico antigo esta preservado em
  `CPC_backup_2026-11-16/BACKUP_GITHUB_CPC_2026-11-16.tar.gz`.
- **`__pycache__/cpc_calcbacia.cpython-312.pyc`** nao pode ser removido
  (`Permission denied`, trava do OneDrive). Ignorado pelo git; apagar a
  mao quando o OneDrive liberar.
