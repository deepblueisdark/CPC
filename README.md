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
│   └── coleta_servidor.sh       Inventaria o CPC do servidor (roda no WSL)
├── SHAPES/SAE/              Shapefile das bacias SAE
├── OPERA/                   Diretorio de trabalho (fora do git)
├── SAIDAS/                  Mapas e NetCDF gerados (fora do git)
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

As escalas e a regiao `[SAE]` sao os valores que rodam em producao.
**Ao rodar no servidor**, descomente em `[DEFAULT]` as linhas
`repositorio` e `shapefile_dir` com os caminhos `/home/<usuario>/...`.

## Pendencias conhecidas

- **Quatro planilhas acumuladoras existem so no servidor** e sao
  necessarias para a cadeia SIMOP rodar localmente:
  `CHUVA_OBSERVADA.xlsx`, `CHUVA_OBSERVADA_CPC_COMPLETA.xlsx`,
  `CHUVA_CPC_COMPLETA_SIN.xlsx`, `chuva_SIN.xlsx`.
  Devem ser copiadas para `OPERA/`.
- **Shapes de SIN e NFE nao estao no repositorio** — ficam em
  `/home/<usuario>/STUFF/SHAPES/{SIN,NFE}` no servidor. Ajuste
  `SHAPEDIR_SIN` / `SHAPEDIR_NFE` ao rodar fora dele.
- **`[SAE]` tem `min_lat` (-8) maior que `max_lat` (-21)**, invertido
  em relacao as outras secoes. Esta assim em producao e as imagens
  diarias saem; copiado verbatim, mas vale confirmar.
