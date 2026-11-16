# CPC — Chuva observada (CPC/NOAA) para bacias operacionais

Sistema de aquisicao, calculo por bacia e plotagem de precipitacao
observada a partir dos dados **CPC Unified Gauge-Based Analysis
(NOAA/NCEP)**, usado nas rotinas operacionais SAE, SIN e NFE.

> **Esta pasta e a unica copia de desenvolvimento valida:**
> `<raiz-do-repositorio>`
>
> Consolidada em 2026-11-16 a partir de duas copias que viviam em
> paralelo e estavam ambas incompletas (ver `HISTORICO_CONSOLIDACAO.md`).

---

## Estrutura

```
CPC/
├── cpc.config               Configuracao unica (caminhos, regioes, escalas)
├── OPERACIONAL/             Aquisicao e calculo por bacia
│   ├── cpc_download.py          Baixa e recorta os NetCDF do CPC/NOAA
│   └── cpc_calcbacia.py         Estatisticas por bacia + mapas diarios
├── PLOTAGEM/                Geracao de mapas
│   ├── cpc_plota.py                      Mapa simples por data
│   ├── cpc_gera_plota_chuva_acumulada.py Acumulado entre duas datas -> NetCDF
│   ├── cpc_plota_clima_unificado.py      Plotagem de climatologia (ATIVO)
│   └── cpc_climatology.py                Climatologia (versao enxuta)
├── SAE/                     Rotina SIMOP — bacia do rio Madeira
├── SIN/                     Rotina SIMOP — bacias do SIN
├── NFE/                     Rotina SIMOP — Tocantins / rio Para
├── scripts/
│   └── roda_SAE_opera.sh    Rotina mensal operacional SAE
├── DADOS/                   Planilhas de entrada (fora do git)
├── SAIDAS/                  Tudo que o sistema gera (fora do git)
│   ├── MAPAS_CPC/               Mapas diarios por bacia
│   ├── CLIMA/                   Mapas e NetCDF de climatologia
│   └── PLANILHAS/               Planilhas de saida
└── _legado/                 Versoes antigas mantidas apenas por referencia
```

## Dependencias

```
cartopy  geopandas  shapely  netCDF4  numpy  pandas  scipy
matplotlib  openpyxl
```

As rotinas SIMOP (`SIN/`, `NFE/`) usam `locale pt_BR.UTF-8`.

## Uso

Os scripts localizam o `cpc.config` sozinhos, a partir da raiz do
projeto — podem ser chamados de qualquer diretorio.

**1. Baixar os dados**

```bash
python OPERACIONAL/cpc_download.py -z AMS
```

**2. Calcular por bacia e gerar os mapas diarios**

```bash
python OPERACIONAL/cpc_calcbacia.py --regiao AMS
```

**3. Rotina mensal operacional SAE**

```bash
bash scripts/roda_SAE_opera.sh
# ou, com outro shapefile:
SHAPEFILE=/caminho/SAE_TODAS_AS_BACIAS_HICON_POLY.shp bash scripts/roda_SAE_opera.sh
```

## Configuracao (`cpc.config`)

| Chave | Para que serve |
|---|---|
| `repositorio` | Pasta dos NetCDF baixados do CPC |
| `prefixo` | Prefixo dos arquivos NetCDF |
| `shapefile_dir` / `shapefile_padrao` | Onde estao os shapefiles das bacias |
| `mapa_dir` / `saida_excel` | Destino dos mapas e da planilha |
| `escala_cores` | Escala diaria (mm/dia) |
| `escala_cores_mensal` | Escala mensal/acumulada (mm/mes) |

Regioes disponiveis: `[AMS]`, `[BRASIL]`, `[SAARA]`, `[SAE]`.

**Ao rodar no servidor Linux** (`<servidor>`, conta `simop`):
descomente, em `[DEFAULT]`, as linhas `repositorio` e `shapefile_dir`
com os caminhos `/mnt/...` e comente as versoes Windows.

## Pendencias conhecidas

- **`SAE_TODAS_AS_BACIAS_HICON_POLY.shp` nao esta no repositorio** —
  mora em `shapefile_dir`. O `roda_SAE_opera.sh` falha com mensagem
  clara se nao encontrar.
- **Planilhas referenciadas e ausentes**, usadas pelas rotinas SIMOP:
  `CHUVA_CPC_COMPLETA_SIN.xlsx`, `CHUVA_OBSERVADA.xlsx`,
  `CHUVA_OBSERVADA_CPC_COMPLETA.xlsx`, `chuva_SIN.xlsx`.
  Devem ser colocadas em `DADOS/`.
- **Valores a confirmar no `cpc.config`**: os limites da regiao `[SAE]`
  e a `escala_cores_mensal` foram definidos na consolidacao para
  destravar o `roda_SAE_opera.sh`; revise se batem com a operacao.
- **Versao do servidor nao comparada** — a copia de
  `<servidor>:~simop` ainda nao foi confrontada com esta.
