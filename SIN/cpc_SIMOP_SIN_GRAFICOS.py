import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, DayLocator
import numpy as np
from datetime import datetime, timedelta
import os

# === CONFIGURAÇÕES ===
arquivo = "CHUVA_CPC_COMPLETA_SIN.xlsx"
sheet = "testechuva"
saida_dir = "graficos_bacias"
os.makedirs(saida_dir, exist_ok=True)

# === LEITURA DA PLANILHA ===
df = pd.read_excel(arquivo, sheet_name=sheet)
df.columns = df.columns.str.strip()

# === TRATAMENTO DA DATA ===
# Força ANO, MES, DIA como numérico e renomeia
df["ANO"] = pd.to_numeric(df["ANO"], errors="coerce")
df["MES"] = pd.to_numeric(df["MES"], errors="coerce")
df["DIA"] = pd.to_numeric(df["DIA"], errors="coerce")
df = df.rename(columns={"ANO": "year", "MES": "month", "DIA": "day"})

# Remove linhas com dados ausentes e cria a coluna data
df = df.dropna(subset=["year", "month", "day"])
df["data"] = pd.to_datetime(df[["year", "month", "day"]], errors="coerce")
df = df.dropna(subset=["data"])

# Cria coluna para agrupamento por dia e mês
df["dia_mes"] = df["data"].dt.strftime("%m-%d")

# Detecta as bacias (colunas float que não são data/hora)
colunas_excluir = ["year", "month", "day", "HORA", "SIN"]
bacias = [col for col in df.select_dtypes(include="number").columns if col not in colunas_excluir]

# === DEFINIÇÃO DAS DATAS ===
hoje = datetime.today()
data_ini_atual = hoje - timedelta(days=90)
data_fim_atual = hoje
datas_hoje = pd.date_range(start=data_ini_atual, end=data_fim_atual)

# === LOOP DE GERAÇÃO DE GRÁFICOS ===
for bacia in bacias:
    df_bacia = df[["data", "dia_mes", bacia]].dropna()
    df_bacia[bacia] = pd.to_numeric(df_bacia[bacia], errors="coerce")

    atual = df_bacia[df_bacia["data"].between(data_ini_atual, data_fim_atual)].set_index("data")
    ant = df_bacia[df_bacia["data"].between(data_ini_atual - pd.DateOffset(years=1),
                                            data_fim_atual - pd.DateOffset(years=1))].set_index("data")
    ant.index = ant.index + pd.DateOffset(years=1)
    ant2 = df_bacia[df_bacia["data"].between(data_ini_atual - pd.DateOffset(years=2),
                                             data_fim_atual - pd.DateOffset(years=2))].set_index("data")
    ant2.index = ant2.index + pd.DateOffset(years=2)

    if atual.empty:
        print(f"[AVISO] Sem dados para {bacia}")
        continue

    # === TENDÊNCIA DO ANO ATUAL ===
    x = np.arange(len(atual))
    y = atual[bacia].values
    mask = ~np.isnan(y)
    tendencia = None
    if mask.sum() > 1:
        z = np.polyfit(x[mask], y[mask], 1)
        p = np.poly1d(z)
        tendencia = p(x[mask])

    # === MÁXIMOS HISTÓRICOS DIÁRIOS ===
    maximos = df_bacia.groupby("dia_mes")[bacia].max()
    maximos = maximos.reindex(datas_hoje.strftime("%m-%d")).fillna(0)
    maximos.index = datas_hoje

    # === LEGENDA ===
    l_atual = f"{bacia} {data_fim_atual.year}–{data_ini_atual.year}"
    l_ant = f"{bacia} {data_fim_atual.year -1}–{data_ini_atual.year -1}"
    l_ant2 = f"{bacia} {data_fim_atual.year -2}–{data_ini_atual.year -2}"

    # === PLOT ===
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(atual.index, atual[bacia], color="red", linewidth=2, label=l_atual)
    ax.plot(ant.index, ant[bacia], color="blue", linewidth=2, label=l_ant)
    ax.plot(ant2.index, ant2[bacia], color="green", linewidth=2, label=l_ant2)
    ax.plot(maximos.index, maximos.values, linestyle="dashed", color="gray", alpha=0.7, label="Máximo 1979–2025")
    if tendencia is not None:
        ax.plot(atual.index[mask], tendencia, color="black", linewidth=1.5,
                label=f"Tendência linear {data_ini_atual.year}-{data_fim_atual.year}")

    ax.set_title(f"Bacia {bacia} – Últimos 90 dias")
    ax.set_ylabel("Precipitação (mm)")
    ax.set_ylim(0, 40)
    ax.set_xlabel("Data")
    ax.grid(True, linestyle=":", axis="y")

    # Eixo X com mais datas visíveis
    ax.xaxis.set_major_locator(DayLocator(interval=5))
    ax.xaxis.set_major_formatter(DateFormatter('%d/%m'))

    # Legenda inferior
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)

    plt.tight_layout()
    plt.savefig(f"{saida_dir}/{bacia}.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"\n✅ Gráficos salvos em: {saida_dir}")
