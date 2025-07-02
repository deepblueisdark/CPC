# cpc_bacias_avancado.py
# Script para calcular estatísticas de chuva por bacia e gerar mapas por ponto classificado para dados CPC

import os
import configparser
import netCDF4 as nc
import numpy as np
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point
from shapely.vectorized import contains
from glob import glob
from datetime import datetime
import argparse
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import argparse
from datetime import datetime


CONFIG_FILE = "cpc.config"
SHAPEFILE_DIR = "/mnt/e/OneDrive/OPERACIONAL/SHAPES/SAE/OPERACIONAL/"
SHAPEFILE_DIR = "/mnt/e/OneDrive/OPERACIONAL/SHAPES/SIN/OPERACIONAL/"

SHAPEFILE_PADRAO = "SIN*.shp"
MAPA_DIR = "./MAPAS_CPC"

def ler_configuracao(regiao="AMS", arquivo_config="cpc.config"):
    """
    Lê o arquivo de configuração e extrai:
    - caminho do repositório dos NetCDFs
    - prefixo dos arquivos
    - escala de cores definida para os mapas
    """
    config = configparser.ConfigParser()
    config.read(arquivo_config)
    escala_raw = config[regiao].get("escala_cores") or config["DEFAULT"].get("escala_cores")
    if escala_raw is None:
        raise ValueError("[ERRO] escala_cores ausente em cpc.config")
    escala_cores = list(map(float, escala_raw.strip().split()))
    caminho = config["DEFAULT"].get("repositorio")
    prefixo = config["DEFAULT"].get("prefixo")
    return caminho, prefixo, escala_cores

def get_colormap():
    """
    Retorna um colormap customizado similar ao estilo IMERG.
    """
    cores = ["#FFFFFF", "#EE82EE", "#0000FF", "#00FF00", "#FFFF00", "#FFA500", "#FF0000"]
    return plt.matplotlib.colors.LinearSegmentedColormap.from_list("ChuvaCPC", cores, N=256)

def gerar_mapa_diario(chuva, lon, lat, shapefile_path, data_str, escala_cores):
    """
    Gera e salva um mapa PNG da chuva diária com contorno das bacias e preenchimento por interpolação.
    """
    os.makedirs(MAPA_DIR, exist_ok=True)
    gdf = gpd.read_file(shapefile_path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    bounds = gdf.total_bounds
    lon_min, lat_min, lon_max, lat_max = bounds

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    cmap = get_colormap()
    norm = BoundaryNorm(escala_cores, ncolors=cmap.N, clip=True)

    p = ax.contourf(lon, lat, chuva, levels=escala_cores, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())

    gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.5, transform=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True)

    ax.set_title(f"Chuva Diária CPC - {data_str}\n{os.path.basename(shapefile_path)}", fontsize=10)
    plt.colorbar(p, ax=ax, orientation='vertical', pad=0.05, label='mm')

    nome_img = f"{MAPA_DIR}/chuva_CPC_{data_str}_{os.path.splitext(os.path.basename(shapefile_path))[0]}.png"
    plt.savefig(nome_img, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[MAPA] Diário salvo: {nome_img}")

def gerar_mapa_classificacao(chuva, lon, lat, shapefile_path, data_str, escala_cores):
    """
    Gera um mapa de classificação de pontos de chuva para cada shapefile:
    - Pontos coloridos por classe de intensidade
    - Pontos fora do polígono (azul) e NaN (vermelho)
    - Título com data e nome do shapefile
    """
    gdf = gpd.read_file(shapefile_path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    bounds = gdf.total_bounds
    lon_min, lat_min, lon_max, lat_max = bounds

    lon_grid, lat_grid = np.meshgrid(lon, lat)
    valores = chuva
    mascara_valida = ~np.isnan(valores)
    mascara_nan = np.isnan(valores)
    geom_union = gdf.unary_union
    mascara_dentro_shape = contains(geom_union, lon_grid, lat_grid)
    validos_dentro = mascara_valida & mascara_dentro_shape
    validos_fora = mascara_valida & ~mascara_dentro_shape

    classes = [
        (0.0, 0.0, 'gray', 'Chuva = 0 mm'),
        (0.0, 0.99, 'lightblue', '0 < Chuva < 1 mm'),
        (1.0, 5.0, 'limegreen', '1–5 mm'),
        (5.0, 10.0, 'yellow', '5–10 mm'),
        (10.0, 15.0, 'orange', '10–15 mm'),
        (15.0, 20.0, 'orangered', '15–20 mm'),
        (20.0, 30.0, 'red', '20–30 mm'),
        (30.0, 40.0, 'darkred', '30–40 mm'),
        (40.0, 50.0, 'purple', '40–50 mm'),
        (50.0, np.inf, 'black', '> 50 mm')
    ]

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon_min, lon_max, lat_min, lat_max])

    for (vmin, vmax, cor, legenda) in classes:
        mask = validos_dentro & (valores >= vmin) & (valores < vmax)
        if np.any(mask):
            ax.scatter(lon_grid[mask], lat_grid[mask], s=8, c=cor, label=legenda, alpha=0.7, transform=ccrs.PlateCarree())

    ax.scatter(lon_grid[validos_fora], lat_grid[validos_fora], s=8, c='blue', label='Fora do shape', alpha=0.5, transform=ccrs.PlateCarree())
    ax.scatter(lon_grid[mascara_nan], lat_grid[mascara_nan], s=8, c='red', label='NaN', alpha=0.5, transform=ccrs.PlateCarree())

    gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=1, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(draw_labels=True)
    ax.set_title(f"Pontos classificados — {data_str}\n{os.path.basename(shapefile_path)}")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8)
    plt.tight_layout()

    os.makedirs(MAPA_DIR, exist_ok=True)
    nome_saida = f"{MAPA_DIR}/classif_CPC_{data_str}_{os.path.splitext(os.path.basename(shapefile_path))[0]}.png"
    plt.savefig(nome_saida, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[MAPA] Classificado salvo: {nome_saida}")

def calcular_estatisticas_bacia(chuva, lon, lat, shapefile_path):
    """
    Calcula estatísticas de chuva por geometria:
    - soma, média, máximo, p90, p95, p99
    - número de pontos válidos, nulos, >=1mm, igual a 0
    Retorna uma lista de resultados por geometria.
    """
    gdf = gpd.read_file(shapefile_path)
    
    resultados = []
    for idx, row in gdf.iterrows():
        try:
            mascara = np.zeros_like(chuva, dtype=bool)
            for i, y in enumerate(lat):
                for j, x in enumerate(lon):
                    if row.geometry.contains(Point(x, y)):
                        mascara[i, j] = True
            recorte = np.where(mascara, chuva, np.nan)
            valores = recorte.flatten()
            total_recorte = np.count_nonzero(mascara)
            n_validos = np.count_nonzero(~np.isnan(valores))
            valores_validos = valores[~np.isnan(valores)]
            n_maior_1mm = np.sum(valores_validos >= 1.0)
            n_zeros = np.sum(valores_validos == 0.0)
            soma = float(np.sum(valores_validos))
            media = float(np.mean(valores_validos)) if n_validos > 0 else np.nan
            maximo = float(np.max(valores_validos)) if n_maior_1mm > 0 else np.nan
            p90 = float(np.percentile(valores_validos, 90)) if n_maior_1mm > 0 else np.nan
            p95 = float(np.percentile(valores_validos, 95)) if n_maior_1mm > 0 else np.nan
            p99 = float(np.percentile(valores_validos, 99)) if n_maior_1mm > 0 else np.nan
            resultados.append([
                shapefile_path, idx, f"shape_{idx}",
                soma, total_recorte, total_recorte - n_validos, n_maior_1mm, n_zeros,
                media, maximo, p90, p95, p99
            ])
        except Exception as e:
            print(f"[ERRO] {shapefile_path} idx {idx}: {e}")
            continue
    return resultados

def processar_chuva_cpc(data_ini, data_fim, regiao="AMS", saida="chuva_diaria_CPC.xlsx", mapas=False, classificacao=False):
    """
    Processa os arquivos CPC em um intervalo de datas:
    - Calcula estatísticas por shapefile
    - Salva planilha com abas "chuva_media" e "chuva_detalhes"
    - (opcional) Gera mapas de classificação por ponto
    """
    repo, prefixo, escala = ler_configuracao(regiao, CONFIG_FILE)
    datas = pd.date_range(start=data_ini, end=data_fim)
    media_geral, detalhes = [], []

    for dia in datas:
        data_str = dia.strftime("%Y%m%d")
        arq = f"{prefixo}_{data_str}.nc"
        caminho_arq = os.path.join(repo, arq)
        if not os.path.isfile(caminho_arq):
            print(f"[AVISO] Arquivo não encontrado: {caminho_arq}")
            continue

        ds = nc.Dataset(caminho_arq)
        chuva = ds.variables['chuva_cpc'][:]
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:] - 360
        ds.close()
        chuva = chuva.squeeze()

        shapefiles = sorted(glob(os.path.join(SHAPEFILE_DIR, SHAPEFILE_PADRAO)))

        for shp in shapefiles:
            print(shp)
            estats = calcular_estatisticas_bacia(chuva, lon, lat, shp)
            for linha in estats:
                nome_shp = os.path.splitext(os.path.basename(linha[0]))[0]
                media_geral.append([nome_shp, linha[1], linha[2], dia.strftime('%Y-%m-%d'), linha[8]])
                detalhes.append([nome_shp, linha[1], linha[2], dia.strftime('%Y-%m-%d')] + linha[3:])

            if classificacao:
                gerar_mapa_classificacao(chuva, lon, lat, shp, dia.strftime('%Y-%m-%d'), escala)
            if mapas:
                gerar_mapa_diario(chuva, lon, lat, shp, dia.strftime('%Y-%m-%d'), escala)

    # Cria DataFrame com todas as estatísticas médias
    df = pd.DataFrame(media_geral, columns=["bacia", "id", "nome", "data", "chuva_mm_media"])
    df['data'] = pd.to_datetime(df['data'])
    df['ANO'] = df['data'].dt.year
    df['MES'] = df['data'].dt.month
    df['DIA'] = df['data'].dt.day

    # Cria tabela pivô com médias
    df_pivot = df.pivot_table(index=['ANO', 'MES', 'DIA'], columns='bacia', values='chuva_mm_media')

    # Identifica todas as bacias existentes
    todas_bacias = sorted(set([linha[0] for linha in detalhes]))

    # Adiciona colunas faltantes com 0.0
    for bacia in todas_bacias:
        if bacia not in df_pivot.columns:
            df_pivot[bacia] = 0.0

    # Reordena colunas
    df_pivot = df_pivot.reset_index()
    colunas_finais = ['ANO', 'MES', 'DIA'] + todas_bacias
    df_pivot = df_pivot[colunas_finais]


    df_detalhes = pd.DataFrame(detalhes, columns=[
        "bacia", "id", "nome", "data", "chuva_mm_soma",
        "n_pontos", "n_nulos", "n_maior_1mm", "n_zeros",
        "chuva_mm_media", "chuva_mm_max", "p90", "p95", "p99"])

    with pd.ExcelWriter(saida) as writer:
        df_pivot.to_excel(writer, index=False, sheet_name="chuva_media")
        df_detalhes.to_excel(writer, index=False, sheet_name="chuva_detalhes")

    print(f"\n✅ Planilha salva como {saida} com abas 'chuva_media' e 'chuva_detalhes'")








# Função que tenta interpretar a data em vários formatos
def parse_data_flexivel(data_str):
    formatos_possiveis = [
        "%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d",
        "%Y.%m.%d", "%d.%m.%Y", "%Y%m%d", "%d%m%Y"
    ]
    for fmt in formatos_possiveis:
        try:
            return datetime.strptime(data_str, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    raise argparse.ArgumentTypeError(f"Formato de data inválido: {data_str}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processa arquivos CPC por bacia e gera mapas.")
    parser.add_argument("data_inicial", type=parse_data_flexivel, help="Data inicial (aceita vários formatos)")
    parser.add_argument("data_final", type=parse_data_flexivel, help="Data final (aceita vários formatos)")
    parser.add_argument("--regiao", default="AMS", help="Região no cpc.config (default: AMS)")
    parser.add_argument("--saida", default="chuva_diaria_CPC.xlsx", help="Arquivo Excel de saída")
    parser.add_argument("--mapas", action="store_true", help="gera mapas preenchidos (contourf)")
    parser.add_argument("--classificacao", action="store_true", help="Gera mapas de classificação por ponto (scatter); se omitido, gera mapas preenchidos (contourf)")
    args = parser.parse_args()

    # Chamada da função principal
    processar_chuva_cpc(
        data_ini=args.data_inicial,
        data_fim=args.data_final,
        regiao=args.regiao,
        saida=args.saida,
        mapas=args.mapas,
        classificacao=args.classificacao
    )

