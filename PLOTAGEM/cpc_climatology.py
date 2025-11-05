import argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geopandas as gpd
import configparser
import os

def gerar_mapa(chuva, lat, lon, titulo, cores, saida_png, shapefile=None, limites=None):
    plt.figure(figsize=(10, 6))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)
    else:
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=proj)

    ax.coastlines(color='gray', linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='dimgray', linestyle='--')
    ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='darkgray', linestyle=':')

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    norm = BoundaryNorm(cores, ncolors=256)
    cmap = plt.get_cmap("turbo_r")
    cs = ax.pcolormesh(lon, lat, chuva, cmap=cmap, norm=norm, transform=proj)

    cb = plt.colorbar(cs, orientation='horizontal', pad=0.05, aspect=40)
    cb.set_label("mm/day")
    cb.set_ticks(cores)
    cb.ax.set_xticklabels([f"{c:.0f}" if c >= 1 else f"{c:.1f}" for c in cores], fontsize=8)

    if shapefile and os.path.exists(shapefile):
        try:
            gdf = gpd.read_file(shapefile)
            if gdf.crs is None:
                gdf.set_crs("EPSG:4326", inplace=True)
            gdf.to_crs("EPSG:4326").boundary.plot(ax=ax, edgecolor='black', linewidth=1, alpha=0.8)
        except Exception as e:
            print(f"[ERRO] Falha ao ler shapefile: {shapefile}\n{e}")

    plt.title(titulo)
    plt.savefig(saida_png, dpi=150)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Gera mapa mensal de climatologia CPC a partir de um NetCDF")
    parser.add_argument("arquivo", help="Arquivo NetCDF da climatologia")
    parser.add_argument("--mes", type=int, help="Mês (1 a 12) a ser plotado. Se omitido, plota todos os 12 meses")
    parser.add_argument("--cores", nargs='+', help="Escala de cores ou chave definida no cpc.config")
    parser.add_argument("--titulo", default="Climatologia CPC", help="Título base do mapa")
    parser.add_argument("--saida", default=None, help="Arquivo de saída PNG")
    parser.add_argument("--shapefile", help="Shapefile opcional para sobrepor no mapa")
    parser.add_argument("--regiao", help="Nome da região no cpc.config para limitar o mapa")
    parser.add_argument("--fator", type=float, default=1.0, help="Fator multiplicativo aplicado à chuva (default=1.0)")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read("cpc.config")

    with nc.Dataset(args.arquivo) as ds:
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
        chuva_total = ds.variables['chuva'][:] * args.fator

    meses = ["Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho",
             "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro"]

    # Definir limites da região, se houver
    limites = None
    if args.regiao and args.regiao in config:
        reg = config[args.regiao]
        limites = (
            float(reg.get("min_lat", lat.min())),
            float(reg.get("max_lat", lat.max())),
            float(reg.get("min_lon", lon.min())),
            float(reg.get("max_lon", lon.max())),
        )

    # Identificar índices dentro dos limites
    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
        lat_mask = (lat >= min_lat) & (lat <= max_lat)
        lon_mask = (lon >= min_lon) & (lon <= max_lon)
    else:
        lat_mask = np.full_like(lat, True, dtype=bool)
        lon_mask = np.full_like(lon, True, dtype=bool)

    lat_sel = lat[lat_mask]
    lon_sel = lon[lon_mask]

    # Avaliação estatística por mês dentro dos limites
    print("\n[INFO] Estatísticas mensais de precipitação dentro da região selecionada:")
    for i in range(12):
        chuva_mes = chuva_total[np.ix_(lat_mask, lon_mask, [i])][:, :, 0]
        dados_validos = chuva_mes[~np.isnan(chuva_mes)]

        if dados_validos.size > 0:
            valor_min = np.min(dados_validos)
            valor_max = np.max(dados_validos)
            print(f" - {meses[i]:<10}: min = {valor_min:.2f} mm, max = {valor_max:.2f} mm")
        else:
            print(f" - {meses[i]:<10}: apenas valores NaN na região.")



    # Escala de cores
    if args.cores:
        if all(c.replace('.', '', 1).isdigit() for c in args.cores):
            cores = list(map(float, args.cores))
        else:
            chave = args.cores[0]
            if chave in config['DEFAULT']:
                cores = list(map(float, config['DEFAULT'][chave].split()))
            else:
                raise ValueError(f"Chave de cores '{chave}' não encontrada no cpc.config")
    else:
        cores = list(map(float, config['DEFAULT'].get("escala_cores", "0.1 1 2 3 5 10 20 50 100").split()))

    # Escala de cores
    if args.cores:
        if all(c.replace('.', '', 1).isdigit() for c in args.cores):
            cores = list(map(float, args.cores))
        else:
            chave = args.cores[0]
            if chave in config['DEFAULT']:
                cores = list(map(float, config['DEFAULT'][chave].split()))
            else:
                raise ValueError(f"Chave de cores '{chave}' não encontrada no cpc.config")
    else:
        cores = list(map(float, config['DEFAULT'].get("escala_cores", "0.1 1 2 3 5 10 20 50 100").split()))

    limites = None
    if args.regiao and args.regiao in config:
        reg = config[args.regiao]
        limites = (
            float(reg.get("min_lat", lat.min())),
            float(reg.get("max_lat", lat.max())),
            float(reg.get("min_lon", lon.min())),
            float(reg.get("max_lon", lon.max())),
        )

    meses = ["Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho",
             "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro"]

    if args.mes:
        idx = args.mes - 1
        nome_mes = meses[idx]
        titulo = f"{args.titulo} - {nome_mes}"
        saida = args.saida or f"climatologia_{nome_mes.lower()}.png"
        gerar_mapa(chuva_total[:, :, idx], lat, lon, titulo, cores, saida,
                   shapefile=args.shapefile, limites=limites)
    else:
        for idx in range(12):
            nome_mes = meses[idx]
            titulo = f"{args.titulo} - {nome_mes}"
            saida = f"climatologia_{nome_mes.lower()}.png"
            gerar_mapa(chuva_total[:, :, idx], lat, lon, titulo, cores, saida,
                       shapefile=args.shapefile, limites=limites)

if __name__ == "__main__":
    main()
