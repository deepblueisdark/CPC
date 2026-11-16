import argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geopandas as gpd
import configparser
import os

# ==========================================================
# Resolucao de caminhos - adicionado na consolidacao 2026-11-16.
# Antes, o cpc.config era procurado apenas no diretorio atual,
# o que impedia rodar o script de fora da raiz do projeto.
# Agora ele e procurado na pasta do script e nas pastas acima.
# ==========================================================
import os as _os


def _raiz_projeto(nome="cpc.config", niveis=4):
    _d = _os.path.dirname(_os.path.abspath(__file__))
    for _ in range(niveis):
        if _os.path.isfile(_os.path.join(_d, nome)):
            return _d
        _pai = _os.path.dirname(_d)
        if _pai == _d:
            break
        _d = _pai
    return _os.getcwd()


RAIZ_PROJETO = _raiz_projeto()
CONFIG_FILE = _os.path.join(RAIZ_PROJETO, "cpc.config")
# ==========================================================


# ==========================
# REGRID
# ==========================
def regrid_data(lat, lon, data, resolucao):
    from scipy.interpolate import griddata

    lat_new = np.arange(lat.min(), lat.max(), resolucao)
    lon_new = np.arange(lon.min(), lon.max(), resolucao)

    lon2d, lat2d = np.meshgrid(lon, lat)
    lon2d_new, lat2d_new = np.meshgrid(lon_new, lat_new)

    pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    valores = data.ravel()

    mask = np.isfinite(valores)

    pontos = pontos[mask]
    valores = valores[mask]

    dados_interp = griddata(pontos, valores,
                            (lon2d_new, lat2d_new),
                            method="linear")

    return lat_new, lon_new, dados_interp


# ==========================
# CORES CONFIG
# ==========================
def ler_cores_config(config, nome_secao):
    secao = config[nome_secao]
    valores = secao.get("valores") or secao.get("cores")
    return [float(v) for v in valores.split()]


# ==========================
# CLASSES %
# ==========================
def obter_classes_percentual():
    classes = [-999, -75, -50, -25, -10, 10, 25, 50, 75, 999]

    cores = ["#000000","#ff0000","#ff9900","#ffff00",
             "#e0e0e0","#99cc66","#00aa44","#33aaff","#7a1fa2"]

    labels = ["Extremo Seco","Muito Seco","Seco","Seco-Normal",
              "Normal","Normal-Chuvoso","Chuvoso","Muito Chuvoso","Extremo"]

    return classes, cores, labels


# ==========================
# PLOT
# ==========================
def gerar_mapa(anomalia, lat, lon, titulo, cores, saida_png,
               shapefile=None, limites=None,
               inverter_cores=False,
               classes=None, labels=None,
               politica=False,
               cor_shapefile="black",
               interpolacao="auto"):

    plt.figure(figsize=(10, 6))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)

    if limites:
        ax.set_extent([limites[2], limites[3], limites[0], limites[1]])
    else:
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])

    ax.coastlines()

    if not politica:
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.3)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # ==========================
    # CATEGÓRICO
    # ==========================
    if classes is not None:

        cmap = ListedColormap(cores)
        norm = BoundaryNorm(np.arange(len(cores)+1)-0.5, len(cores))

        dados_class = np.full_like(anomalia, np.nan)

        for i in range(len(classes)-1):
            mask = (anomalia >= classes[i]) & (anomalia < classes[i+1])
            dados_class[mask] = i

        cs = ax.pcolormesh(
            lon, lat, dados_class,
            cmap=cmap, norm=norm,
            shading=interpolacao,
            transform=proj
        )

        cb = plt.colorbar(cs, orientation='horizontal', pad=0.05)
        cb.set_ticks(range(len(labels)))
        cb.set_ticklabels(labels)
        cb.ax.tick_params(labelsize=8)

    # ==========================
    # CONTÍNUO
    # ==========================
    else:
        cmap = plt.get_cmap("RdBu_r")
        if inverter_cores:
            cmap = cmap.reversed()

        norm = BoundaryNorm(cores, 256)

        cs = ax.pcolormesh(
            lon, lat, anomalia,
            cmap=cmap, norm=norm,
            shading=interpolacao,
            transform=proj
        )

        cb = plt.colorbar(cs, orientation='horizontal', pad=0.05)
        cb.set_ticks(cores)
        cb.ax.set_xticklabels([f"{c:.0f}" for c in cores], fontsize=8)

    # shapefile
    if shapefile and os.path.exists(shapefile):
        gdf = gpd.read_file(shapefile)
        gdf.boundary.plot(ax=ax, edgecolor=cor_shapefile, linewidth=1)

    plt.title(titulo)
    plt.savefig(saida_png, dpi=150, bbox_inches="tight")
    plt.close()


# ==========================
# MAIN
# ==========================
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("arquivo")
    parser.add_argument("--saida", default="anomalia.png")
    parser.add_argument("--saida-dir", default=".")
    parser.add_argument("--titulo", default="Anomalia de Precipitação")
    parser.add_argument("--shapefile")
    parser.add_argument("--regiao")
    parser.add_argument("--auto-escala", action="store_true")
    parser.add_argument("--cores")
    parser.add_argument("--inverter-cores", action="store_true")
    parser.add_argument("--classes", action="store_true")
    parser.add_argument("--politica", action="store_true")
    parser.add_argument("--cor-shapefile", default="black")
    parser.add_argument("--interpolacao", default="auto")
    parser.add_argument("--resolucao", type=float)

    args = parser.parse_args()

    os.makedirs(args.saida_dir, exist_ok=True)

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    with nc.Dataset(args.arquivo) as ds:
        lat = ds.variables["lat"][:]
        lon = ds.variables["lon"][:]
        anom = ds.variables["anomalia"][:]
        anom_perc = ds.variables["anomalia_percentual"][:]

    # limpeza
    anom = np.array(anom, float)
    anom[~np.isfinite(anom)] = np.nan

    anom_perc = np.array(anom_perc, float)
    anom_perc[~np.isfinite(anom_perc)] = np.nan

    # ==========================
    # REGRID
    # ==========================
    if args.resolucao:

        lat_orig = lat.copy()
        lon_orig = lon.copy()

        lat_new, lon_new, anom = regrid_data(lat_orig, lon_orig, anom, args.resolucao)
        _, _, anom_perc = regrid_data(lat_orig, lon_orig, anom_perc, args.resolucao)

        lat = lat_new
        lon = lon_new

    # ==========================
    # CORES
    # ==========================
    if args.cores:
        cores = ler_cores_config(config, args.cores)
    elif args.auto_escala:
        vmax = np.nanpercentile(np.abs(anom), 98)
        cores = np.linspace(-vmax, vmax, 15).tolist()
    else:
        cores = [-300,-200,-150,-100,-50,-20,-10,0,10,20,50,100,150,200,300]

    # ==========================
    # REGIÃO
    # ==========================
    limites = None
    if args.regiao and args.regiao in config:
        r = config[args.regiao]
        limites = (float(r["min_lat"]), float(r["max_lat"]),
                   float(r["min_lon"]), float(r["max_lon"]))

    # ==========================
    # MAPA ABSOLUTO
    # ==========================
    saida_png = os.path.join(args.saida_dir, args.saida)

    gerar_mapa(
        anom, lat, lon,
        args.titulo,
        cores,
        saida_png,
        shapefile=args.shapefile,
        limites=limites,
        inverter_cores=args.inverter_cores,
        politica=args.politica,
        cor_shapefile=args.cor_shapefile,
        interpolacao=args.interpolacao
    )

    print("[OK] mapa absoluto")

    # ==========================
    # MAPA %
    # ==========================
    saida_png_perc = saida_png.replace(".png", "_percentual.png")

    if args.classes:
        classes, cores_perc, labels = obter_classes_percentual()

        gerar_mapa(
            anom_perc, lat, lon,
            args.titulo + " (%)",
            cores_perc,
            saida_png_perc,
            shapefile=args.shapefile,
            limites=limites,
            classes=classes,
            labels=labels,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile,
            interpolacao=args.interpolacao
        )

    else:
        gerar_mapa(
            anom_perc, lat, lon,
            args.titulo + " (%)",
            cores,
            saida_png_perc,
            shapefile=args.shapefile,
            limites=limites,
            inverter_cores=args.inverter_cores,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile,
            interpolacao=args.interpolacao
        )

    print("[OK] mapa percentual")


if __name__ == "__main__":
    main()