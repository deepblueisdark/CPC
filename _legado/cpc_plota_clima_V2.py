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
from scipy.interpolate import griddata

# ==========================================================
# INTERPOLAÇÃO
# ==========================================================
def interpolar_grade(lat, lon, campo, resolucao):
    lon2d, lat2d = np.meshgrid(lon, lat)

    pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    valores = campo.ravel()

    mask = np.isfinite(valores)
    pontos = pontos[mask]
    valores = valores[mask]

    novo_lon = np.arange(lon.min(), lon.max(), resolucao)
    novo_lat = np.arange(lat.min(), lat.max(), resolucao)

    lon_grid, lat_grid = np.meshgrid(novo_lon, novo_lat)

    campo_interp = griddata(pontos, valores, (lon_grid, lat_grid), method="linear")

    return lat_grid, lon_grid, campo_interp

# ==========================================================
# UNIDADE
# ==========================================================
def fator_conversao_unidade(unidade_origem, unidade_destino):
    if unidade_origem is None:
        return 1.0

    uo = unidade_origem.lower()

    if uo == "mm/day":
        f1 = 1.0
    elif uo == "mm/month":
        f1 = 1 / 30.0
    elif uo == "mm/year":
        f1 = 1 / 365.0
    else:
        return 1.0

    if unidade_destino == "mm/dia":
        f2 = 1.0
    elif unidade_destino == "mm/mes":
        f2 = 30.0
    elif unidade_destino == "mm/ano":
        f2 = 365.0
    else:
        f2 = 1.0

    return f1 * f2

# ==========================================================
# MAPA
# ==========================================================
def gerar_mapa(chuva, lat, lon, titulo, cores, saida_png,
               shapefile=None, limites=None, label_unidade="mm",
               colormap="turbo",
               interpolacao=False, resolucao=0.1,
               politica=None, cor_shapefile="red"):

    if interpolacao:
        lat, lon, chuva = interpolar_grade(lat, lon, chuva, resolucao)

    plt.figure(figsize=(12, 7))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)

    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)
    else:
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=proj)

    ax.coastlines()

    # ===== POLÍTICA =====
    if politica is None:
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.3)

    elif politica[0].lower() == "off":
        pass

    else:
        cor = politica[0]
        lw_b = float(politica[1]) if len(politica) > 1 else 0.5
        lw_s = float(politica[2]) if len(politica) > 2 else 0.3
        ax.add_feature(cfeature.BORDERS, edgecolor=cor, linewidth=lw_b)
        ax.add_feature(cfeature.STATES, edgecolor=cor, linewidth=lw_s)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    norm = BoundaryNorm(cores, ncolors=256)
    cmap = plt.get_cmap(colormap)

    cs = ax.pcolormesh(lon, lat, chuva, cmap=cmap, norm=norm, transform=proj)

    cb = plt.colorbar(
    cs,
    orientation='horizontal',
    pad=0.08,
    aspect=40,
    shrink=0.9
    )

    cb.set_label(label_unidade)

    # 🔥 reduzir quantidade de ticks
    passo = max(1, len(cores)//10)
    ticks = cores##[::passo]

    cb.set_ticks(ticks)

    # 🔥 formatar labels
    labels = [f"{c:.0f}" if c >= 1 else f"{c:.1f}" for c in ticks]
    cb.ax.set_xticklabels(labels)

    # 🔥 fonte menor
    cb.ax.tick_params(labelsize=6)

    # ===== SHAPEFILE =====
    if shapefile and os.path.exists(shapefile):
        gdf = gpd.read_file(shapefile)

        if gdf.crs is None:
            print("[AVISO] Shapefile sem CRS. Assumindo EPSG:4326.")
            gdf.set_crs("EPSG:4326", inplace=True)

        gdf = gdf.to_crs("EPSG:4326")

        gdf.boundary.plot(
            ax=ax,
            edgecolor=cor_shapefile,
            linewidth=1
        )

    plt.title(titulo)
    plt.tight_layout()
    plt.savefig(saida_png, dpi=150, bbox_inches="tight")
    plt.close()

# ==========================================================
# MAIN
# ==========================================================
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("arquivo")

    parser.add_argument("--tipo", choices=["soma", "media", "ambos"], default="soma")
    parser.add_argument("--regiao")
    parser.add_argument("--saida")
    parser.add_argument("--figura", default="mapa")
    parser.add_argument("--titulo", default="Mapa CPC")

    parser.add_argument("--shapefile")

    parser.add_argument("--cores", nargs="+")
    parser.add_argument("--colormap", default="turbo")

    parser.add_argument("--interpolacao", action="store_true")
    parser.add_argument("--resolucao", type=float, default=0.1)

    parser.add_argument("--politica", nargs="+")
    parser.add_argument("--cor-shapefile", default="red")

    parser.add_argument("--clima", action="store_true")

    parser.add_argument("--unidade",
        choices=["mm/dia","mm/mes","mm/ano"],
        default="mm/dia")

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read("cpc.config")

    # ===== REGIÃO =====
    limites = None
    if args.regiao and args.regiao in config:
        reg = config[args.regiao]
        limites = (
            float(reg.get("min_lat")),
            float(reg.get("max_lat")),
            float(reg.get("min_lon")),
            float(reg.get("max_lon")),
        )

    with nc.Dataset(args.arquivo) as ds:
        lat = ds.variables["lat"][:]
        lon = ds.variables["lon"][:]

        tipos = ["soma", "media"] if args.tipo == "ambos" else [args.tipo]
        mapas = {}

        for tipo in tipos:
            if args.clima:
                varname = f"chuva_{tipo}_mes"
            else:
                varname = f"chuva_{tipo}"

            var = ds.variables[varname]
            mapas[tipo] = var[:]

            if tipo == "media":
                unidade_origem = getattr(var, "units", None)

    # ===== CORES =====
    if args.cores:
        if all(c.replace('.', '', 1).isdigit() for c in args.cores):
            cores = list(map(float, args.cores))
        else:
            cores = list(map(float, config["DEFAULT"][args.cores[0]].split()))
    else:
        cores = list(map(float, config["DEFAULT"].get(
            "escala_cores", "0.1 1 2 5 10 20 50 100"
        ).split()))

    # ===== NOMES DOS MESES =====
    # nomes_meses = {
    #     1: "Janeiro", 2: "Fevereiro", 3: "Março", 4: "Abril",
    #     5: "Maio", 6: "Junho", 7: "Julho", 8: "Agosto",
    #     9: "Setembro", 10: "Outubro", 11: "Novembro", 12: "Dezembro"
    # }
    # ===== NOMES DOS MESES =====
    nomes_meses = {
        1: "JANEIRO", 2: "FEVEREIRO", 3: "MARÇO", 4: "ABRIL",
        5: "MAIO", 6: "JUNHO", 7: "JULHO", 8: "AGOSTO",
        9: "SETEMBRO", 10: "OUTUBRO", 11: "NOVEMBRO", 12: "DEZEMBRO"
    }
    # ===== LOOP =====
    for tipo in tipos:
        dados = mapas[tipo]

        if tipo == "media":
            fator = fator_conversao_unidade(unidade_origem, args.unidade)
            dados = dados * fator
            label = args.unidade
        else:
            label = "mm"

        if dados.ndim == 2:
            gerar_mapa(
                dados, lat, lon,
                f"{args.titulo} ({tipo})",
                cores,
                f"{args.figura}_{tipo}.png",
                shapefile=args.shapefile,
                limites=limites,
                colormap=args.colormap,
                interpolacao=args.interpolacao,
                resolucao=args.resolucao,
                politica=args.politica,
                cor_shapefile=args.cor_shapefile,
                label_unidade=label
            )
        else:
            for i in range(dados.shape[0]):
                mes_nome = nomes_meses.get(i+1, f"M{i+1:02d}")
                titulo = f"{mes_nome} - {args.titulo}"

                gerar_mapa(
                    dados[i], lat, lon,
                    titulo,
                    cores,
                    f"{args.figura}_{tipo}_m{i+1:02d}.png",
                    shapefile=args.shapefile,
                    limites=limites,
                    colormap=args.colormap,
                    interpolacao=args.interpolacao,
                    resolucao=args.resolucao,
                    politica=args.politica,
                    cor_shapefile=args.cor_shapefile,
                    label_unidade=label
                )

if __name__ == "__main__":
    main()