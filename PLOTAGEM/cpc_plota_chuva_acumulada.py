import os
import argparse
import configparser
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import geopandas as gpd

def ler_data(data_str):
    formatos = ["%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d"]
    for fmt in formatos:
        try:
            return datetime.strptime(data_str, fmt)
        except ValueError:
            continue
    raise ValueError(f"Formato de data inválido: {data_str}")

def ler_configuracao(arquivo="cpc.config", regiao="AMS"):
    config = configparser.ConfigParser()
    config.read(arquivo)

    if regiao not in config:
        raise ValueError(f"Região '{regiao}' não encontrada no arquivo {arquivo}")

    secao = config[regiao]
    repositorio = config["DEFAULT"].get("repositorio")
    prefixo = config["DEFAULT"].get("prefixo")
    cores_str = config["DEFAULT"].get("escala_cores", "0.1 1 2 3 4 5 7 10 15 20 25 30 40 50 75")
    cores = list(map(float, cores_str.split()))

    min_lat = float(secao["min_lat"])
    max_lat = float(secao["max_lat"])
    min_lon = float(secao["min_lon"])
    max_lon = float(secao["max_lon"])

    return config, repositorio, prefixo, cores, (min_lat, max_lat, min_lon, max_lon)

def carregar_dados_periodo(repositorio, prefixo, data_ini, data_fim):
    delta = timedelta(days=1)
    lista_chuvas = []
    lat = lon = None

    dt = data_ini
    while dt <= data_fim:
        nome = f"{repositorio}/{prefixo}_{dt.strftime('%Y%m%d')}.nc"
        if os.path.exists(nome):
            with nc.Dataset(nome) as ds:
                if "chuva_cpc" not in ds.variables:
                    raise KeyError(f"Variável 'chuva_cpc' não encontrada em {nome}")
                chuva = ds.variables["chuva_cpc"][0, :, :]
                if lat is None:
                    lat = ds.variables["lat"][:]
                    lon = ds.variables["lon"][:]
                lista_chuvas.append(chuva)
        else:
            print(f"[AVISO] Arquivo não encontrado: {nome}")
        dt += delta

    if not lista_chuvas:
        raise ValueError("Nenhum arquivo de chuva foi carregado.")

    pilha = np.stack(lista_chuvas)
    return lat, lon, pilha

def salvar_netcdf_saida(lat, lon, chuva_soma, chuva_media, data_ini, data_fim, saida="chuva_acumulada.nc"):
    with nc.Dataset(saida, "w") as ds:
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        ds.createVariable("lat", "f4", ("lat",))[:] = lat
        ds.createVariable("lon", "f4", ("lon",))[:] = lon
        ds.createVariable("chuva_soma", "f4", ("lat", "lon"))[:] = chuva_soma
        ds.createVariable("chuva_media", "f4", ("lat", "lon"))[:] = chuva_media

        ds.data_inicial = data_ini.strftime("%Y-%m-%d")
        ds.data_final = data_fim.strftime("%Y-%m-%d")
        ds.descricao = "Soma e média de precipitação do CPC no período"

def interpolar_grade(lat, lon, campo, resolucao=0.1):
    lon2d, lat2d = np.meshgrid(lon, lat)
    pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    valores = campo.ravel()

    lon_interp = np.arange(lon.min(), lon.max(), resolucao)
    lat_interp = np.arange(lat.min(), lat.max(), resolucao)
    lon_i, lat_i = np.meshgrid(lon_interp, lat_interp)

    campo_interp = griddata(pontos, valores, (lon_i, lat_i), method="linear")
    return lat_interp, lon_interp, campo_interp

def gerar_mapa(chuva, lat, lon, titulo, cores, saida_png, limites,
               shapefile=None, colormap="turbo", politica=["black", "0.5", "0.3"]):
    min_lat, max_lat, min_lon, max_lon = limites
    plt.figure(figsize=(10, 6))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)

    ax.coastlines()

    # ---- NOVO BLOCO: controle das fronteiras políticas ----
    if not (len(politica) == 1 and politica[0].lower() == "off"):
        try:
            cor = politica[0]
            lw_borda = float(politica[1]) if len(politica) > 1 else 0.5
            lw_estado = float(politica[2]) if len(politica) > 2 else 0.3
            ax.add_feature(cfeature.BORDERS, edgecolor=cor, linewidth=lw_borda)
            ax.add_feature(cfeature.STATES, edgecolor=cor, linewidth=lw_estado)
        except Exception as e:
            print(f"[AVISO] Erro ao aplicar opção --politica: {politica} -> {e}")
    # -------------------------------------------------------

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


    norm = BoundaryNorm(cores, ncolors=256)
    cmap = plt.get_cmap(colormap)
    cs = ax.pcolormesh(lon, lat, chuva, cmap=cmap, norm=norm, transform=proj)

    cb = plt.colorbar(cs, orientation='horizontal', pad=0.05, aspect=40)
    cb.set_label("mm")
    cb.set_ticks(cores)
    cb.ax.set_xticklabels([f"{c:.0f}" if c >= 1 else f"{c:.1f}" for c in cores], fontsize=8)

    if shapefile and os.path.exists(shapefile):
        try:
            gdf = gpd.read_file(shapefile)
            if gdf.crs is None:
                gdf.set_crs("EPSG:4326", inplace=True)
            gdf.to_crs("EPSG:4326").boundary.plot(ax=ax, edgecolor='red', linewidth=1, alpha=0.8)
        except Exception as e:
            print(f"[ERRO] Falha ao ler shapefile: {shapefile}\n{e}")

    plt.title(titulo)
    plt.savefig(saida_png, dpi=150)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gera NetCDF com soma/média da chuva e mapas.")
    parser.add_argument("data_inicial", help="Data inicial (formato flexível)")
    parser.add_argument("data_final", help="Data final (formato flexível)")
    parser.add_argument("--regiao", default="AMS", help="Região definida no cpc.config")
    parser.add_argument("--saida", default="chuva_acumulada.nc", help="Arquivo NetCDF de saída")
    parser.add_argument("--interpolacao", action="store_true", help="Aplica interpolação para suavizar os mapas")
    parser.add_argument("--resolucao", type=float, default=0.1, help="Resolução da interpolação (default: 0.1)")
    parser.add_argument("--shapefile", help="Arquivo shapefile a ser plotado no mapa")
    parser.add_argument("--cores", nargs="+",
                        help="Escala de cores ou chave do cpc.config (ex: --cores 0 10 20 ou --cores escala_verao)")
    parser.add_argument("--figura", default="", help="Prefixo para nome dos arquivos gerados")
    parser.add_argument("--titulo", default="", help="Título exibido na figura antes do intervalo de datas")
    parser.add_argument("--colormap", default="turbo", help="Nome do colormap do Matplotlib (padrão: turbo)")
    parser.add_argument("--politica", nargs="+", default=["black", "0.5", "0.3"],
                        help="Desenha fronteiras políticas. Use 'off' para desligar ou 'cor largura_borda largura_estado' (ex: --politica black 0.5 0.3)")
    args = parser.parse_args()
    # Validação opcional do colormap
    try:
        plt.get_cmap(args.colormap)
    except ValueError:
        print(f"[AVISO] Colormap '{args.colormap}' não encontrado. Usando 'turbo'.")
        args.colormap = "turbo"

    data_ini = ler_data(args.data_inicial)
    data_fim = ler_data(args.data_final)

    config, repositorio, prefixo, cores, limites = ler_configuracao(regiao=args.regiao)

    if args.cores:
        if all(c.replace('.', '', 1).isdigit() for c in args.cores):
            cores = list(map(float, args.cores))
        else:
            chave = args.cores[0]
            if chave in config["DEFAULT"]:
                cores = list(map(float, config["DEFAULT"][chave].split()))
            else:
                raise ValueError(f"A chave '{chave}' não foi encontrada no cpc.config.")

    lat, lon, pilha = carregar_dados_periodo(repositorio, prefixo, data_ini, data_fim)

    chuva_soma = np.sum(pilha, axis=0)
    chuva_media = np.mean(pilha, axis=0)

    salvar_netcdf_saida(lat, lon, chuva_soma, chuva_media, data_ini, data_fim, saida=args.saida)

    data_ini_str = data_ini.strftime("%Y%m%d")
    data_fim_str = data_fim.strftime("%Y%m%d")
    nome_soma = f"{args.figura}_soma_{data_ini_str}_{data_fim_str}.png"
    nome_media = f"{args.figura}_media_{data_ini_str}_{data_fim_str}.png"

    if args.interpolacao:
        lat_i, lon_i, soma_i = interpolar_grade(lat, lon, chuva_soma, resolucao=args.resolucao)
        lat_j, lon_j, media_j = interpolar_grade(lat, lon, chuva_media, resolucao=args.resolucao)
        gerar_mapa(soma_i, lat_i, lon_i, f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d}",
                   cores, nome_soma, limites, shapefile=args.shapefile, colormap=args.colormap,politica=args.politica)
        gerar_mapa(media_j, lat_j, lon_j, f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d} (média)",
                   cores, nome_media, limites, shapefile=args.shapefile, colormap=args.colormap,politica=args.politica)
    else:
        gerar_mapa(chuva_soma, lat, lon, f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d}",
                   cores, nome_soma, limites, shapefile=args.shapefile, colormap=args.colormap,politica=args.politica)
        gerar_mapa(chuva_media, lat, lon, f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d} (média)",
                   cores, nome_media, limites, shapefile=args.shapefile, colormap=args.colormap,politica=args.politica)
