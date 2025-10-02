import os
import argparse
import configparser
from datetime import datetime, timedelta
from urllib.request import urlretrieve
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import cm
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib import colormaps as mpl_colormaps


def carregar_config(config_path='cpc.config', regiao='AMS'):
    config = configparser.ConfigParser()
    config.read(config_path)

    if regiao not in config:
        raise ValueError(f"Região '{regiao}' não encontrada em {config_path}")

    padrao = config['DEFAULT']
    reg = config[regiao]

    repositorio = padrao.get('repositorio', './HISTORICO/')
    prefixo = padrao.get('prefixo', 'CPC_GAUGE')

    escala_str = padrao.get("escala_cores", "0 1 2 3 4 5 7 10 15 20 25 30 40 50 75")
    escala_cores = list(map(float, escala_str.strip().split()))

    return {
        'repositorio': repositorio,
        'prefixo': prefixo,
        'min_lat': float(reg['min_lat']),
        'max_lat': float(reg['max_lat']),
        'min_lon': float(reg['min_lon']),
        'max_lon': float(reg['max_lon']),
        'regiao': regiao,
        'escala_cores': escala_cores,
    }


import netCDF4 as nc

def gerar_colormap_personalizado(niveis):
    if niveis[0] <= 0:
        niveis = [0.1] + niveis[1:]
    cmap_base = mpl_colormaps["rainbow"].resampled(len(niveis) - 1)
    cmap = ListedColormap(cmap_base(np.linspace(0, 1, len(niveis) - 1)))
    cmap.set_bad((1, 1, 1, 0))  # branco totalmente transparente        # NaN
    cmap.set_under("lightgray")  # 0.0
    norm = BoundaryNorm(niveis, cmap.N)
    return cmap, norm


def plot_mapa(chuva, lat, lon, path_png, title,
              min_lat=None, max_lat=None, min_lon=None, max_lon=None,
              escala_cores=None):

    if lon.ndim == 1 and lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lon2d, lat2d = lon, lat

    if lon2d.shape != chuva.shape:
        try:
            ordem = np.argsort(lon)
            lon = lon[ordem]
            chuva = chuva[:, ordem]
            lon2d, lat2d = np.meshgrid(lon, lat)
        except Exception as e:
            raise ValueError("Dimensões incompatíveis entre lon/lat e a matriz de chuva.") from e

    chuva_masked = np.ma.masked_invalid(chuva)
    cmap, norm = gerar_colormap_personalizado(escala_cores or [0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75])

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': proj})

    img = ax.pcolormesh(lon2d, lat2d, chuva_masked, cmap=cmap, norm=norm, shading='auto', transform=proj)

    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.5)

    if all(v is not None for v in [min_lon, max_lon, min_lat, max_lat]):
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)

    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.colorbar(img, ax=ax, orientation="vertical", label="mm/dia")
    plt.tight_layout()
    plt.savefig(path_png)
    plt.close()
    print(f"[IMAGEM] {path_png}")

from scipy.interpolate import griddata

def visualizar_netcdf(file_path, config, resolucao_interp=None, metodo_interp='linear'):
    ds = nc.Dataset(file_path, 'r')
    chuva = ds.variables['chuva_cpc'][:]
    lat = ds.variables['lat'][:]
    if chuva.ndim == 3:
        chuva = chuva[0, :, :]
    lon = ds.variables['lon'][:]
    ds.close()

    if resolucao_interp:
        lon2d, lat2d = np.meshgrid(lon, lat)
        pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
        valores = chuva.ravel()

        new_lon = np.arange(lon.min(), lon.max() + resolucao_interp, resolucao_interp)
        new_lat = np.arange(lat.min(), lat.max() + resolucao_interp, resolucao_interp)
        new_lon2d, new_lat2d = np.meshgrid(new_lon, new_lat)

        chuva_interp = griddata(pontos, valores, (new_lon2d, new_lat2d), method=metodo_interp)
        lon, lat, chuva = new_lon, new_lat, chuva_interp

    nome_base = os.path.splitext(os.path.basename(file_path))[0]
    titulo = f"CPC - {config['regiao']} {nome_base[-8:]}"

    plot_mapa(chuva, lat, lon, f"{nome_base}.png", titulo,
              min_lat=config['min_lat'], max_lat=config['max_lat'],
              min_lon=config['min_lon'], max_lon=config['max_lon'],
              escala_cores=config['escala_cores'])

def plot_netcdf_por_data(data_inicial, data_final, config, resolucao_interp=None, metodo_interp='linear'):
    dt_ini = datetime.strptime(data_inicial, "%d/%m/%Y")
    dt_fim = datetime.strptime(data_final, "%d/%m/%Y") if data_final else dt_ini

    for data in [dt_ini + timedelta(days=i) for i in range((dt_fim - dt_ini).days + 1)]:
        nome_arquivo = f"{config['prefixo']}_{data.strftime('%Y%m%d')}.nc"
        caminho_arquivo = os.path.join(config['repositorio'], nome_arquivo)
        if not os.path.isfile(caminho_arquivo):
            print(f"[AVISO] Arquivo não encontrado: {caminho_arquivo}")
            continue
        visualizar_netcdf(caminho_arquivo, config, resolucao_interp, metodo_interp)

def main():
    parser = argparse.ArgumentParser(description="Plotar imagens do CPC binário ou NetCDF.")
    parser.add_argument("--interp", type=float, default=None, help="Interpolar NetCDF para nova resolução (ex: 0.1)")
    parser.add_argument("--method", type=str, default="linear", choices=["linear", "nearest", "cubic"], help="Método de interpolação (linear, nearest, cubic)")
    parser.add_argument("data_inicial", help="Data inicial (DD/MM/AAAA)")
    parser.add_argument("data_final", nargs="?", help="Data final (opcional)")
    parser.add_argument("-r", "--regiao", default="AMS", help="Região definida em cpc.config")
    parser.add_argument("-n", "--netcdf", action="store_true", help="Usar arquivos NetCDF no lugar dos binários")
    args = parser.parse_args()

    config = carregar_config(regiao=args.regiao)

    if args.netcdf:
        plot_netcdf_por_data(args.data_inicial, args.data_final, config, resolucao_interp=args.interp, metodo_interp=args.method)
    else:
        # Modo binário
        for data in [datetime.strptime(args.data_inicial, "%d/%m/%Y") + timedelta(days=i)
                     for i in range((datetime.strptime(args.data_final, "%d/%m/%Y") - datetime.strptime(args.data_inicial, "%d/%m/%Y")).days + 1 if args.data_final else 1)]:
            data_str = data.strftime('%Y%m%d')
            nome_bin = f"PRCP_CU_GAUGE_V1.0GLB_0.50deg.lnx.{data_str}.RT"
            url = f"https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_GLB/RT/{data.year}/{nome_bin}"
            local_bin = os.path.join("./temp", nome_bin)

            if not os.path.exists("./temp"):
                os.makedirs("./temp")

            try:
                print(f"[DOWNLOAD] {url}")
                urlretrieve(url, local_bin)
            except Exception as e:
                print(f"[ERRO] Falha no download: {e}")
                continue

            try:
                with open(local_bin, 'rb') as f:
                    raw = np.fromfile(f, dtype='<f').reshape(2, 360, 720)[0] * 0.1
                raw[raw < 0] = np.nan
                lat = np.arange(-89.75, 90, 0.5)
                lon = np.arange(0.25, 360.25, 0.5)

                nome_base = f"CPC_GAUGE_BIN_{data_str}"
                titulo = f"CPC BINÁRIO - {config['regiao']} {data.strftime('%d/%m/%Y')}"

                lon_corrigido = np.where(lon > 180, lon - 360, lon)
                ordem = np.argsort(lon_corrigido)
                lon_corrigido = lon_corrigido[ordem]
                raw_ordenado = raw[:, ordem]

                plot_mapa(raw_ordenado, lat, lon_corrigido, f"{nome_base}_global.png", titulo,
                          escala_cores=config['escala_cores'])

                # mapa recortado
                plot_mapa(raw, lat, lon, f"{nome_base}_{config['regiao']}.png", titulo,
                          min_lat=config['min_lat'], max_lat=config['max_lat'],
                          min_lon=config['min_lon'], max_lon=config['max_lon'],
                          escala_cores=config['escala_cores'])
            except Exception as e:
                print(f"[ERRO] Falha ao processar {nome_bin}: {e}")

if __name__ == "__main__":
    main()
