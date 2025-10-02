import os
import re
import numpy as np
import netCDF4 as nc
import configparser
from urllib import request
from datetime import datetime, timedelta
import argparse

CONFIG_FILE = "cpc.config"
TEMP_DIR = "./tmp_cpc_bin"

def ler_config(zona="AMS", arquivo="cpc.config", args=None):
    config = configparser.ConfigParser()
    config.read(arquivo)

    if zona not in config:
        print(f"[ERRO] Zona '{zona}' não está definida em {arquivo}")
        exit(1)

    repositorio = args.repositorio or config["DEFAULT"].get("repositorio", "./dados/")
    prefixo = args.prefixo or config["DEFAULT"].get("prefixo", "CPC")

    try:
        min_lat = float(config[zona]["min_lat"])
        max_lat = float(config[zona]["max_lat"])
        min_lon = float(config[zona]["min_lon"])
        max_lon = float(config[zona]["max_lon"])
    except Exception as e:
        print(f"[ERRO] Limites geográficos da zona '{zona}' estão incompletos: {e}")
        exit(1)

    return repositorio, prefixo, min_lat, max_lat, min_lon, max_lon

def baixar_cpc_netcdf(start_date, end_date, repositorio, prefixo, zona, min_lat, max_lat, min_lon, max_lon):
    os.makedirs(TEMP_DIR, exist_ok=True)
    os.makedirs(repositorio, exist_ok=True)

    base_url = "https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_GLB/RT/"
    processados = set()

    for ano in range(start_date.year, end_date.year + 1):
        ano_url = f"{base_url}{ano}/"
        try:
            with request.urlopen(ano_url) as response:
                html = response.read().decode("utf-8")
        except Exception as e:
            print(f"[ERRO] Não foi possível acessar {ano_url}: {e}")
            continue

        arquivos = re.findall(r'PRCP_CU_GAUGE_V1\.0GLB_0\.50deg\.lnx\.(\d{8})\.RT', html)
        for data_str in arquivos:
            if data_str in processados:
                continue

            data_arquivo = datetime.strptime(data_str, "%Y%m%d").date()
            if start_date <= data_arquivo <= end_date:
                ano, mes, dia = data_str[:4], data_str[4:6], data_str[6:8]
                arq_nc = f"{repositorio}/{prefixo}_{ano}{mes}{dia}.nc"

                if os.path.isfile(arq_nc):
                    print(f"[OK] NetCDF já existe: {arq_nc}")
                    processados.add(data_str)
                else:
                    nome_arquivo = f"PRCP_CU_GAUGE_V1.0GLB_0.50deg.lnx.{data_str}.RT"
                    url_download = f"{ano_url}{nome_arquivo}"
                    caminho_tmp = f"{TEMP_DIR}/{nome_arquivo}"
                    try:
                        request.urlretrieve(url_download, caminho_tmp)
                        print(f"[↓] Binário baixado: {nome_arquivo}")
                        sucesso = cria_cpc_netcdf(
                            caminho_tmp, repositorio, prefixo,
                            min_lat=min_lat, max_lat=max_lat,
                            min_lon=min_lon, max_lon=max_lon
                        )
                        if sucesso:
                            print(f"[✔] NetCDF criado: {arq_nc}")
                            processados.add(data_str)
                        else:
                            print(f"[✘] Falha ao criar NetCDF: {arq_nc}")
                        os.remove(caminho_tmp)
                    except Exception as e:
                        print(f"[ERRO] Falha ao baixar ou processar {nome_arquivo}: {e}")

def cria_cpc_netcdf(filename, netcdf_path="", prefixo="", min_lat=-45, max_lat=20, min_lon=280, max_lon=330):
    import netCDF4 as nc
    import numpy as np
    import os
    from datetime import datetime
    import re

    # Extrair data do nome do arquivo
    data_str = re.search(r"(\d{8})", filename).group(1)
    ano, mes, dia = data_str[:4], data_str[4:6], data_str[6:8]
    nc_path = f"{netcdf_path}/{prefixo}_{ano}{mes}{dia}.nc"

    if os.path.exists(nc_path):
        print(f"[CPC CRIA NETCDF] Arquivo {nc_path} já existe. Não será recriado.")
        return nc_path

    # Leitura direta do arquivo binário
    TAMFILE = 2073600
    if os.path.getsize(filename) != TAMFILE:
        print(f"[ERRO] Tamanho incorreto do arquivo binário: {filename}")
        return None

    with open(filename, "rb") as f:
        data = np.frombuffer(f.read(), dtype='<f').reshape((2, 360, 720))
    chuva = data[0] * 0.1
    chuva = np.where(chuva < 0, np.nan, chuva)

    # Coordenadas globais
    lat = np.arange(-89.75, 90, 0.5)
    lon = np.arange(0.25, 360.25, 0.5)
    lat_mask = (lat >= min_lat) & (lat <= max_lat)
    lon_mask = (lon >= min_lon) & (lon <= max_lon)
    lat_subset = lat[lat_mask]
    lon_subset = lon[lon_mask]
    chuva_subset = chuva[lat_mask, :][:, lon_mask]

    # Criar NetCDF
    ncdf = nc.Dataset(nc_path, 'w', format="NETCDF4")

    # Criar dimensões
    ncdf.createDimension('time', None)
    ncdf.createDimension('lat', len(lat_subset))
    ncdf.createDimension('lon', len(lon_subset))

    # Criar variáveis
    lat_var = ncdf.createVariable('lat', 'f4', ('lat',))
    lon_var = ncdf.createVariable('lon', 'f4', ('lon',))
    time_var = ncdf.createVariable('time', 'f8', ('time',))
    chuva_var = ncdf.createVariable('chuva_cpc', 'f4', ('time', 'lat', 'lon'), fill_value=-999.0)

    # Preencher variáveis
    lat_var[:] = lat_subset
    lat_var.units = "degrees_north"
    lat_var.long_name = "latitude"

    lon_var[:] = lon_subset
    lon_var.units = "degrees_east"
    lon_var.long_name = "longitude"

    ref = datetime(1979, 1, 1)
    curr = datetime.strptime(data_str, "%Y%m%d")
    delta_horas = (curr - ref).days * 24.0
    time_var[:] = [delta_horas]
    time_var.units = "hours since 1979-01-01 00:00:0.0"
    time_var.calendar = "standard"
    time_var.long_name = "time"

    chuva_var[0, :, :] = chuva_subset
    chuva_var.units = "mm/day"
    chuva_var.long_name = "CPC_UNI_precipitação diária"

    # Atributos globais
    ncdf.title = "CPC Unified Gauge-Based Daily Precipitation"
    ncdf.institution = "NOAA CPC"
    ncdf.source = "https://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/"
    ncdf.history = f"Criado em {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    ncdf.references = "https://climatedataguide.ucar.edu"

    ncdf.close()
    print(f"[✔] NetCDF  criado: {nc_path}")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Baixa e converte arquivos CPC binários para NetCDF.")
    parser.add_argument("-i", "--inicio", help="Data inicial no formato DD/MM/AAAA")
    parser.add_argument("-f", "--fim", help="Data final no formato DD/MM/AAAA")
    parser.add_argument("-r", "--repositorio", help="Diretório de saída dos NetCDF (sobrescreve o config)")
    parser.add_argument("-p", "--prefixo", help="Prefixo dos arquivos NetCDF (sobrescreve o config)")
    parser.add_argument("-z", "--zona", help="Nome da zona de recorte do cpc.config", default="AMS")

    args = parser.parse_args()

    data_ini = datetime.strptime(args.inicio, "%d/%m/%Y").date() if args.inicio else datetime.now().date() - timedelta(days=7)
    data_fim = datetime.strptime(args.fim, "%d/%m/%Y").date() if args.fim else datetime.now().date()

    if data_ini > data_fim:
        print("[ERRO] A data inicial não pode ser posterior à final.")
        exit()

    zona = args.zona.upper()
    repositorio, prefixo, min_lat, max_lat, min_lon, max_lon = ler_config(zona, CONFIG_FILE, args)
    baixar_cpc_netcdf(data_ini, data_fim, repositorio, prefixo, zona, min_lat, max_lat, min_lon, max_lon)
