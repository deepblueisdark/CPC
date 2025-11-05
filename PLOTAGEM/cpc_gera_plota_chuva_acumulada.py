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
    """
    Lê os arquivos diários de chuva entre data_ini e data_fim.

    Retorna:
        lat, lon, pilha (time, lat, lon), datas (lista de datetime para cada slice da pilha)
    """
    delta = timedelta(days=1)
    lista_chuvas = []
    datas = []
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
                datas.append(dt)
        else:
            print(f"[AVISO] Arquivo não encontrado: {nome}")
        dt += delta

    if not lista_chuvas:
        raise ValueError("Nenhum arquivo de chuva foi carregado.")

    pilha = np.stack(lista_chuvas)
    return lat, lon, pilha, datas

def salvar_netcdf_saida(
    lat,
    lon,
    chuva_soma,
    chuva_media,
    data_ini,
    data_fim,
    saida="chuva_acumulada.nc",
    unidade_media="mm/day",
    soma_mes=None,
    media_mes=None,
    cont_mes=None
):
    """
    Salva NetCDF com:
      - chuva_soma (lat, lon)
      - chuva_media (lat, lon)
    E, opcionalmente (quando soma_mes/media_mes não forem None):
      - chuva_soma_mes (month, lat, lon)
      - chuva_media_mes (month, lat, lon)
      - n_dias_mes (month)
    """
    with nc.Dataset(saida, "w") as ds:
        # Dimensões básicas
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        var_lat = ds.createVariable("lat", "f4", ("lat",))
        var_lon = ds.createVariable("lon", "f4", ("lon",))
        var_lat[:] = lat
        var_lon[:] = lon
        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"

        # Variáveis de saída principais
        var_soma = ds.createVariable("chuva_soma", "f4", ("lat", "lon"))
        var_soma[:] = chuva_soma
        var_soma.units = "mm"
        var_soma.long_name = "Soma da precipitação no período"

        var_media = ds.createVariable("chuva_media", "f4", ("lat", "lon"))
        var_media[:] = chuva_media
        var_media.units = unidade_media
        var_media.long_name = "Média da precipitação no período"

        # Climatologia mensal (opcional)
        if soma_mes is not None and media_mes is not None and cont_mes is not None:
            ds.createDimension("month", 12)
            var_month = ds.createVariable("month", "i4", ("month",))
            var_month[:] = np.arange(1, 13)
            var_month.long_name = "Mês do ano (1-12)"

            var_soma_mes = ds.createVariable("chuva_soma_mes", "f4", ("month", "lat", "lon"))
            var_soma_mes[:] = soma_mes
            var_soma_mes.units = "mm"
            var_soma_mes.long_name = (
                "Soma da precipitação em todos os dias de cada mês ao longo do período"
            )

            var_media_mes = ds.createVariable("chuva_media_mes", "f4", ("month", "lat", "lon"))
            var_media_mes[:] = media_mes
            var_media_mes.units = "mm/day"
            var_media_mes.long_name = (
                "Média diária da precipitação em cada mês ao longo do período"
            )

            var_n_dias = ds.createVariable("n_dias_mes", "i4", ("month",))
            var_n_dias[:] = cont_mes
            var_n_dias.long_name = "Número de dias usados em cada mês na climatologia"

        ds.data_inicial = data_ini.strftime("%Y-%m-%d")
        ds.data_final = data_fim.strftime("%Y-%m-%d")
        ds.descricao = (
            "Soma e média de precipitação do CPC no período; "
            "inclui climatologia mensal se calculada."
        )

def interpolar_grade(lat, lon, campo, resolucao=0.1):
    lon2d, lat2d = np.meshgrid(lon, lat)
    pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    valores = campo.ravel()

    lon_interp = np.arange(lon.min(), lon.max(), resolucao)
    lat_interp = np.arange(lat.min(), lat.max(), resolucao)
    lon_i, lat_i = np.meshgrid(lon_interp, lat_interp)

    campo_interp = griddata(pontos, valores, (lon_i, lat_i), method="linear")
    return lat_interp, lon_interp, campo_interp

def gerar_mapa(
    chuva,
    lat,
    lon,
    titulo,
    cores,
    saida_png,
    limites,
    shapefile=None,
    colormap="turbo",
    politica=["black", "0.5", "0.3"],
    cor_shapefile=None
):
    min_lat, max_lat, min_lon, max_lon = limites
    plt.figure(figsize=(10, 6))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)

    ax.coastlines()

    # Controle de fronteiras políticas via --politica
    if not (len(politica) == 1 and politica[0].lower() == "off"):
        try:
            cor = politica[0]
            lw_borda = float(politica[1]) if len(politica) > 1 else 0.5
            lw_estado = float(politica[2]) if len(politica) > 2 else 0.3
            ax.add_feature(cfeature.BORDERS, edgecolor=cor, linewidth=lw_borda)
            ax.add_feature(cfeature.STATES, edgecolor=cor, linewidth=lw_estado)
        except Exception as e:
            print(f"[AVISO] Erro ao aplicar opção --politica: {politica} -> {e}")

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
            gdf = gdf.to_crs("EPSG:4326")
            edge_color = cor_shapefile if cor_shapefile is not None else 'red'
            gdf.boundary.plot(ax=ax, edgecolor=edge_color, linewidth=1, alpha=0.8)
        except Exception as e:
            print(f"[ERRO] Falha ao ler shapefile: {shapefile}\n{e}")

    plt.title(titulo)
    plt.savefig(saida_png, dpi=150)
    plt.close()

def calcular_climatologia_mensal(pilha, datas):
    """
    Calcula soma e média diária por mês do ano ao longo do período.

    Entradas:
      - pilha: array (time, lat, lon) com chuva diária (mm/dia)
      - datas: lista de datetime com mesmo comprimento de time

    Saídas:
      - soma_mes: (12, lat, lon) -> soma de todos os dias daquele mês (mm)
      - media_mes: (12, lat, lon) -> média diária daquele mês (mm/dia)
      - cont_mes: (12,) -> número de dias usados em cada mês
    """
    pilha = np.asarray(pilha)
    nt, nlat, nlon = pilha.shape

    soma_mes = np.zeros((12, nlat, nlon), dtype=np.float32)
    cont_mes = np.zeros(12, dtype=np.int32)

    for i, dt in enumerate(datas):
        m_idx = dt.month - 1  # 0..11
        soma_mes[m_idx] += pilha[i]
        cont_mes[m_idx] += 1

    media_mes = np.full_like(soma_mes, np.nan, dtype=np.float32)
    for m in range(12):
        if cont_mes[m] > 0:
            media_mes[m] = soma_mes[m] / cont_mes[m]

    return soma_mes, media_mes, cont_mes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gera NetCDF com soma/média da chuva e mapas.")
    parser.add_argument("data_inicial", help="Data inicial (formato flexível)")
    parser.add_argument("data_final", help="Data final (formato flexível)")
    parser.add_argument("--regiao", default="AMS", help="Região definida no cpc.config")
    parser.add_argument("--saida", default="chuva_acumulada.nc", help="Arquivo NetCDF de saída")
    parser.add_argument("--interpolacao", action="store_true", help="Aplica interpolação para suavizar os mapas")
    parser.add_argument("--resolucao", type=float, default=0.1, help="Resolução da interpolação (default: 0.1)")
    parser.add_argument("--shapefile", help="Arquivo shapefile a ser plotado no mapa")
    parser.add_argument(
        "--cores",
        nargs="+",
        help="Escala de cores ou chave do cpc.config (ex: --cores 0 10 20 ou --cores escala_verao)"
    )
    parser.add_argument("--figura", default="", help="Prefixo para nome dos arquivos gerados")
    parser.add_argument("--titulo", default="", help="Título exibido na figura antes do intervalo de datas")
    parser.add_argument("--colormap", default="turbo", help="Nome do colormap do Matplotlib (padrão: turbo)")
    parser.add_argument(
        "--politica",
        nargs="+",
        default=["black", "0.5", "0.3"],
        help="Fronteiras políticas: 'off' para desligar ou 'cor largura_borda largura_estado' "
             "(ex: --politica black 0.5 0.3)"
    )
    parser.add_argument(
        "--unidade",
        choices=["mm/dia", "mm/mes", "mm/ano"],
        default="mm/dia",
        help="Unidade da chuva média no NetCDF (mm/dia, mm/mes ou mm/ano). Padrão: mm/dia."
    )
    parser.add_argument(
        "--cor-shapefile",
        default=None,
        help="Cor da borda do shapefile (ex: --cor-shapefile yellow). Se omitido, usa 'red'."
    )
    parser.add_argument(
        "--clima",
        action="store_true",
        help="Se presente, também calcula climatologia mensal (média diária e soma por mês) ao longo do período."
    )

    args = parser.parse_args()

    # Validação do colormap
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

    lat, lon, pilha, datas = carregar_dados_periodo(repositorio, prefixo, data_ini, data_fim)

    chuva_soma = np.sum(pilha, axis=0)          # mm no período
    chuva_media_dia = np.mean(pilha, axis=0)    # mm/dia

    # Ajuste da unidade da média conforme --unidade
    if args.unidade == "mm/dia":
        chuva_media = chuva_media_dia
        unidade_nc = "mm/day"
    elif args.unidade == "mm/mes":
        chuva_media = chuva_media_dia * 30.0
        unidade_nc = "mm/month"
    else:  # mm/ano
        chuva_media = chuva_media_dia * 365.0
        unidade_nc = "mm/year"

    # Climatologia mensal (opcional)
    soma_mes = media_mes = cont_mes = None
    if args.clima:
        print("[INFO] Calculando climatologia mensal (soma e média diária por mês)...")
        soma_mes, media_mes, cont_mes = calcular_climatologia_mensal(pilha, datas)

    salvar_netcdf_saida(
        lat,
        lon,
        chuva_soma,
        chuva_media,
        data_ini,
        data_fim,
        saida=args.saida,
        unidade_media=unidade_nc,
        soma_mes=soma_mes,
        media_mes=media_mes,
        cont_mes=cont_mes
    )

    data_ini_str = data_ini.strftime("%Y%m%d")
    data_fim_str = data_fim.strftime("%Y%m%d")
    nome_soma = f"{args.figura}_soma_{data_ini_str}_{data_fim_str}.png"
    nome_media = f"{args.figura}_media_{data_ini_str}_{data_fim_str}.png"

    if args.interpolacao:
        lat_i, lon_i, soma_i = interpolar_grade(lat, lon, chuva_soma, resolucao=args.resolucao)
        lat_j, lon_j, media_j = interpolar_grade(lat, lon, chuva_media, resolucao=args.resolucao)
        gerar_mapa(
            soma_i, lat_i, lon_i,
            f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d}",
            cores, nome_soma, limites,
            shapefile=args.shapefile,
            colormap=args.colormap,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile
        )
        gerar_mapa(
            media_j, lat_j, lon_j,
            f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d} (média)",
            cores, nome_media, limites,
            shapefile=args.shapefile,
            colormap=args.colormap,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile
        )
    else:
        gerar_mapa(
            chuva_soma, lat, lon,
            f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d}",
            cores, nome_soma, limites,
            shapefile=args.shapefile,
            colormap=args.colormap,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile
        )
        gerar_mapa(
            chuva_media, lat, lon,
            f"{args.titulo} {data_ini:%Y-%m-%d} a {data_fim:%Y-%m-%d} (média)",
            cores, nome_media, limites,
            shapefile=args.shapefile,
            colormap=args.colormap,
            politica=args.politica,
            cor_shapefile=args.cor_shapefile
        )
