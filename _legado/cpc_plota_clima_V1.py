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

def gerar_mapa(chuva, lat, lon, titulo, cores, saida_png,
               shapefile=None, limites=None, label_unidade="mm",
               colormap="turbo", invert_cmap=False):
    """
    Gera um único mapa (2D) de chuva.
    """
    plt.figure(figsize=(10, 6))
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)
    else:
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=proj)

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.3)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    norm = BoundaryNorm(cores, ncolors=256)

    # colormap com opção de inversão
    try:
        cmap = plt.get_cmap(colormap)
    except ValueError:
        print(f"[AVISO] Colormap '{colormap}' não encontrado. Usando 'turbo'.")
        cmap = plt.get_cmap("turbo")
    if invert_cmap:
        cmap = cmap.reversed()

    cs = ax.pcolormesh(lon, lat, chuva, cmap=cmap, norm=norm, transform=proj)

    cb = plt.colorbar(cs, orientation='horizontal', pad=0.05, aspect=40)
    cb.set_label(label_unidade)
    cb.set_ticks(cores)
    tick_labels = [f"{c:.1f}" if c < 1 else f"{c:.0f}" for c in cores]
    cb.ax.set_xticklabels(tick_labels, fontsize=8)

    if shapefile and os.path.exists(shapefile):
        try:
            gdf = gpd.read_file(shapefile)
            if gdf.crs is None:
                gdf.set_crs("EPSG:4326", inplace=True)
            gdf.to_crs("EPSG:4326").boundary.plot(ax=ax, edgecolor='black', linewidth=1, alpha=0.8)
        except Exception as e:
            print(f"[ERRO] Falha ao ler shapefile: {shapefile}\n{e}")

    plt.title(titulo)
    plt.savefig(saida_png, dpi=150, bbox_inches="tight")
    plt.close()

def gerar_painel_mensal(dados, lat, lon, nomes_meses, cores, titulo_base,
                        saida_png, shapefile=None, limites=None,
                        label_unidade="mm", colormap="turbo", invert_cmap=False):
    """
    Gera um painel 3x4 com os 12 meses (ou até nmes) de 'dados' (month, lat, lon).
    """
    proj = ccrs.PlateCarree()
    fig, axes = plt.subplots(3, 4, figsize=(14, 8), subplot_kw={"projection": proj})
    axes = axes.ravel()

    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
    else:
        min_lat, max_lat = float(lat.min()), float(lat.max())
        min_lon, max_lon = float(lon.min()), float(lon.max())

    norm = BoundaryNorm(cores, ncolors=256)
    try:
        cmap = plt.get_cmap(colormap)
    except ValueError:
        print(f"[AVISO] Colormap '{colormap}' não encontrado. Usando 'turbo'.")
        cmap = plt.get_cmap("turbo")
    if invert_cmap:
        cmap = cmap.reversed()

    # Carrega shapefile uma vez só
    gdf = None
    if shapefile and os.path.exists(shapefile):
        try:
            gdf = gpd.read_file(shapefile)
            if gdf.crs is None:
                gdf.set_crs("EPSG:4326", inplace=True)
            gdf = gdf.to_crs("EPSG:4326")
        except Exception as e:
            print(f"[ERRO] Falha ao ler shapefile para painel: {shapefile}\n{e}")
            gdf = None

    nmes = dados.shape[0]
    cs = None

    for i in range(12):
        ax = axes[i]
        if i < nmes:
            campo = dados[i, :, :]
            ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)
            ax.coastlines(linewidth=0.4)
            ax.add_feature(cfeature.BORDERS, linewidth=0.3)
            ax.add_feature(cfeature.STATES, linewidth=0.2)

            cs = ax.pcolormesh(lon, lat, campo, cmap=cmap, norm=norm, transform=proj)

            mes_num = i + 1
            nome_mes = nomes_meses.get(mes_num, f"M{mes_num:02d}")
            ax.set_title(nome_mes, fontsize=9, pad=2)

            if gdf is not None:
                gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.6, alpha=0.7)

            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_visible(False)

    # Espaçamento entre subplots
    plt.subplots_adjust(wspace=0.05, hspace=0.08, top=0.93, bottom=0.16)

    if cs is not None:
        cb = fig.colorbar(cs, ax=axes, orientation="horizontal",
                          fraction=0.06, pad=0.08, aspect=40)
        cb.set_label(label_unidade, fontsize=9)

        # Todos os valores da escala na barra
        cb.set_ticks(cores)
        tick_labels = [f"{t:.1f}" if t < 1 else f"{t:.0f}" for t in cores]
        cb.ax.set_xticklabels(tick_labels, fontsize=8)

    fig.suptitle(titulo_base, fontsize=13, y=0.98)
    plt.savefig(saida_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

def fator_conversao_unidade(unidade_origem, unidade_destino_cli):
    """
    Converte de unidade_origem (mm/day, mm/month, mm/year)
    para unidade_destino_cli (mm/dia, mm/mes, mm/ano) via mm/dia como base.
    Retorna fator multiplicativo (novo = antigo * fator).
    """
    if unidade_origem is None:
        print("[AVISO] Unidade de origem da chuva_media desconhecida; não será aplicada conversão.")
        return 1.0

    uo = unidade_origem.strip().lower()

    # Fator para converter da unidade de origem para mm/dia
    if uo == "mm/day":
        f_origem_para_dia = 1.0
    elif uo == "mm/month":
        f_origem_para_dia = 1.0 / 30.0
    elif uo == "mm/year":
        f_origem_para_dia = 1.0 / 365.0
    else:
        print(f"[AVISO] Unidade de origem '{unidade_origem}' não reconhecida; sem conversão.")
        return 1.0

    # Fator para converter de mm/dia para unidade destino
    if unidade_destino_cli == "mm/dia":
        f_dia_para_dest = 1.0
    elif unidade_destino_cli == "mm/mes":
        f_dia_para_dest = 30.0
    elif unidade_destino_cli == "mm/ano":
        f_dia_para_dest = 365.0
    else:
        print(f"[AVISO] Unidade de destino '{unidade_destino_cli}' não reconhecida; sem conversão.")
        return 1.0

    return f_origem_para_dia * f_dia_para_dest

def main():
    parser = argparse.ArgumentParser(
        description="Gera mapas de soma e/ou média a partir de um NetCDF do CPC "
                    "(inclui climatologia mensal se --clima for usado)."
    )
    parser.add_argument("arquivo", help="Arquivo NetCDF de entrada")
    parser.add_argument("--tipo", choices=["soma", "media", "ambos"], default="soma",
                        help="Tipo de mapa a gerar")
    parser.add_argument("--cores", nargs='+', help="Escala de cores ou chave definida no cpc.config")
    parser.add_argument("--titulo", default="Mapa CPC", help="Título base a ser exibido nos mapas")
    parser.add_argument("--saida", default=None, help="Prefixo do nome dos arquivos PNG de saída")
    parser.add_argument("--saida-dir", default=".",
                        help="Diretório onde os PNG serão salvos (será criado se não existir).")
    parser.add_argument("--shapefile", help="Shapefile opcional para sobrepor no mapa")
    parser.add_argument("--regiao", help="Nome da região no cpc.config para limitar o mapa")
    parser.add_argument(
        "--unidade",
        choices=["mm/dia", "mm/mes", "mm/ano"],
        help="Reescala a chuva_media para mm/dia, mm/mes ou mm/ano "
             "(apenas para tipo 'media' ou 'ambos')."
    )
    parser.add_argument(
        "--auto-escala",
        action="store_true",
        help="Gera automaticamente a escala de cores a partir dos dados "
             "(ignora --cores e cpc.config)."
    )
    parser.add_argument(
        "--clima",
        action="store_true",
        help="Usa climatologia mensal: lê chuva_soma_mes/chuva_media_mes e gera mapas para todos os meses."
    )
    parser.add_argument(
        "--painel",
        action="store_true",
        help="Gera também um painel 3x4 com os 12 meses para cada tipo selecionado."
    )
    parser.add_argument(
        "--invert-cmap",
        action="store_true",
        help="Inverte a escala de cores (ex: azul para valores altos e vermelho para baixos)."
    )
    parser.add_argument(
        "--cmap",
        default="turbo",
        help="Nome do colormap do Matplotlib (ex: turbo, viridis, terrain, Spectral_r)."
    )

    args = parser.parse_args()

    saida_dir = args.saida_dir
    os.makedirs(saida_dir, exist_ok=True)

    config = configparser.ConfigParser()
    config.read("cpc.config")

    with nc.Dataset(args.arquivo) as ds:
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]

        tipos = ["soma", "media"] if args.tipo == "ambos" else [args.tipo]
        mapas = {}
        unidade_media_origem = None
        mes_coord = None

        if args.clima:
            # Climatologia mensal: esperamos variáveis (month, lat, lon)
            for tipo in tipos:
                if tipo == "soma":
                    varname = "chuva_soma_mes"
                else:  # "media"
                    varname = "chuva_media_mes"

                if varname not in ds.variables:
                    raise KeyError(
                        f"Variável '{varname}' não encontrada no arquivo. "
                        f"Gere o NetCDF com --clima no script de soma/média ou remova --clima aqui."
                    )

                var = ds.variables[varname]
                mapas[tipo] = var[:]  # (month, lat, lon)

                if tipo == "media":
                    unidade_media_origem = getattr(var, "units", None)

            # Coordenada de mês (opcional, senão assume 1..N)
            if "month" in ds.variables:
                mes_coord = ds.variables["month"][:]
            else:
                algum_tipo = tipos[0]
                nmes = mapas[algum_tipo].shape[0]
                mes_coord = np.arange(1, nmes + 1)

        else:
            # Modo "normal" 2D: chuva_soma / chuva_media
            for tipo in tipos:
                varname = f"chuva_{tipo}"
                if varname not in ds.variables:
                    varname = "chuva"  # fallback genérico
                if varname not in ds.variables:
                    raise KeyError(f"Variável '{varname}' não encontrada no arquivo.")
                var = ds.variables[varname]
                mapas[tipo] = var[:]  # 2D

                if tipo == "media" and varname == "chuva_media":
                    unidade_media_origem = getattr(var, "units", None)

        # Ajuste de unidade para média (se solicitado)
        if "media" in mapas and args.unidade:
            if unidade_media_origem is None:
                print("[AVISO] Unidade de 'chuva_media' não encontrada; ignorando --unidade.")
            else:
                fator = fator_conversao_unidade(unidade_media_origem, args.unidade)
                mapas["media"] = mapas["media"] * fator

    # Definição da escala de cores base
    cores_base = None
    if not args.auto_escala:
        if args.cores:
            if all(c.replace('.', '', 1).isdigit() for c in args.cores):
                cores_base = list(map(float, args.cores))
            else:
                chave = args.cores[0]
                if chave in config['DEFAULT']:
                    cores_base = list(map(float, config['DEFAULT'][chave].split()))
                else:
                    raise ValueError(f"Chave de cores '{chave}' não encontrada no cpc.config")
        else:
            cores_base = list(map(float, config['DEFAULT'].get(
                "escala_cores", "0.1 1 2 3 5 10 20 50 100"
            ).split()))
    else:
        if args.cores:
            print("[AVISO] --auto-escala ativo: ignorando valores de --cores.")
        # cores_base será definida por tipo mais abaixo

    # Limites da região (se definidos)
    limites = None
    if args.regiao and args.regiao in config:
        reg = config[args.regiao]
        limites = (
            float(reg.get("min_lat", lat.min())),
            float(reg.get("max_lat", lat.max())),
            float(reg.get("min_lon", lon.min())),
            float(reg.get("max_lon", lon.max())),
        )

    # Mês -> nome (para título)
    nomes_meses = {
        1: "Jan", 2: "Fev", 3: "Mar", 4: "Abr",
        5: "Mai", 6: "Jun", 7: "Jul", 8: "Ago",
        9: "Set", 10: "Out", 11: "Nov", 12: "Dez"
    }

    prefixo = args.saida or "mapa"

    for tipo in tipos:
        dados = mapas[tipo]

        # Define escala de cores por tipo
        if args.auto_escala:
            arr = dados
            if isinstance(arr, np.ma.MaskedArray):
                arr = arr.filled(np.nan)
            arr = np.array(arr, dtype=float)
            arr[~np.isfinite(arr)] = np.nan
            arr[arr < 0] = np.nan

            if np.all(np.isnan(arr)):
                print(f"[AVISO] Todos os dados do tipo '{tipo}' são NaN após filtragem; "
                      "usando escala padrão do cpc.config.")
                if cores_base is None:
                    cores_tipo = [0, 1, 2, 5, 10, 20, 50, 100]
                else:
                    cores_tipo = cores_base
            else:
                vmin = float(np.nanpercentile(arr, 2))
                vmax = float(np.nanpercentile(arr, 98))
                if vmin == vmax:
                    vmin -= 0.5
                    vmax += 0.5
                niveis = np.linspace(vmin, vmax, 13)
                cores_tipo = niveis.tolist()
        else:
            cores_tipo = cores_base

        # Rótulo de unidade
        if tipo == "soma":
            label_unidade = "mm"
        else:  # tipo == "media"
            label_unidade = args.unidade if args.unidade else "mm"

        # 2D (modo normal) ou 3D (clima mensal)
        if not args.clima or dados.ndim == 2:
            chuva = dados
            titulo = f"{args.titulo} ({tipo})"
            saida = args.saida or f"mapa_{tipo}.png"
            if args.tipo == "ambos" and args.saida:
                saida = f"{args.saida}_{tipo}.png"
            saida = os.path.join(saida_dir, saida)

            gerar_mapa(
                chuva, lat, lon, titulo, cores_tipo, saida,
                shapefile=args.shapefile,
                limites=limites,
                label_unidade=label_unidade,
                colormap=args.cmap,
                invert_cmap=args.invert_cmap
            )
        else:
            # Climatologia mensal: dados (month, lat, lon)
            nmes = dados.shape[0]
            if mes_coord is not None and len(mes_coord) == nmes:
                meses_num = [int(m) for m in mes_coord]
            else:
                meses_num = list(range(1, nmes + 1))

            # Mapas individuais por mês
            for idx in range(nmes):
                chuva = dados[idx, :, :]
                mes_num = meses_num[idx]
                nome_mes = nomes_meses.get(mes_num, f"M{mes_num:02d}")
                titulo = f"{args.titulo} ({tipo}) - {nome_mes}"

                fname = f"{prefixo}_{tipo}_m{mes_num:02d}.png"
                saida = os.path.join(saida_dir, fname)

                gerar_mapa(
                    chuva, lat, lon, titulo, cores_tipo, saida,
                    shapefile=args.shapefile,
                    limites=limites,
                    label_unidade=label_unidade,
                    colormap=args.cmap,
                    invert_cmap=args.invert_cmap
                )

            # Painel 3x4 com os 12 meses
            if args.painel:
                painel_titulo = f"{args.titulo} ({tipo}) - Climatologia mensal"
                painel_nome = os.path.join(saida_dir, f"{prefixo}_{tipo}_painel.png")
                gerar_painel_mensal(
                    dados, lat, lon, nomes_meses, cores_tipo,
                    painel_titulo, painel_nome,
                    shapefile=args.shapefile,
                    limites=limites,
                    label_unidade=label_unidade,
                    colormap=args.cmap,
                    invert_cmap=args.invert_cmap
                )

if __name__ == "__main__":
    main()
