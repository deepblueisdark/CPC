"""
cpc_plota_clima_unificado.py

Versão unificada que combina:
  - Da V1: --auto-escala, --painel (3x4 meses), --saida-dir, fallback de variável
            (chuva_soma/chuva_media -> "chuva"), validação de variáveis ausentes
            (KeyError com mensagem clara), leitura de coordenada "month" do NetCDF,
            --invert-cmap, tratamento de erros (try/except) em shapefile e colormap.
  - Da V2: interpolação de grade (--interpolacao/--resolucao), controle de
            fronteiras/estados via --politica (cor, espessura, ou "off"),
            --cor-shapefile customizável, nomes de meses por extenso em maiúsculo.

Uso básico (igual às versões anteriores):
    python cpc_plota_clima_unificado.py arquivo.nc --tipo ambos --clima --painel
"""

import argparse
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geopandas as gpd
import configparser
from scipy.interpolate import griddata

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


# ==========================================================
# INTERPOLAÇÃO (V2)
# ==========================================================
def interpolar_grade(lat, lon, campo, resolucao):
    """Interpola um campo 2D (lat, lon) para uma grade mais fina."""
    lon2d, lat2d = np.meshgrid(lon, lat)

    pontos = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    valores = np.asarray(campo).ravel()

    mask = np.isfinite(valores)
    pontos = pontos[mask]
    valores = valores[mask]

    novo_lon = np.arange(lon.min(), lon.max(), resolucao)
    novo_lat = np.arange(lat.min(), lat.max(), resolucao)

    lon_grid, lat_grid = np.meshgrid(novo_lon, novo_lat)
    campo_interp = griddata(pontos, valores, (lon_grid, lat_grid), method="linear")

    return lat_grid, lon_grid, campo_interp


# ==========================================================
# UNIDADE (V1/V2 - idênticas, mantendo os avisos da V1)
# ==========================================================
def fator_conversao_unidade(unidade_origem, unidade_destino_cli):
    """
    Converte de unidade_origem (mm/day, mm/month, mm/year)
    para unidade_destino_cli (mm/dia, mm/mes, mm/ano) via mm/dia como base.
    """
    if unidade_origem is None:
        print("[AVISO] Unidade de origem da chuva_media desconhecida; não será aplicada conversão.")
        return 1.0

    uo = unidade_origem.strip().lower()

    if uo == "mm/day":
        f_origem_para_dia = 1.0
    elif uo == "mm/month":
        f_origem_para_dia = 1.0 / 30.0
    elif uo == "mm/year":
        f_origem_para_dia = 1.0 / 365.0
    else:
        print(f"[AVISO] Unidade de origem '{unidade_origem}' não reconhecida; sem conversão.")
        return 1.0

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


# ==========================================================
# HELPERS DE ESTILO (colormap, política, shapefile)
# ==========================================================
def resolver_colormap(colormap, invert_cmap):
    try:
        cmap = plt.get_cmap(colormap)
    except ValueError:
        print(f"[AVISO] Colormap '{colormap}' não encontrado. Usando 'turbo'.")
        cmap = plt.get_cmap("turbo")
    if invert_cmap:
        cmap = cmap.reversed()
    return cmap


def aplicar_politica(ax, politica):
    """Controla fronteiras/estados: None=padrão, 'off'=desliga, ou [cor, lw_borda, lw_estado]."""
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


def carregar_shapefile(shapefile):
    if not (shapefile and os.path.exists(shapefile)):
        return None
    try:
        gdf = gpd.read_file(shapefile)
        if gdf.crs is None:
            print("[AVISO] Shapefile sem CRS. Assumindo EPSG:4326.")
            gdf.set_crs("EPSG:4326", inplace=True)
        return gdf.to_crs("EPSG:4326")
    except Exception as e:
        print(f"[ERRO] Falha ao ler shapefile: {shapefile}\n{e}")
        return None


# ==========================================================
# MAPA INDIVIDUAL
# ==========================================================
def gerar_mapa(chuva, lat, lon, titulo, cores, saida_png,
               shapefile=None, limites=None, label_unidade="mm",
               colormap="turbo", invert_cmap=False,
               interpolacao=False, resolucao=0.1,
               politica=None, cor_shapefile="black"):

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
    aplicar_politica(ax, politica)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    norm = BoundaryNorm(cores, ncolors=256)
    cmap = resolver_colormap(colormap, invert_cmap)

    cs = ax.pcolormesh(lon, lat, chuva, cmap=cmap, norm=norm, transform=proj)

    cb = plt.colorbar(cs, orientation="horizontal", pad=0.08, aspect=40, shrink=0.9)
    cb.set_label(label_unidade)
    cb.set_ticks(cores)
    labels = [f"{c:.0f}" if c >= 1 else f"{c:.1f}" for c in cores]
    cb.ax.set_xticklabels(labels)
    cb.ax.tick_params(labelsize=6)

    gdf = carregar_shapefile(shapefile)
    if gdf is not None:
        gdf.boundary.plot(ax=ax, edgecolor=cor_shapefile, linewidth=1, alpha=0.8)

    plt.title(titulo)
    plt.tight_layout()
    plt.savefig(saida_png, dpi=150, bbox_inches="tight")
    plt.close()


# ==========================================================
# PAINEL 3x4 (V1)
# ==========================================================
def gerar_painel_mensal(dados, lat, lon, nomes_meses, cores, titulo_base,
                         saida_png, shapefile=None, limites=None,
                         label_unidade="mm", colormap="turbo", invert_cmap=False,
                         interpolacao=False, resolucao=0.1,
                         politica=None, cor_shapefile="black"):

    proj = ccrs.PlateCarree()
    fig, axes = plt.subplots(3, 4, figsize=(14, 8), subplot_kw={"projection": proj})
    axes = axes.ravel()

    if limites:
        min_lat, max_lat, min_lon, max_lon = limites
    else:
        min_lat, max_lat = float(lat.min()), float(lat.max())
        min_lon, max_lon = float(lon.min()), float(lon.max())

    norm = BoundaryNorm(cores, ncolors=256)
    cmap = resolver_colormap(colormap, invert_cmap)
    gdf = carregar_shapefile(shapefile)

    nmes = dados.shape[0]
    cs = None

    for i in range(12):
        ax = axes[i]
        if i < nmes:
            campo = dados[i, :, :]
            plat, plon = lat, lon

            if interpolacao:
                plat, plon, campo = interpolar_grade(lat, lon, campo, resolucao)

            ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=proj)
            ax.coastlines(linewidth=0.4)
            aplicar_politica(ax, politica)

            cs = ax.pcolormesh(plon, plat, campo, cmap=cmap, norm=norm, transform=proj)

            mes_num = i + 1
            nome_mes = nomes_meses.get(mes_num, f"M{mes_num:02d}")
            ax.set_title(nome_mes, fontsize=9, pad=2)

            if gdf is not None:
                gdf.boundary.plot(ax=ax, edgecolor=cor_shapefile, linewidth=0.6, alpha=0.7)

            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_visible(False)

    plt.subplots_adjust(wspace=0.05, hspace=0.08, top=0.93, bottom=0.16)

    if cs is not None:
        cb = fig.colorbar(cs, ax=axes, orientation="horizontal", fraction=0.06, pad=0.08, aspect=40)
        cb.set_label(label_unidade, fontsize=9)
        cb.set_ticks(cores)
        labels = [f"{c:.0f}" if c >= 1 else f"{c:.1f}" for c in cores]
        cb.ax.set_xticklabels(labels, fontsize=8)

    fig.suptitle(titulo_base, fontsize=13, y=0.98)
    plt.savefig(saida_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ==========================================================
# MAIN
# ==========================================================
def main():
    parser = argparse.ArgumentParser(
        description="Gera mapas de soma e/ou média a partir de um NetCDF do CPC "
                     "(inclui climatologia mensal, painel 3x4, interpolação, "
                     "controle de fronteiras/shapefile e escala automática)."
    )
    parser.add_argument("arquivo", help="Arquivo NetCDF de entrada")
    parser.add_argument("--tipo", choices=["soma", "media", "ambos"], default="soma",
                         help="Tipo de mapa a gerar")
    parser.add_argument("--regiao", help="Nome da região no cpc.config para limitar o mapa")
    parser.add_argument("--saida", default=None, help="Prefixo do nome dos arquivos PNG de saída")
    parser.add_argument("--saida-dir", default=".",
                         help="Diretório onde os PNG serão salvos (será criado se não existir).")
    parser.add_argument("--titulo", default="Mapa CPC", help="Título base a ser exibido nos mapas")

    parser.add_argument("--shapefile", help="Shapefile opcional para sobrepor no mapa")
    parser.add_argument("--cor-shapefile", default="black", help="Cor da linha do shapefile.")

    parser.add_argument("--cores", nargs="+", help="Escala de cores ou chave definida no cpc.config")
    parser.add_argument("--colormap", "--cmap", dest="colormap", default="turbo",
                         help="Nome do colormap do Matplotlib (ex: turbo, viridis, terrain, Spectral_r).")
    parser.add_argument("--invert-cmap", action="store_true",
                         help="Inverte a escala de cores.")
    parser.add_argument("--auto-escala", action="store_true",
                         help="Gera automaticamente a escala de cores a partir dos dados "
                              "(ignora --cores e cpc.config).")

    parser.add_argument("--interpolacao", action="store_true",
                         help="Interpola os dados para uma grade mais fina antes de plotar.")
    parser.add_argument("--resolucao", type=float, default=0.1,
                         help="Resolução (graus) da grade interpolada, se --interpolacao ativo.")

    parser.add_argument("--politica", nargs="+",
                         help="Controla fronteiras/estados: 'off' desliga, ou "
                              "'<cor> [lw_borda] [lw_estado]' para customizar.")

    parser.add_argument("--clima", action="store_true",
                         help="Usa climatologia mensal: lê chuva_soma_mes/chuva_media_mes.")
    parser.add_argument("--painel", action="store_true",
                         help="Gera também um painel 3x4 com os 12 meses (requer --clima).")

    parser.add_argument("--unidade", choices=["mm/dia", "mm/mes", "mm/ano"], default=None,
                         help="Reescala a chuva_media para mm/dia, mm/mes ou mm/ano.")

    args = parser.parse_args()

    saida_dir = args.saida_dir
    os.makedirs(saida_dir, exist_ok=True)

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    # ===== REGIÃO =====
    limites = None

    with nc.Dataset(args.arquivo) as ds:
        lat = ds.variables["lat"][:]
        lon = ds.variables["lon"][:]

        if args.regiao and args.regiao in config:
            reg = config[args.regiao]
            limites = (
                float(reg.get("min_lat", lat.min())),
                float(reg.get("max_lat", lat.max())),
                float(reg.get("min_lon", lon.min())),
                float(reg.get("max_lon", lon.max())),
            )

        tipos = ["soma", "media"] if args.tipo == "ambos" else [args.tipo]
        mapas = {}
        unidade_media_origem = None
        mes_coord = None

        if args.clima:
            for tipo in tipos:
                varname = "chuva_soma_mes" if tipo == "soma" else "chuva_media_mes"
                if varname not in ds.variables:
                    raise KeyError(
                        f"Variável '{varname}' não encontrada no arquivo. "
                        f"Gere o NetCDF com --clima no script de soma/média ou remova --clima aqui."
                    )
                var = ds.variables[varname]
                mapas[tipo] = var[:]
                if tipo == "media":
                    unidade_media_origem = getattr(var, "units", None)

            if "month" in ds.variables:
                mes_coord = ds.variables["month"][:]
            else:
                nmes = mapas[tipos[0]].shape[0]
                mes_coord = np.arange(1, nmes + 1)
        else:
            for tipo in tipos:
                varname = f"chuva_{tipo}"
                if varname not in ds.variables:
                    varname = "chuva"  # fallback genérico
                if varname not in ds.variables:
                    raise KeyError(f"Variável '{varname}' não encontrada no arquivo.")
                var = ds.variables[varname]
                mapas[tipo] = var[:]
                if tipo == "media":
                    unidade_media_origem = getattr(var, "units", None)

        if "media" in mapas and args.unidade:
            if unidade_media_origem is None:
                print("[AVISO] Unidade de 'chuva_media' não encontrada; ignorando --unidade.")
            else:
                fator = fator_conversao_unidade(unidade_media_origem, args.unidade)
                mapas["media"] = mapas["media"] * fator

    # ===== CORES BASE =====
    cores_base = None
    if not args.auto_escala:
        if args.cores:
            if all(c.replace(".", "", 1).lstrip("-").isdigit() for c in args.cores):
                cores_base = list(map(float, args.cores))
            else:
                chave = args.cores[0]
                if chave in config["DEFAULT"]:
                    cores_base = list(map(float, config["DEFAULT"][chave].split()))
                else:
                    raise ValueError(f"Chave de cores '{chave}' não encontrada no cpc.config")
        else:
            cores_base = list(map(float, config["DEFAULT"].get(
                "escala_cores", "0.1 1 2 5 10 20 50 100"
            ).split()))
    elif args.cores:
        print("[AVISO] --auto-escala ativo: ignorando valores de --cores.")

    # ===== NOMES DOS MESES (por extenso, maiúsculo - V2) =====
    nomes_meses = {
        1: "JANEIRO", 2: "FEVEREIRO", 3: "MARÇO", 4: "ABRIL",
        5: "MAIO", 6: "JUNHO", 7: "JULHO", 8: "AGOSTO",
        9: "SETEMBRO", 10: "OUTUBRO", 11: "NOVEMBRO", 12: "DEZEMBRO",
    }

    prefixo = args.saida or "mapa"

    for tipo in tipos:
        dados = mapas[tipo]

        # ----- escala de cores por tipo -----
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
                cores_tipo = cores_base if cores_base is not None else [0, 1, 2, 5, 10, 20, 50, 100]
            else:
                vmin = float(np.nanpercentile(arr, 2))
                vmax = float(np.nanpercentile(arr, 98))
                if vmin == vmax:
                    vmin -= 0.5
                    vmax += 0.5
                cores_tipo = np.linspace(vmin, vmax, 13).tolist()
        else:
            cores_tipo = cores_base

        # ----- rótulo de unidade -----
        if tipo == "soma":
            label_unidade = "mm"
        else:
            label_unidade = args.unidade if args.unidade else "mm"

        # ----- 2D (modo normal) ou 3D (climatologia) -----
        if not args.clima or dados.ndim == 2:
            titulo = f"{args.titulo} ({tipo})"
            fname = args.saida or f"mapa_{tipo}.png"
            if args.tipo == "ambos" and args.saida:
                fname = f"{args.saida}_{tipo}.png"
            saida = os.path.join(saida_dir, fname)

            gerar_mapa(
                dados, lat, lon, titulo, cores_tipo, saida,
                shapefile=args.shapefile, limites=limites,
                label_unidade=label_unidade, colormap=args.colormap,
                invert_cmap=args.invert_cmap,
                interpolacao=args.interpolacao, resolucao=args.resolucao,
                politica=args.politica, cor_shapefile=args.cor_shapefile,
            )
        else:
            nmes = dados.shape[0]
            if mes_coord is not None and len(mes_coord) == nmes:
                meses_num = [int(m) for m in mes_coord]
            else:
                meses_num = list(range(1, nmes + 1))

            for idx in range(nmes):
                mes_num = meses_num[idx]
                nome_mes = nomes_meses.get(mes_num, f"M{mes_num:02d}")
                titulo = f"{nome_mes} - {args.titulo} ({tipo})"
                fname = f"{prefixo}_{tipo}_m{mes_num:02d}.png"
                saida = os.path.join(saida_dir, fname)

                gerar_mapa(
                    dados[idx], lat, lon, titulo, cores_tipo, saida,
                    shapefile=args.shapefile, limites=limites,
                    label_unidade=label_unidade, colormap=args.colormap,
                    invert_cmap=args.invert_cmap,
                    interpolacao=args.interpolacao, resolucao=args.resolucao,
                    politica=args.politica, cor_shapefile=args.cor_shapefile,
                )

            if args.painel:
                painel_titulo = f"{args.titulo} ({tipo}) - Climatologia mensal"
                painel_nome = os.path.join(saida_dir, f"{prefixo}_{tipo}_painel.png")
                gerar_painel_mensal(
                    dados, lat, lon, nomes_meses, cores_tipo,
                    painel_titulo, painel_nome,
                    shapefile=args.shapefile, limites=limites,
                    label_unidade=label_unidade, colormap=args.colormap,
                    invert_cmap=args.invert_cmap,
                    interpolacao=args.interpolacao, resolucao=args.resolucao,
                    politica=args.politica, cor_shapefile=args.cor_shapefile,
                )


if __name__ == "__main__":
    main()
