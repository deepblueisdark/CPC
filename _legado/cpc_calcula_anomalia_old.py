import xarray as xr
import numpy as np
import calendar 

# ==========================
# ARQUIVOS
# ==========================
TEMPO_NC = "JUL2026_CPC.nc"
CLIMA_NC = "CLIMA1990_2020.nc"
SAIDA_NC = "ANOMALIA_JUL2026.nc"

# ==========================
# ABRIR DADOS
# ==========================
ds_obs = xr.open_dataset(TEMPO_NC)
ds_clim = xr.open_dataset(CLIMA_NC)

# ==========================
# EXTRAIR MÊS
# ==========================
data_inicial = ds_obs.attrs["data_inicial"]
ano = int(data_inicial.split("-")[0])
mes = int(data_inicial.split("-")[1])

print(f"Mês identificado: {mes}")

# ==========================
# CLIMATOLOGIA
# ==========================
clim_media = ds_clim["chuva_media_mes"].sel(month=mes)

# dias reais do mês (CORRETO)
dias_mes = calendar.monthrange(ano, mes)[1]

# converter para mm/mês
clim_soma_equivalente = clim_media * dias_mes

# ==========================
# OBSERVADO
# ==========================
obs_soma = ds_obs["chuva_soma"]

# ==========================
# ALINHAR
# ==========================
clim_soma_equivalente = clim_soma_equivalente.interp(
    lat=obs_soma.lat,
    lon=obs_soma.lon,
    method="linear"
)

# ==========================
# ANOMALIA ABSOLUTA
# ==========================
anomalia = obs_soma - clim_soma_equivalente

# ==========================
# ANOMALIA PERCENTUAL (ROBUSTA)
# ==========================
# evita valores absurdos em regiões muito secas
limiar_min = 1.0  # mm

anomalia_perc = xr.where(
    clim_soma_equivalente > limiar_min,
    (obs_soma / clim_soma_equivalente - 1) * 100,
    np.nan
)

# opcional: limitar extremos (evita distorção visual)
anomalia_perc = anomalia_perc.clip(min=-500, max=500)

# ==========================
# DATASET DE SAÍDA
# ==========================
ds_out = xr.Dataset(
    {
        "chuva_soma_obs": obs_soma,
        "chuva_climatologia": clim_soma_equivalente,
        "anomalia": anomalia,
        "anomalia_percentual": anomalia_perc
    },
    coords={
        "lat": ds_obs["lat"],
        "lon": ds_obs["lon"]
    }
)

# ==========================
# ATRIBUTOS GLOBAIS
# ==========================
ds_out.attrs["descricao"] = (
    "Anomalia de precipitação CPC baseada em climatologia mensal "
    "(media diaria convertida para soma mensal)"
)
ds_out.attrs["metodo"] = "obs - climatologia; percentual = (obs/clim - 1)*100"
ds_out.attrs["mes_referencia"] = mes
ds_out.attrs["data_inicial"] = ds_obs.attrs["data_inicial"]
ds_out.attrs["data_final"] = ds_obs.attrs["data_final"]

# ==========================
# ATRIBUTOS VARIÁVEIS
# ==========================
ds_out["chuva_soma_obs"].attrs["units"] = "mm"
ds_out["chuva_soma_obs"].attrs["long_name"] = "Precipitação observada no período"

ds_out["chuva_climatologia"].attrs["units"] = "mm"
ds_out["chuva_climatologia"].attrs["long_name"] = "Climatologia mensal equivalente"

ds_out["anomalia"].attrs["units"] = "mm"
ds_out["anomalia"].attrs["long_name"] = "Anomalia absoluta (obs - climatologia)"

ds_out["anomalia_percentual"].attrs["units"] = "%"
ds_out["anomalia_percentual"].attrs["long_name"] = "Anomalia percentual de precipitação"

# ==========================
# SALVAR
# ==========================
ds_out.to_netcdf(SAIDA_NC)

print(f"Arquivo salvo em: {SAIDA_NC}")
