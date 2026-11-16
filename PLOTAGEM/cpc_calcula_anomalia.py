import argparse
import calendar

import numpy as np
import xarray as xr


def calcular_anomalia(
    tempo_nc, clima_nc, saida_nc, limiar_min, limite_percentual, metodo_interpolacao
):
    """Calcula anomalias absoluta e percentual de chuva mensal."""
    with xr.open_dataset(tempo_nc) as ds_obs, xr.open_dataset(clima_nc) as ds_clim:
        data_inicial = ds_obs.attrs["data_inicial"]
        data_final = ds_obs.attrs["data_final"]
        ano, mes, _ = map(int, data_inicial.split("-")[:3])

        print(f"Mês identificado: {mes}")

        clim_media = ds_clim["chuva_media_mes"].sel(month=mes)
        dias_mes = calendar.monthrange(ano, mes)[1]
        clim_soma_equivalente = clim_media * dias_mes

        obs_soma = ds_obs["chuva_soma"]

        # Os arquivos observados e de climatologia podem usar grades distintas.
        # A climatologia precisa estar nas mesmas coordenadas da observação antes
        # do cálculo; usar join='exact' interrompia a execução nesses casos.
        if not (
            obs_soma["lat"].equals(clim_soma_equivalente["lat"])
            and obs_soma["lon"].equals(clim_soma_equivalente["lon"])
        ):
            print(
                "[INFO] Reamostrando a climatologia para a grade observada "
                f"(método: {metodo_interpolacao})."
            )
            clim_soma_equivalente = clim_soma_equivalente.interp(
                lat=obs_soma["lat"],
                lon=obs_soma["lon"],
                method=metodo_interpolacao,
            )

        anomalia = obs_soma - clim_soma_equivalente
        anomalia_perc = xr.where(
            clim_soma_equivalente > limiar_min,
            (obs_soma / clim_soma_equivalente - 1) * 100,
            np.nan,
        ).clip(min=-limite_percentual, max=limite_percentual)

        ds_out = xr.Dataset(
            {
                "chuva_soma_obs": obs_soma,
                "chuva_climatologia": clim_soma_equivalente,
                "anomalia": anomalia,
                "anomalia_percentual": anomalia_perc,
            }
        )

        ds_out.attrs.update(
            {
                "descricao": (
                    "Anomalia de precipitação CPC baseada em climatologia mensal "
                    "(média diária convertida para soma mensal)"
                ),
                "metodo": "obs - climatologia; percentual = (obs/clim - 1)*100",
                "mes_referencia": mes,
                "data_inicial": data_inicial,
                "data_final": data_final,
            }
        )

        ds_out["chuva_soma_obs"].attrs.update(
            units="mm", long_name="Precipitação observada no período"
        )
        ds_out["chuva_climatologia"].attrs.update(
            units="mm", long_name="Climatologia mensal equivalente"
        )
        ds_out["anomalia"].attrs.update(
            units="mm", long_name="Anomalia absoluta (obs - climatologia)"
        )
        ds_out["anomalia_percentual"].attrs.update(
            units="%", long_name="Anomalia percentual de precipitação"
        )

        ds_out.to_netcdf(saida_nc)

    print(f"Arquivo salvo em: {saida_nc}")


def main():
    parser = argparse.ArgumentParser(
        description="Calcula anomalias de precipitação CPC em relação à climatologia mensal."
    )
    parser.add_argument("tempo_nc", help="NetCDF de chuva observada acumulada")
    parser.add_argument("clima_nc", help="NetCDF da climatologia mensal")
    parser.add_argument(
        "--saida", default="anomalia.nc", help="Arquivo NetCDF de saída (padrão: anomalia.nc)"
    )
    parser.add_argument(
        "--limiar-min",
        type=float,
        default=1.0,
        help="Climatologia mínima (mm) para calcular o percentual (padrão: 1.0)",
    )
    parser.add_argument(
        "--limite-percentual",
        type=float,
        default=500.0,
        help="Valor absoluto máximo da anomalia percentual (padrão: 500)",
    )
    parser.add_argument(
        "--metodo-interpolacao",
        choices=["linear", "nearest"],
        default="linear",
        help="Método para adequar a climatologia à grade observada (padrão: linear)",
    )
    args = parser.parse_args()

    if args.limiar_min < 0:
        parser.error("--limiar-min deve ser maior ou igual a zero")
    if args.limite_percentual <= 0:
        parser.error("--limite-percentual deve ser maior que zero")

    calcular_anomalia(
        args.tempo_nc,
        args.clima_nc,
        args.saida,
        args.limiar_min,
        args.limite_percentual,
        args.metodo_interpolacao,
    )


if __name__ == "__main__":
    main()
