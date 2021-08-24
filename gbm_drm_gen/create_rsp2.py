import numpy as np

from .drmgen import DRMGen
from responsum.response import RSP2
from responsum.utils.time_interval import TimeInterval


def create_rsp2(file_name: str,
                response_generator: DRMGen,
                ra: float,
                dec: float,
                tstart: float,
                tstop: float,
                delta_time: float = 3,
                overwrite: bool=False) -> None:
    
    """

    """

    # convert to MET

    met_tstart = response_generator.met_at(tstart)

    met_tstop = response_generator.met_at(tstop)
    
    
    met_time_bins = np.arange(met_tstart, met_tstop, delta_time).tolist()

    met_time_bins.append(met_time_bins[-1] + delta_time)

    met_time_bins = np.array(met_time_bins)
    
    met_start = met_time_bins[:-1]

    met_stop = met_time_bins[1:]


    time_bins = np.arange(tstart, tstop, delta_time).tolist()

    time_bins.append(time_bins[-1] + delta_time)

    time_bins = np.array(time_bins)
    
    start = time_bins[:-1]

    stop = time_bins[1:]

    
    list_of_matrices = []
    
    for a, b in zip(start, stop):

        coverage_interval = TimeInterval(a,b)

        mean_time = 0.5 * (a + b)

        response_generator.set_time(mean_time)

        rsp = response_generator.to_3ML_response(ra, dec, coverage_interval=coverage_interval)

        list_of_matrices.append(rsp.matrix)
        
    rsp2 = RSP2(rsp.monte_carlo_energies, rsp.ebounds, list_of_matrices, "Fermi", "GBM", met_start, met_stop)


    rsp2.writeto(file_name, overwrite=overwrite)
