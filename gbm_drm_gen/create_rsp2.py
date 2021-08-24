from drmgen import DRMGen
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
    
    
    time_bins = np.arange(tstart, tstop, delta_time)

    start = time_bins[:-1]

    stop = time_bins[1:]

    list_of_matrices = []
    
    for a, b in zip(start, stop):

        coverage_interval = TimeInterval(a,b)

        mean_time = 0.5 * (a + b)

        response_generator.set_time(mean_time)

        rsp = response_generator.to_3ML_response(ra, dec, coverage_interval=coverage_interval)

        list_of_matrices.append(rsp)
        
    rsp2 = RSP2(rsp.monte_carlo_energies, rsp.ebounds, list_of_matrices, "Fermi", "GBM", start, stop)


    rsp2.to_fits(file_name, overwrite=overwrite)
