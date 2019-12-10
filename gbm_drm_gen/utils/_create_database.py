import h5py
import bitshuffle.h5
from gbm_drm_gen.detdatabase import DetDatabase


# This code creates an HDF5 library from the FITS library


f = h5py.File("balrog_db.h5", "w", libver="latest")

dets = ["n%d" % i for i in range(10)]
dets.extend(["na", "nb", "b0", "b1"])

for det in dets:

    dd = DetDatabase(det)

    grp = f.create_group(det)

    for k, v in dd._rsp_dict.iteritems():
        grp.create_dataset(k, data=v, compression="lzf")

    grp.create_dataset("at_scat_data", data=dd.at_scat_data, compression="lzf")

    grp.create_dataset("Azimuth", data=dd.Azimuth, compression="lzf")
    grp.create_dataset(
        "double_phi_cent",
        data=dd.double_phi_cent,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "e_in",
        data=dd.e_in,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "energ_hi",
        data=dd.energ_hi,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "energ_lo",
        data=dd.energ_lo,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "epx_hi",
        data=dd.epx_hi,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "epx_lo",
        data=dd.epx_lo,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "ichan",
        data=dd.ichan,
        # compression='lzf'
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "ienerg",
        data=dd.ienerg,
        # compression='lzf'
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "lat_cent",
        data=dd.lat_cent,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "lat_edge",
        data=dd.lat_edge,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset(
        "LEND",
        data=dd.LEND,
        compression="lzf"
        # compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4
    )
    grp.create_dataset("LIST", data=dd.LIST, compression="lzf")
    grp.create_dataset("LPTR", data=dd.LPTR, compression="lzf")
    grp.create_dataset("milliaz", data=dd.milliaz, compression="lzf")
    grp.create_dataset("millizen", data=dd.millizen, compression="lzf")
    grp.create_dataset("phi_cent", data=dd.phi_cent, compression="lzf")
    grp.create_dataset("phi_edge", data=dd.phi_edge, compression="lzf")
    grp.create_dataset("theta_cent", data=dd.theta_cent, compression="lzf")
    grp.create_dataset("theta_edge", data=dd.theta_edge, compression="lzf")
    grp.create_dataset("X", data=dd.X, compression="lzf")
    grp.create_dataset("Y", data=dd.Y, compression="lzf")
    grp.create_dataset("Z", data=dd.Z, compression="lzf")
    grp.create_dataset("Zenith", data=dd.Zenith, compression="lzf")

    f.close()
