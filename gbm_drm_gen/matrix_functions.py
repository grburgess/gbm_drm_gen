import numpy as np
import numba as nb


@nb.njit(fastmath=True)
def geocords(theta_geo, phi_geo, theta_source, phi_source):

    gx = np.empty(3)
    gy = np.empty(3)
    gz = np.empty(3)
    sl = np.empty(3)

    gz[0] = np.sin(theta_geo * dtr) * np.cos(phi_geo * dtr)
    gz[1] = np.sin(theta_geo * dtr) * np.sin(phi_geo * dtr)
    gz[2] = np.cos(theta_geo * dtr)

    gzr = np.sqrt(gz[0] * gz[0] + gz[1] * gz[1] + gz[2] * gz[2])
    gz[0] = gz[0] / gzr
    gz[1] = gz[1] / gzr
    gz[2] = gz[2] / gzr

    sl[0] = np.sin(theta_source * dtr) * np.cos(phi_source * dtr)
    sl[1] = np.sin(theta_source * dtr) * np.sin(phi_source * dtr)
    sl[2] = np.cos(theta_source * dtr)

    slr = np.sqrt(sl[0] * sl[0] + sl[1] * sl[1] + sl[2] * sl[2])
    sl[0] = sl[0] / slr
    sl[1] = sl[1] / slr
    sl[2] = sl[2] / slr

    gy[0] = gz[1] * sl[2] - gz[2] * sl[1]
    gy[1] = gz[2] * sl[0] - gz[0] * sl[2]
    gy[2] = gz[0] * sl[1] - gz[1] * sl[0]

    gyr = np.sqrt(gy[0] * gy[0] + gy[1] * gy[1] + gy[2] * gy[2])

    gy[0] = gy[0] / gyr
    gy[1] = gy[1] / gyr
    gy[2] = gy[2] / gyr

    gx[0] = gy[1] * gz[2] - gy[2] * gz[1]
    gx[1] = gy[2] * gz[0] - gy[0] * gz[2]
    gx[2] = gy[0] * gz[1] - gy[1] * gz[0]

    gxr = np.sqrt(gx[0] * gx[0] + gx[1] * gx[1] + gx[2] * gx[2])
    gx[0] = gx[0] / gxr
    gx[1] = gx[1] / gxr
    gx[2] = gx[2] / gxr

    return gx, gy, gz, sl


@nb.njit(fastmath=True)
def geo_to_space(theta_u, phi_u, gx, gy, gz):

    dtr = np.arccos(-1.0) / 180.0

    xg = np.sin(theta_u * dtr) * np.cos(phi_u * dtr)
    yg = np.sin(theta_u * dtr) * np.sin(phi_u * dtr)
    zg = np.cos(theta_u * dtr)

    dirx = xg * gx[0] + yg * gy[0] + zg * gz[0]
    diry = xg * gx[1] + yg * gy[1] + zg * gz[1]
    dirz = xg * gx[2] + yg * gy[2] + zg * gz[2]

    r = np.sqrt(dirx * dirx + diry * diry + dirz * dirz)
    dirx = dirx / r
    diry = diry / r
    dirz = dirz / r

    az = np.arctan2(diry, dirx) / dtr

    if az < 0.0:
        az = az + 360.0

    el = 90.0 - np.arccos(dirz) / dtr

    return dirx, diry, dirz, az, el


@nb.njit(fastmath=True)
def calc_sphere_dist(ra1, dec1, ra2, dec2):
    dtr = np.arccos(-1.0) / 180.0
    y = np.sqrt(
        (np.cos(dec2 * dtr) * np.sin((ra1 - ra2) * dtr)) ** 2
        + (
            np.cos(dec1 * dtr) * np.sin(dec2 * dtr)
            - np.sin(dec1 * dtr) * np.cos(dec2 * dtr) * np.cos((ra1 - ra2) * dtr)
        )
        ** 2
    )
    x = np.sin(dec1 * dtr) * np.sin(dec2 * dtr) + np.cos(dec1 * dtr) * npcos(
        dec2 * dtr
    ) * np.cos((ra1 - ra2) * dtr)
    dist = np.arctan2(y, x) / dtr


@nb.njit(fastmath=True)
def highres_ephoton_interpolator(
    ebin_edge_in, nobins_in, ein, nvbins, matrix, edif_edge_lo, edif_edge_hi, nhbins
):

    new_epx_lo = np.empty((nobins_in, 64))
    new_epx_hi = np.empty((nobins_in, 64))

    diff_matrix = np.empty((nobins_in, 64))

    ivfind = 1

    for i in range(nobins_in):
        for j in range(ivfind - 1, nvbins):
            if (ebin_edge_in[i] >= e_in[j]) and ebin_edge_in[i] < e_in[j + 1]:
                ivfind = j

                mu = (np.log(ebin_edge_in[i]) - np.log(e_in[ivfind])) / (
                    np.loglog(e_in[ivfind + 1]) - np.log(e_in[ivfind])
                )
                if (mu < 0.0) and (mu > -1e-5):
                    mu = 0.0
                elif (mu > 1.0) and (mu < 1.00001):
                    mu = 1.0

                for k in range(nhbins):

                    new_epx_lo[i, k] = (
                        edif_edge_lo[ivfind, k] / ein[ivfind] * (1 - mu)
                        + edif_edge_lo[ivfind + 1, k] / ein[ivfind + 1] * mu
                    ) * ebin_edge_in[i]
                    new_epx_hi[i, k] = (
                        edif_edge_hi[ivfind, k] / ein[ivfind] * (1 - mu)
                        + edif_edge_hi[ivfind + 1, k] / ein[ivfind + 1] * mu
                    ) * ebin_edge_in[i]
                    diff_matrix[i, k] = (
                        matrix[ivfind, k] * (1 - mu) + matrix[ivfind + 1, k] * mu
                    )
    return new_epx_lo, new_epx_hi, highres_ephoton_interpolator


@nb.njit(fastmath=True)
def atscat_highres_ephoton_interpolator(
    ebin_edge_in, nobins_in, e_in, nvbins, matrix, nobins_out
):

    new_matrix = np.zeros((nobins_out, nobins_in))

    ivfind = 1

    for i in range(nobins_in):
        for j in range(ivfind - 1, nvbins):
            if (ebin_edge_in[i] >= e_in[j]) and ebin_edge_in[i] < e_in[j + 1]:
                ivfind = j
                mu = (np.log(ebin_edge_in[i]) - np.log(e_in[ivfind])) / (
                    np.loglog(e_in[ivfind + 1]) - np.log(e_in[ivfind])
                )
                if (mu < 0.0) and (mu > -1e-5):
                    mu = 0.0
                elif (mu > 1.0) and (mu < 1.00001):
                    mu = 1.0

                for k in range(nobins_out):
                    new_matrix[i, k] = (
                        matrix[ivfind, k] * (1 - mu) + matrix[ivfind + 1, k] * mu
                    )

    return new_matrix


@nb.njit(fastmath=True)
def echan_integrator(diff_matrix, edif_edge_lo, edif_edge_hi, nhbins, ebin_edge_out):

    nobins_in = diff_matrix.shape[0]

    nobins_out = len(ebin_edge_out) - 1

    # check that this is the right order from FORTRAN
    binned_matrix = np.zeros((nobins_in, nobins_out))
    # binned_matrix = np.zeros((nobins_out,nobins_in))

    row_tot = np.zeros(nobins_out + 1)
    diff_matrix_vec = np.empty(nhbins)

    edif_edgeh = np.empty(nhbins + 1)
    edif_cent = np.empty(nhbins)
    # first is a loop over the photon energies
    # row_entry = 0.
    # ihover =0
    for jcdif in range(1, nobins_in + 1):

        for ivh in range(1, nhbins + 1):

            diff_matrix_vec[ivh - 1] = diff_matrix[jcdif - 1, ivh - 1] / (
                edif_edge_hi[jcdif - 1, ivh - 1] - edif_edge_lo[jcdif - 1, ivh - 1]
            )
            edif_edgeh[ivh - 1] = edif_edge_hi[jcdif - 1, ivh - 1]
            edif_cent[ivh - 1] = (
                edif_edge_lo[jcdif - 1, ivh - 1] + edif_edge_hi[jcdif - 1, ivh - 1]
            ) / 2.0

        edif_edgeh[nhbins] = edif_edge_hi[jcdif - 1, nhbins - 1] + (
            edif_edge_hi[jcdif - 1, nhbins - 1,] - edif_edge_hi[jcdif - 1, nhbins - 2]
        )

        ihlow = 0
        ihhigh = 0

        for ihbin in range(1, nobins_out + 1):

            hlow = ebin_edge_out[ihbin - 1]
            hhigh = ebin_edge_out[ihbin]
            hwide = hhigh - hlow

            if ihlow == 0:

                ihlfind = 1
                if True:  # (hlow > edif_cent[0]) and (hlow <= edif_cent[-1]):

                    while ihlfind < nhbins:
                        if (hlow > edif_cent[ihlfind - 1]) and (
                            hlow <= edif_cent[ihlfind]
                        ):
                            ihlow = ihlfind
                            break

                        ihlfind += 1

            if hlow <= edif_cent[0]:

                #           print('sec 1')

                if hhigh > edif_cent[0]:

                    # locate ihhigh

                    ihfind = 1
                    while ihfind < nhbins:
                        if (hhigh > edif_cent[ihfind - 1]) and (
                            hhigh < edif_cent[ihfind]
                        ):
                            ihhigh = ihfind
                            break

                        ihfind += 1

                    nhpoints = (ihhigh) + 2
                    hchunk = hwide / float(nhpoints - 1)

                    for icbin in range(1, nhpoints + 1):
                        euse = hlow + hchunk * float(icbin - 1)

                        if euse <= edif_cent[0]:
                            #                    print('sec 1 a')

                            row_entry = diff_matrix_vec[0] * euse / edif_cent[0]

                        else:

                            #                   print('sec 1 b')

                            icdif = 1

                            while icdif < (ihhigh + 1):

                                if (euse > edif_cent[icdif - 1]) and (
                                    euse <= edif_cent[icdif]
                                ):
                                    # interpolation
                                    row_entry = diff_matrix_vec[icdif - 1] + (
                                        diff_matrix_vec[icdif]
                                        - diff_matrix_vec[icdif - 1]
                                    ) * (euse - edif_cent[icdif - 1]) / (
                                        edif_cent[icdif] - edif_cent[icdif - 1]
                                    )
                                    break
                                icdif += (
                                    1  ############# BB ADDED #########################
                                )
                                # sum up horizontal
                        row_tot[ihbin - 1] += row_entry
                        row_entry = 0.0

                    ##

                    row_tot[ihbin - 1] *= hwide / float(nhpoints)

                    # hwide = horizontal bin width
                    # nhpoints = # of samples used
                    # convert from counts/(unit energy) to counts/bin
                else:

                    row_tot[ihbin - 1] = (
                        diff_matrix_vec[0]
                        * ((hlow + hhigh) / 2.0)
                        / edif_cent[0]
                        * hwide
                    )

            #       if row_tot[ihbin] > 0: print(row_tot[ihbin], ihbin)

            if ihlow >= nhbins:
                #        print('sec 2')

                if hlow > edif_edgeh[nhbins]:
                    row_tot[ihbin - 1] = -1.0
                    ihover = ihbin

                else:

                    if hhigh <= edif_edgeh[nhbins]:

                        row_tot[ihbin - 1] = (
                            diff_matrix_vec[nhbins - 1]
                            * (edif_edgeh[nhbins] - (hlow + hhigh) / 2.0)
                            / (edif_edgeh[nhbins] - edif_cent[nhbins - 1])
                            * hwide
                        )

                    else:

                        row_tot[ihbin - 1] = (
                            ((edif_edgeh[nhbins] - hlow) ** 2)
                            * diff_matrix_vec[nhbins - 1]
                            / (2.0 * (edif_edgeh[nhbins] - edif_cent[nhbins - 1]))
                        )

            #     if row_tot[ihbin] > 0: print(row_tot[ihbin], ihbin)

            elif ihlow >= 1:  # could be zero??

                if hhigh > edif_edgeh[nhbins]:

                    hwide = (
                        edif_edgeh[nhbins] - hlow
                    )  # total width adjusted for active response range
                    nhpoints = (nhbins) - (ihlow) + 2

                    hchunk = hwide / float(nhpoints - 1)

                    for icbin in range(1, nhpoints + 1):

                        euse = hlow + hchunk * float(icbin - 1)

                        icdif = ihlow

                        while icdif < nhbins:  # again check the index

                            if (euse > edif_cent[icdif - 1]) and (
                                euse <= edif_cent[icdif]
                            ):

                                mu = (euse - edif_cent[icdif - 1]) / (
                                    edif_cent[icdif] - edif_cent[icdif - 1]
                                )
                                mu2 = (1 - np.cos(mu * np.pi)) / 2.0
                                row_entry = (
                                    diff_matrix_vec[icdif - 1]
                                    + (
                                        diff_matrix_vec[icdif]
                                        - diff_matrix_vec[icdif - 1]
                                    )
                                    * mu2
                                )

                                break
                            icdif += 1

                        row_tot[ihbin - 1] += row_entry
                        row_entry = 0.0
                        # del row_entry

                    row_tot[ihbin - 1] *= hwide / float(nhpoints)
                    ihlow = nhbins

                else:

                    for ihfind in range(ihlow, nhbins):
                        # SEARCH
                        if (hhigh > edif_cent[ihfind - 1]) and (
                            hhigh <= edif_cent[ihfind]
                        ):
                            ihhigh = ihfind
                    if hhigh > edif_cent[nhbins - 1]:
                        ihhigh = nhbins

                    nhpoints = ihhigh - ihlow + 2

                    if nhpoints < 9:
                        nhpoints = 9

                    hchunk = hwide / float(nhpoints - 1)

                    for icbin in range(1, nhpoints + 1):

                        euse = hlow + hchunk * float(icbin - 1)

                        icdif = ihlow

                        while icdif < ihhigh + 1:  # check

                            # print(icdif, nhbins -2, ihhigh)
                            if icdif <= nhbins - 1:

                                #   print(edif_cent[icdif], euse, edif_cent[icdif+1] )
                                if (euse > edif_cent[icdif - 1]) and (
                                    euse <= edif_cent[icdif]
                                ):

                                    row_entry = diff_matrix_vec[icdif - 1] + (
                                        diff_matrix_vec[icdif]
                                        - diff_matrix_vec[icdif - 1]
                                    ) * (euse - edif_cent[icdif - 1]) / (
                                        edif_cent[icdif] - edif_cent[icdif - 1]
                                    )
                                    break

                            else:

                                row_entry = (
                                    diff_matrix_vec[icdif - 1]
                                    * (hhigh - edif_cent[nhbins - 1])
                                    / (edif_cent[nhbins - 1] - edif_cent[nhbins - 2])
                                )
                                flag = False

                            icdif += 1

                        row_tot[ihbin - 1] += row_entry
                        row_entry = 0.0
                        # del row_entry

                    # print(row_tot[ihbin])
                    row_tot[ihbin - 1] *= hwide / float(
                        nhpoints
                    )  ################# BB ADDED [ihbin] #############
                    ihlow = ihhigh

            if row_tot[ihbin - 1] == -1:

                break

            if ihbin == nobins_out:

                ihover = nobins_out + 1  # + 1

        for ivhsum in range(1, ihover):
            binned_matrix[jcdif - 1, ivhsum - 1] += row_tot[ivhsum - 1]
            row_tot[ivhsum - 1] = 0.0

        row_tot[ihover - 1] = 0.0

    return binned_matrix
