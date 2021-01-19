import numba as nb
import numpy as np


@nb.njit(fastmath=False, parallel=False)
def closest(point_assume_f, grid_points_list_f):

    N = grid_points_list_f.shape[0]

    i1 = 0

    maxi = 11.0

    for i in range(N):

        distance = np.arccos(
            grid_points_list_f[i, 0] * point_assume_f[0]
            + grid_points_list_f[i, 1] * point_assume_f[1]
            + grid_points_list_f[i, 2] * point_assume_f[2]
        )
        if distance < maxi:
            i1 = i
            maxi = distance
    return i1

@nb.njit(fastmath=False, parallel=False)
def trfind(point_assume_f, grid_points_list_f):

    N = grid_points_list_f.shape[0]

    i1 = 0
    i2 = 0
    i3 = 0

    mini = 9.0
    midi = 10.0
    maxi = 11.0

    for i in range(N):

        distance = np.arccos(
            grid_points_list_f[i, 0] * point_assume_f[0]
            + grid_points_list_f[i, 1] * point_assume_f[1]
            + grid_points_list_f[i, 2] * point_assume_f[2]
        )
        if distance < maxi:
            if distance > mini:

                if distance < midi:
                    maxi = midi
                    i3 = i2
                    midi = distance
                    i2 = i

                else:
                    maxi = distance
                    i3 = i

            else:

                maxi = midi
                midi = mini
                mini = distance
                i3 = i2
                i2 = i1
                i1 = i

    weights = calc_weights_numba(
        grid_points_list_f[i1],
        grid_points_list_f[i2],
        grid_points_list_f[i3],
        point_assume_f,
    )
    weights.sort()
    return weights[2], weights[1], weights[0], i1, i2, i3


@nb.njit(fastmath=True)
def calc_weights_numba(p1, p2, p3, p_find):
    """
    ###################### Weights from https://codeplea.com/triangular-interpolation ############
    p1_lat = np.arcsin(p1[2])
    p1_lon = np.arctan2(p1[1],p1[0])

    p2_lat = np.arcsin(p2[2])
    p2_lon = np.arctan2(p2[1],p2[0])

    p3_lat = np.arcsin(p3[2])
    p3_lon = np.arctan2(p3[1],p3[0])

    pf_lat = np.arcsin(p_find[2])
    pf_lon = np.arctan2(p_find[1],p_find[0])

    W1 = ((p2_lat-p3_lat)*(pf_lon-p3_lon)+(p3_lon-p2_lon)*(pf_lat-p3_lat))/((p2_lat-p3_lat)*(p1_lon-p3_lon)+(p3_lon-p2_lon)*(p1_lat-p3_lat))
    W2 = ((p3_lat-p1_lat)*(pf_lon-p3_lon)+(p1_lon-p3_lon)*(pf_lat-p3_lat))/((p2_lat-p3_lat)*(p1_lon-p3_lon)+(p3_lon-p2_lon)*(p1_lat-p3_lat))
    W3 = 1-W1-W2
    """

    ###################### Weights from ftran code. But NOT set to 0 when they are negative ############

    w = np.zeros(3)

    w[0] = (
        p_find[0] * (p1[1] * p2[2] - p2[1] * p1[2])
        - p_find[1] * (p1[0] * p2[2] - p2[0] * p1[2])
        + p_find[2] * (p1[0] * p2[1] - p2[0] * p1[1])
    )
    w[1] = (
        p_find[0] * (p2[1] * p3[2] - p3[1] * p2[2])
        - p_find[1] * (p2[0] * p3[2] - p3[0] * p2[2])
        + p_find[2] * (p2[0] * p3[1] - p3[0] * p2[1])
    )
    w[2] = (
        p_find[0] * (p3[1] * p1[2] - p1[1] * p3[2])
        - p_find[1] * (p3[0] * p1[2] - p1[0] * p3[2])
        + p_find[2] * (p3[0] * p1[1] - p1[0] * p3[1])
    )
    return np.abs(w)


@nb.njit(fastmath=True)
def geocoords(theta_geo, phi_geo, theta_source, phi_source):

    gx = np.empty(3)
    gy = np.empty(3)
    gz = np.empty(3)
    sl = np.empty(3)

    gz[0] = np.sin(theta_geo) * np.cos(phi_geo)
    gz[1] = np.sin(theta_geo) * np.sin(phi_geo)
    gz[2] = np.cos(theta_geo)

    gzr = np.sqrt(gz[0] * gz[0] + gz[1] * gz[1] + gz[2] * gz[2])
    gz[0] = gz[0] / gzr
    gz[1] = gz[1] / gzr
    gz[2] = gz[2] / gzr

    sl[0] = np.sin(theta_source) * np.cos(phi_source)
    sl[1] = np.sin(theta_source) * np.sin(phi_source)
    sl[2] = np.cos(theta_source)

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

    dtr = np.pi / 180.0

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
def calc_sphere_dist(ra1, dec1, ra2, dec2, dtr):

    y = np.sqrt(
        (np.cos(dec2 * dtr) * np.sin((ra1 - ra2) * dtr)) ** 2
        + (
            np.cos(dec1 * dtr) * np.sin(dec2 * dtr)
            - np.sin(dec1 * dtr) * np.cos(dec2 * dtr) *
            np.cos((ra1 - ra2) * dtr)
        )
        ** 2
    )
    x = np.sin(dec1 * dtr) * np.sin(dec2 * dtr) + np.cos(dec1 * dtr) * np.cos(
        dec2 * dtr
    ) * np.cos((ra1 - ra2) * dtr)
    return np.arctan2(y, x) / dtr


@nb.njit(fastmath=True)
def highres_ephoton_interpolator(
    ebin_edge_in, ein, matrix, edif_edge_lo, edif_edge_hi, nhbins
):

    nobins_in = len(ebin_edge_in)

    nvbins = len(ein) - 1

    new_epx_lo = np.zeros((nobins_in, 64))
    new_epx_hi = np.zeros((nobins_in, 64))

    diff_matrix = np.zeros((nobins_in, 64))

    ivfind = 1

    for i in range(nobins_in):

        for j in range(ivfind, 70):

            if (ebin_edge_in[i] >= ein[j - 1]) and ebin_edge_in[i] < ein[j]:

                ivfind = j

                mu = (np.log(ebin_edge_in[i]) - np.log(ein[ivfind - 1])) / (
                    np.log(ein[ivfind]) - np.log(ein[ivfind - 1])
                )
                if (mu < 0.0) and (mu > -1e-5):
                    mu = 0.0
                elif (mu > 1.0) and (mu < 1.00001):
                    mu = 1.0

                for k in range(nhbins):
                    # print(ivfind, k)

                    new_epx_lo[i, k] = (
                        edif_edge_lo[ivfind - 1, k] /
                        ein[ivfind - 1] * (1 - mu)
                        + edif_edge_lo[ivfind, k] / ein[ivfind] * mu
                    ) * ebin_edge_in[i]

                    new_epx_hi[i, k] = (
                        edif_edge_hi[ivfind - 1, k] /
                        ein[ivfind - 1] * (1 - mu)
                        + edif_edge_hi[ivfind, k] / ein[ivfind] * mu
                    ) * ebin_edge_in[i]

                    diff_matrix[i, k] = (
                        matrix[ivfind - 1, k] *
                        (1 - mu) + matrix[ivfind, k] * mu
                    )

    return new_epx_lo, new_epx_hi, diff_matrix


@nb.njit(fastmath=True)
def atscat_highres_ephoton_interpolator(ebin_edge_in, ein, matrix):

    nobins_in = len(ebin_edge_in) - 1
    nvbins = len(ein) - 1
    nobins_out = matrix.shape[1]

    new_matrix = np.zeros((nobins_in, nobins_out))

    ivfind = 0

    # max_i = np.searchsorted(ebin_edge_in, ein[-1])
    # min_i = np.searchsorted(ebin_edge_in, ein[0])

    for i in range(nobins_in):

        for j in range(ivfind, nvbins):
            if (ebin_edge_in[i] >= ein[j]) and (ebin_edge_in[i] < ein[j + 1]):

                ivfind = j

                mu = (np.log(ebin_edge_in[i]) - np.log(ein[ivfind])) / (
                    np.log(ein[ivfind + 1]) - np.log(ein[ivfind])
                )
                if (mu < 0.0) and (mu > -1e-5):
                    mu = 0.0
                elif (mu > 1.0) and (mu < 1.00001):
                    mu = 1.0

                for k in range(nobins_out):
                    new_matrix[i, k] = (
                        matrix[ivfind, k] * (1 - mu) +
                        matrix[ivfind + 1, k] * mu
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
    diff_matrix_vec = np.zeros(nhbins)

    edif_edgeh = np.zeros(nhbins + 1)
    edif_cent = np.zeros(nhbins)
    # first is a loop over the photon energies
    # row_entry = 0.
    # ihover =0
    for jcdif in range(1, nobins_in + 1):

        total = 0
        for ivh in range(1, nhbins + 1):

            edif_cent[ivh - 1] = (
                edif_edge_lo[jcdif - 1, ivh - 1] +
                edif_edge_hi[jcdif - 1, ivh - 1]
            ) / 2.0

            total += edif_cent[ivh - 1]
            if edif_cent[ivh - 1] > 0:

                diff_matrix_vec[ivh - 1] = diff_matrix[jcdif - 1, ivh - 1] / (
                    edif_edge_hi[jcdif - 1, ivh - 1] -
                    edif_edge_lo[jcdif - 1, ivh - 1]
                )
                edif_edgeh[ivh - 1] = edif_edge_hi[jcdif - 1, ivh - 1]

                edif_edgeh[nhbins] = edif_edge_hi[jcdif - 1, nhbins - 1] + (
                    edif_edge_hi[
                        jcdif - 1,
                        nhbins - 1,
                    ]
                    - edif_edge_hi[jcdif - 1, nhbins - 2]
                )

        ihlow = 0
        ihhigh = 0
        if total == 0:
            continue

        for ihbin in range(1, nobins_out + 1):

            hlow = ebin_edge_out[ihbin - 1]
            hhigh = ebin_edge_out[ihbin]
            hwide = hhigh - hlow

            if ihlow == 0:

                ihlow = np.searchsorted(edif_cent, hlow)

            if hlow <= edif_cent[0]:

                #           print('sec 1')

                if hhigh > edif_cent[0]:

                    ihhigh = np.searchsorted(edif_cent, hhigh)

                    nhpoints = (ihhigh) + 2
                    hchunk = hwide / float(nhpoints - 1)

                    for icbin in range(1, nhpoints + 1):
                        euse = hlow + hchunk * float(icbin - 1)

                        if euse <= edif_cent[0]:
                            #                    print('sec 1 a')

                            row_entry = diff_matrix_vec[0] * \
                                euse / edif_cent[0]

                        else:

                            #                   print('sec 1 b')

                            icdif = 1

                            icdif = np.searchsorted(edif_cent, euse)
                            if icdif < (ihhigh + 1):
                                row_entry = diff_matrix_vec[icdif - 1] + (
                                    diff_matrix_vec[icdif] -
                                    diff_matrix_vec[icdif - 1]
                                ) * (euse - edif_cent[icdif - 1]) / (
                                    edif_cent[icdif] - edif_cent[icdif - 1]
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

            elif ihlow >= 1:  # could be zero??

                if hhigh > edif_edgeh[nhbins]:

                    hwide = (
                        edif_edgeh[nhbins] - hlow
                    )  # total width adjusted for active response range
                    nhpoints = (nhbins) - (ihlow) + 2

                    hchunk = hwide / float(nhpoints - 1)

                    for icbin in range(1, nhpoints + 1):

                        euse = hlow + hchunk * float(icbin - 1)

                        icdif = np.searchsorted(
                            edif_cent[ihlow:], euse) + ihlow

                        if icdif < nhbins:  # again check the index

                            mu = (euse - edif_cent[icdif - 1]) / (
                                edif_cent[icdif] - edif_cent[icdif - 1]
                            )
                            mu2 = (1 - np.cos(mu * np.pi)) / 2.0
                            row_entry = (
                                diff_matrix_vec[icdif - 1]
                                + (diff_matrix_vec[icdif] -
                                   diff_matrix_vec[icdif - 1])
                                * mu2
                            )

                        row_tot[ihbin - 1] += row_entry
                        row_entry = 0.0
                        # del row_entry

                    row_tot[ihbin - 1] *= hwide / float(nhpoints)
                    ihlow = nhbins

                else:

                    if hhigh > edif_cent[nhbins - 1]:
                        ihhigh = nhbins

                    else:
                        ihhigh = np.searchsorted(
                            edif_cent[ihlow:], hhigh) + ihlow

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

                            icdif += 1

                        row_tot[ihbin - 1] += row_entry
                        row_entry = 0.0
                        # del row_entry

                    # print(row_tot[ihbin])
                    row_tot[ihbin - 1] *= hwide / float(
                        nhpoints
                    )  # BB ADDED [ihbin] #############
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
