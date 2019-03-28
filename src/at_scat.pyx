#!python
# cython: boundscheck=False
# cython: cdivision=True 
# cython: wraparound=False



cimport numpy as np
import numpy as np
import ftran



DTYPE = np.float32
ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE2_t


cpdef np.ndarray[DTYPE_t, ndim=2] get_at_scat(np.ndarray[DTYPE_t, ndim=1] gx, np.ndarray[DTYPE_t, ndim=1] gy, np.ndarray[DTYPE_t, ndim=1] gz, int il_low, int il_high, float l_frac,  int nobins_out, np.ndarray[DTYPE_t, ndim=1] out_edge, database):



    cdef int itr, i1,i2,i3, ist, N
    cdef float dirx, diry, dirz, az, el
    cdef float dist1, dist2, dist3
    
    cdef float sf, plat,plon
    cdef np.ndarray[DTYPE_t,ndim=1] P, dist_array, 
    cdef np.ndarray[DTYPE2_t,ndim=1] theta_u
    cdef np.ndarray[DTYPE_t,ndim=2] direct_diff_matrix,  tmp_out
    cdef np.ndarray[DTYPE2_t,ndim=2]  phi_u
    cdef np.ndarray[DTYPE_t,ndim=3] out_matrix = np.zeros((len(database.theta_cent)*len(database.double_phi_cent)*2,database.ienerg,nobins_out),dtype=np.float32)


    
    cdef int num_theta = len(database.theta_cent)
    cdef int num_phi   = len(database.double_phi_cent)
    cdef int i,j,k, num_loops
    
    num_loops = num_theta*num_phi*2

    theta_u = database.theta_cent
    phi_u   = database.double_phi_cent


    out_matrix = np.zeros((num_loops,database.ienerg,nobins_out),dtype=np.float32)

    itr = 0
    for i in xrange(num_theta):
        for j in xrange(num_phi):
            for k in xrange(2):



                    dirx, diry, dirz, az, el = ftran.geo_to_space(theta_u[i],phi_u[j,k],gx,gy,gz)

                    sf   = np.arctan(1.)/45.
                    plat = el * sf
                    plon = az * sf
                    P = np.array([np.cos(plat) * np.cos(plon),
                                  np.cos(plat) * np.sin(plon),
                                  np.sin(plat)],np.float32)

                    # Find a new interpolated matrix
                    ist = 0

                    N=len(database.LEND)
                    b1,b2,b3,i1,i2,i3=ftran.trfind(ist,
                                                   P,
                                                   N,
                                                   database.X,
                                                   database.Y,
                                                   database.Z,
                                                   database.LIST,
                                                   database.LPTR,
                                                   database.LEND)
                    i_array= [i1,i2,i3]

                    dist1 = ftran.calc_sphere_dist(az,90.-el,
                                                   database.Azimuth[i1-1],
                                                   database.Zenith[i1-1])
                    dist2 = ftran.calc_sphere_dist(az,90.-el,
                                                   database.Azimuth[i2-1],
                                                   database.Zenith[i2-1])
                    dist3 = ftran.calc_sphere_dist(az,90.-el,
                                                   database.Azimuth[i3-1],
                                                   database.Zenith[i3-1])

                    dist_array = np.array([dist1,dist2,dist3],dtype=np.float32)
                    i1 = i_array[np.argmin(dist_array)]


                    tmpdrm = database.get_rsp(database.millizen[i1-1],
                                              database.milliaz[i1-1])


                    # intergrate the new drm
                    direct_diff_matrix = ftran.echan_integrator(tmpdrm,
                                                                database.epx_lo,
                                                                database.epx_hi,
                                                                database.ichan,
                                                                out_edge)




                    # Now let FORTRAN add the at scat to the direct
                    tmp_out = ftran.sum_at_scat(direct_diff_matrix,
                                                database.at_scat_data[:,:,il_low,i,j],
                                                database.at_scat_data[:,:,il_high,i,j],
                                                l_frac)

                    out_matrix[i] = tmp_out

    return out_matrix.sum(axis=0)
