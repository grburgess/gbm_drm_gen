import numpy as np
import numba as nb


@nb.njit(fastmath=True)
def det(x1,y1,z1,x2,y2,z2,x0,y0,z0):

    return x0*(y1*z2-y2*z1) - y0*(x1*z2-x2*z1) + z0*(x1*y2-x2*y1)


def trfind(nst, p, n, x, y, z, llist, lptr, lend):


    not_found = True

    if (n0 < 1.) or (n < n0):

        # I think I need to draw a random number here
        n0 = np.random.randint(0, n-1)

        # check that this is the right thing to do 
        
        
        pass


    
    while not_found:

        lp = lend[n0]
        n1 = llist[lp]
        lp = lptr[lp]
        nf = llist[lp]
        n1 = nf

        # Find a pair of adjacent neighbors N1,N2 of N0 that define
        # a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.

        if nl > 0 :


            # while loop
            flag = True
            while flag:
            
                if det(x[n0],y[n0],z[n0],x[n1],y[n1],z[n1],xp,yp,zp ) < 0. :
                    lp = lptr[lp]
                    n1 = llist[lp]
                    if n1 == nl:

                        # break here
                        pass
                    
        else:

            pass


spec = [
    ('n', nb.int32),
    
    ('x', nb.float32[:]),
    ('y', nb.float32[:]),
    ('z', nb.float32[:]),

    ('llist', nb.int32[:]),
    ('lptr', nb.int32[:]),
    ('lend', nb.int32[:]),
    
    ('xp', nb.float32),
    ('yp', nb.float32),
    ('zp', nb.float32),
    ('nst', nb.int32),

    ('i1', nb.int32),
    ('i2', nb.int32),
    ('i3', nb.int32),

    ('b1', nb.float32),
    ('b2', nb.float32),
    ('b3', nb.float32),
    
    ]


@nb.jitclass(spec)
class Tri(object):

    def __init__(x, y, z, llist, lptr, lend):

        self.n = len(lend)
        self.x = x
        self.y = y
        self.z = z
        self.llist = llist
        self.lptr = lptr
        self.lend = lend

        self.xp = 0.
        self.yp = 0.
        self.zp = 0.
        self.nst = 0


        self.i1 = 0
        self.i2 = 0
        self.i3 = 0


        self.n0 = 0
        self.lp = 0.
        self.nl = 0
        self.nf = 0
        self.n1 = 0
        

        
    def triangulate(self, p, nst):

        self.xp = p[0]
        self.yp = p[1]
        self.zp = p[2]
        self.nst = nst


        not_found = True

        self.n0 = nst

        if (self.n0 < 1) or (self.n < self.n0):

            # need to check the bounds here
            self.n0 = np.random.randint(1, self.n)

        while not_found:

            self.lp = self.lend[self.n0]
            self.n1 = self.llist[self.lp]
            self.lp = lptr[self.lp]
            self.nf = llist[self.lp]
            self.n1 = self.nf


            if self.nl > 0:
                pass


    def control_6(self):

        self.n2 = self.nf

    def control_7(self):

        self.n3 = self.n0
        self.n1s = self.n1
        self.n2s = self.n2 

    def control_8(self):

        self.b3 = self.det(n1, n2)

        if self.b3 < 0. :

            self.lp = lsptijfoidjf


            if self.llist(self.lp) < 0:

                pass # control 9

            self.lp = self.lptr(self.lp)
            self.n4 = np.abs(self.llist(self.lp))

            if self.det(self.n0, self.n4) < 0. :

                
                



    def lstptr(self, lpl, nb):

        lp = self.lptr[lpl]
        nd = self.llist[lp]

        while (nd != nb) and (lp ):

            lp = self.lptr[lp]

        return lp
            
                
        
    def det(self, i, j):

        return self.xp * (y[i]*z[j] - y[j]*z[i]) - self.yp * (x[i]*z[j] - x[j]*z[i] ) + self.zp * (x[i]*y[j] - x[j]*y[i])
