import numpy as np
import numba as nb



spec = [
    ("n", nb.int32),
    ("x", nb.float32[:]),
    ("y", nb.float32[:]),
    ("z", nb.float32[:]),
    ("list", nb.int32[:]),
    ("lptr", nb.int32[:]),
    ("lend", nb.int32[:]),
    ("xp", nb.float32),
    ("yp", nb.float32),
    ("zp", nb.float32),
    ("nst", nb.int32),
    ("i1", nb.int32),
    ("i2", nb.int32),
    ("i3", nb.int32),
    ("b1", nb.float32),
    ("b2", nb.float32),
    ("b3", nb.float32),
    ("n1", nb.int32),
    ("n2", nb.int32),
    ("n3", nb.int32),
    ("nl", nb.int32),
    ("nf", nb.int32),
    ("n0", nb.int32),
    ("n1s", nb.int32),
    ("n2s", nb.int32),
    ("not_found", nb.bool),
    ("eps", nb.float32),
    ("tol", nb.float32),
]


@nb.jitclass(spec)
class Tri(object):
    def __init__(self, x, y, z, llist, lptr, lend):

        self.n = len(lend)
        self.x = x
        self.y = y
        self.z = z
        self.list = llist
        self.lptr = lptr
        self.lend = lend

        self.xp = 0.0
        self.yp = 0.0
        self.zp = 0.0
        self.nst = 0

        self.i1 = 0
        self.i2 = 0
        self.i3 = 0

        self.n0 = 0
        self.lp = 0
        self.nl = 0
        self.nf = 0
        self.n1 = 0
        self.n2 = 0
        self.n3 = 0

        self.eps = np.finfo(np.float32).eps
        self.tol = 100.0 * self.eps

    def triangulate(self, p, nst):

        self.xp = p[0]
        self.yp = p[1]
        self.zp = p[2]
        self.nst = nst

        self.not_found = True

        self.n0 = nst

        if (self.n0 < 1) or (self.n < self.n0):

            # need to check the bounds here
            self.n0 = np.random.randint(1, self.n)

        self.control_2()

    def control_2(self):

        while self.not_found:

            self.lp = self.lend[self.n0]
            self.n1 = self.list[self.lp]
            self.lp = self.lptr[self.lp]
            self.nf = self.list[self.lp]
            self.n1 = self.nf

            if self.nl > 0:

                self.control_3()

            else:

                self.nl = -self.nl

                if self.det(self.n0, self.nf) < 0.0:

                    self.n1 = self.n0
                    self.n2 = self.nf

                    self.control_9()

                if self.det(self.nl, self.n0) < 0.0:

                    self.n1 = self.nl
                    self.n2 = self.n0

                    self.control_9()

            self.control_4()

        return

    def control_3(self):

        if self.det(self.n0, self.n1) < 0.0:

            self.lp = self.lptr[self.lp]
            self.n1 = self.list[self.lp]

            if self.n1 == self.nl:

                self.control_6()

            self.control_3()

    def control_4(self):

        self.lp = self.lptr[self.lp]

        self.n2 = np.abs(self.list[self.lp])

        if self.det(self.n0, self.n2) < 0.0:

            self.control_7()

        self.n1 = self.n2

        if self.n1 != self.nl:

            self.control_4()

        if self.det(self.n0, self.n4) < 0.0:

            self.control_6()

        if np.abs(
            self.x[self.n0] * self.xp
            + self.y[self.n0] * self.yp
            + self.z[self.n0 * self.zp]
        ) < 1 - (4.0 * self.eps):

            # need some kind of STORE thing here which shoves stuff in a common block

            flag = self.det(self.n1, self.n0) < 0.0

            while flag:

                self.lp = self.lptr[self.lp]
                self.n1 = np.abs(self.list[self.lp])

                if self.n1 == self.nl:

                    self.i1 = 0
                    self.i2 = 0
                    self.i3 = 0

                    # return
                    self.not_found = False

                    self.control_2()

                flag = self.det(self.n1, self.n0) < 0.0

        self.n0 = self.n1
        self.control_2()

    def control_6(self):

        self.n2 = self.nf

    def control_7(self):

        self.n3 = self.n0
        self.n1s = self.n1
        self.n2s = self.n2

    def control_8(self):

        self.b3 = self.det(self.n1, self.n2)

        if self.b3 < 0.0:

            self.lp = self.lstptr(self.lend[self.n2], self.n1)

            if self.list[self.lp] < 0:

                pass  # control 9

            self.lp = self.lptr[self.lp]
            self.n4 = np.abs(self.list[self.lp])

            if self.det(self.n0, self.n4) < 0.0:

                self.n3 = self.n2
                self.n2 = self.n4
                self.n1s = self.n1

                if (self.n2 != self.n2s) and (self.n2 != self.n0):

                    self.control_8()

            else:

                self.n3 = self.n1
                self.n1 = self.n4
                self.n2s = self.n2

                if (self.n1 != self.n1s) and (self.n1 != self.n0):

                    self.control_8()

            self.n0 = np.random.randint(1, self.n)

            self.control_2()

        if self.b3 >= self.eps:

            self.b1 = self.det(self.n2, self.n3)
            self.b3 = self.det(self.n3, self.n1)

            if (self.b1 < -self.tot) or (self.b2 < -self.tol):

                self.n0 = np.random.randint(1, self.n)

                self.control_2()

        self.i1 = self.n1
        self.i2 = self.n2
        self.i3 = self.n3

        self.b1 = np.max(self.b1, 0.0)
        self.b2 = np.max(self.b2, 0.0)

        self.not_found = False
        self.control_2()

    def control_9(self):

        self.n1s = self.n1
        self.n2s = self.n2
        self.nl = 0

        self.control_10()

    def control_10(self):

        self.lp = self.lend[self.n2]
        self.lp = self.lptr[self.lp]

        self.next = self.list[self.lp]

        if self.det(self.n2, self.next):

            s12 = (
                self.x[self.n1] * self.x[self.n2]
                + self.y[self.n1] * self.y[self.n2]
                + self.z[self.n1] * self.z[self.n2]
            )

            self.q = [
                self.x[self.n1] - s12 * self.x[self.n2],
                self.y[self.n1] - s12 * self.y[self.n2],
                self.z[self.n1] - s12 * self.z[self.n2],
            ]

            if self.xp * self.q[0] + self.yp * self.q[1] + self.zp ** self.q[2] >= 0.0:

                self.control_11()

            if self.xp * self.q[0] + self.yp * self.q[1] + self.zp ** self.q[2] >= 0.0:

                self.control_11()

            self.nl = self.n2

        self.n1 = self.n2
        self.n2 = self.next

        if self.n2 != self.n1s:
            self.control_10()

        self.i1 = self.n1s
        self.i2 = self.n1s
        self.i3 = 0

        #### RETURN
        self.not_found = False
        self.control_2()

    def control_11(self):

        self.nf = self.n2

        if self.nl == 0:

            self.control_12()

    def control_12(self):

        self.lp = self.lend[self.n1]
        self.next = -self.list[self.lp]

        if 0.0 < self.det(self.next, self.n1):

            s12 = (
                self.x[self.n1] * self.x[self.n2]
                + self.y[self.n1] * self.y[self.n2]
                + self.z[self.n1] * self.z[self.n2]
            )

            self.q = [
                self.x[self.n1] - s12 * self.x[self.n2],
                self.y[self.n1] - s12 * self.y[self.n2],
                self.z[self.n1] - s12 * self.z[self.n2],
            ]

            if self.xp * self.q[0] + self.yp * self.q[1] + self.zp ** self.q[2] >= 0.0:

                self.control_13()

            if self.xp * self.q[0] + self.yp * self.q[1] + self.zp ** self.q[2] >= 0.0:

                self.control_13()

            self.nf = self.n1

        self.n2 = self.n1
        self.n1 = self.next

        if self.n1 != self.n1s:

            self.control_12()

        self.i1 = self.n1
        self.i2 = self.n1
        self.i3 = 0

        ###### return
        self.not_found = False
        self.control_2()

    def control_13(self):

        self.nl = self.n1

        self.i1 = self.nf
        self.i2 = self.nl
        self.i3 = 0

        # return
        self.not_found = False
        self.control_2()

    #       INTEGER LP, ND
    # C
    #       LP = LPTR(LPL)
    #     1 ND = LIST(LP)
    #         IF (ND .EQ. NB) GO TO 2
    #         LP = LPTR(LP)
    #         IF (LP .NE. LPL) GO TO 1
    # C
    #     2 LSTPTR = LP
    #       RETURN
    #       END
    #       INTEGER FUNCTION NBCNT (LPL,LPTR)
    #       INTEGER LPL, LPTR(*)
    # C

    def lstptr(self, lpl, nb):

        lp = self.lptr[lpl]
        nd = self.list[lp]

        while (nd != nb) and (lp != lpl):

            lp = self.lptr[lp]

        return lp

    def det(self, i, j):

        return (
            self.xp * (self.y[i] * self.z[j] - self.y[j] * self.z[i])
            - self.yp * (self.x[i] * self.z[j] - self.x[j] * self.z[i])
            + self.zp * (self.x[i] * self.y[j] - self.x[j] * self.y[i])
        )
