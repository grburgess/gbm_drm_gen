subroutine trfind ( nst, p, n, x, y, z, list, lptr, lend, b1, b2, b3, i1, &
  i2, i3 )

!*****************************************************************************80
!
!! TRFIND locates a point relative to a triangulation.
!
!  Discussion:
!
!    This subroutine locates a point P relative to a triangulation
!    created by TRMESH.  If P is contained in
!    a triangle, the three vertex indexes and barycentric 
!    coordinates are returned.  Otherwise, the indexes of the
!    visible boundary nodes are returned.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 772: STRIPACK,
!    Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 3, September 1997, pages 416-434.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, index of a node at which TRFIND begins
!    its search.  Search time depends on the proximity of this node to P.
!
!    Input, real ( kind = 8 ) P(3), the x, y, and z coordinates (in that order)
!    of the point P to be located.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the coordinates of the
!    triangulation nodes (unit vectors).
!
!    Input, integer ( kind = 4 ) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N), the 
!    data structure defining the triangulation, created by TRMESH.
!
!    Output, real ( kind = 8 ) B1, B2, B3, the unnormalized barycentric
!    coordinates of the central projection of P onto the underlying planar
!    triangle if P is in the convex hull of the nodes.  These parameters 
!    are not altered if I1 = 0.
!
!    Output, integer ( kind = 4 ) I1, I2, I3, the counterclockwise-ordered 
!    vertex indexes of a triangle containing P if P is contained in a triangle.
!    If P is not in the convex hull of the nodes, I1 and I2 are the rightmost 
!    and leftmost (boundary) nodes that are visible from P, and I3 = 0.  (If 
!    all boundary nodes are visible from P, then I1 and I2 coincide.)
!    I1 = I2 = I3 = 0 if P and all of the nodes are coplanar (lie on a 
!    common great circle.
!
!  Local parameters:
!
!    EPS =      Machine precision
!    IX,IY,IZ = Integer seeds for JRAND
!    LP =       LIST pointer
!    N0,N1,N2 = Nodes in counterclockwise order defining a
!               cone (with vertex N0) containing P, or end-
!               points of a boundary edge such that P Right
!               N1->N2
!    N1S,N2S =  Initially-determined values of N1 and N2
!    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
!    NEXT =     Candidate for I1 or I2 when P is exterior
!    NF,NL =    First and last neighbors of N0, or first
!               (rightmost) and last (leftmost) nodes
!               visible from P when P is exterior to the
!               triangulation
!    PTN1 =     Scalar product <P,N1>
!    PTN2 =     Scalar product <P,N2>
!    Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in
!               the boundary traversal when P is exterior
!    S12 =      Scalar product <N1,N2>
!    TOL =      Tolerance (multiple of EPS) defining an upper
!               bound on the magnitude of a negative bary-
!               centric coordinate (B1 or B2) for P in a
!               triangle -- used to avoid an infinite number
!               of restarts with 0 <= B3 < EPS and B1 < 0 or
!               B2 < 0 but small in magnitude
!    XP,YP,ZP = Local variables containing P(1), P(2), and P(3)
!    X0,Y0,Z0 = Dummy arguments for DET
!    X1,Y1,Z1 = Dummy arguments for DET
!    X2,Y2,Z2 = Dummy arguments for DET
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) det
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), save :: ix = 1
  integer ( kind = 4 ), save :: iy = 2
  integer ( kind = 4 ), save :: iz = 3
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(6*(n-2))
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(6*(n-2))
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1s
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2s
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nst
  real ( kind = 8 ) p(3)
  real ( kind = 8 ) ptn1
  real ( kind = 8 ) ptn2
  real ( kind = 8 ) q(3)
  real ( kind = 8 ) s12
  real ( kind = 8 ) store
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) yp
  real ( kind = 8 ) z(n)
  real ( kind = 8 ) z0
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) zp
!
!  Statement function:
!
!  DET(X1,...,Z0) >= 0 if and only if (X0,Y0,Z0) is in the
!  (closed) left hemisphere defined by the plane containing (0,0,0),
!  (X1,Y1,Z1), and (X2,Y2,Z2), where left is defined relative to an 
!  observer at (X1,Y1,Z1) facing (X2,Y2,Z2).
!
  det (x1,y1,z1,x2,y2,z2,x0,y0,z0) = x0*(y1*z2-y2*z1) &
       - y0*(x1*z2-x2*z1) + z0*(x1*y2-x2*y1)
!
!  Initialize variables.
!
  xp = p(1)
  yp = p(2)
  zp = p(3)
  n0 = nst

  if ( n0 < 1 .or. n < n0 ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Compute the relative machine precision EPS and TOL.
!
  eps = epsilon ( eps )
  tol = 100.0D+00 * eps
!
!  Set NF and NL to the first and last neighbors of N0, and initialize N1 = NF.
!
2 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if ( 0 < nl ) then
!
!  N0 is an interior node.  Find N1.
!
3   continue

    if ( det ( x(n0),y(n0),z(n0),x(n1),y(n1),z(n1),xp,yp,zp ) < 0.0D+00 ) then
      lp = lptr(lp)
      n1 = list(lp)
      if ( n1 == nl ) then
        go to 6
      end if
      go to 3
    end if

  else
!
!  N0 is a boundary node.  Test for P exterior.
!
    nl = -nl
!
!  Is P to the right of the boundary edge N0->NF?
!
    if ( det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf), xp,yp,zp) < 0.0D+00 ) then
      n1 = n0
      n2 = nf
      go to 9
    end if
!
!  Is P to the right of the boundary edge NL->N0?
!
    if ( det(x(nl),y(nl),z(nl),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0D+00 ) then
      n1 = nl
      n2 = n0
      go to 9
    end if

  end if
!
!  P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the
!  next neighbor of N0 (following N1).
!
4 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    if ( det(x(n0),y(n0),z(n0),x(n2),y(n2),z(n2),xp,yp,zp) < 0.0D+00 ) then
      go to 7
    end if

    n1 = n2

    if ( n1 /= nl ) then
      go to 4
    end if

  if ( det ( x(n0), y(n0), z(n0), x(nf), y(nf), z(nf), xp, yp, zp ) &
    < 0.0D+00 ) then
    go to 6
  end if
!
!  P is left of or on arcs N0->NB for all neighbors NB
!  of N0.  Test for P = +/-N0.
!
  if ( store ( abs ( x(n0 ) * xp + y(n0) * yp + z(n0) * zp) ) &
    < 1.0D+00 - 4.0D+00 * eps ) then
!
!  All points are collinear iff P Left NB->N0 for all
!  neighbors NB of N0.  Search the neighbors of N0.
!  Note:  N1 = NL and LP points to NL.
!
    do

      if ( det(x(n1),y(n1),z(n1),x(n0),y(n0),z(n0),xp,yp,zp) < 0.0D+00 ) then
        exit
      end if

      lp = lptr(lp)
      n1 = abs ( list(lp) )

      if ( n1 == nl ) then
        i1 = 0
        i2 = 0
        i3 = 0
        return
      end if

    end do

  end if
!
!  P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1 and start over.
!
  n0 = n1
  go to 2
!
!  P is between arcs N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in a wedge defined by geodesics N0-N1 and
!  N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge-hopping loop:
!
8 continue

  b3 = det ( x(n1),y(n1),z(n1),x(n2),y(n2),z(n2),xp,yp,zp )

  if ( b3 < 0.0D+00 ) then
!
!  Set N4 to the first neighbor of N2 following N1 (the
!  node opposite N2->N1) unless N1->N2 is a boundary arc.
!
    lp = lstptr ( lend(n2), n1, list, lptr )

    if ( list(lp) < 0 ) then
      go to 9
    end if

    lp = lptr(lp)
    n4 = abs ( list(lp) )
!
!  Define a new arc N1->N2 which intersects the geodesic N0-P.
!
    if ( det ( x(n0),y(n0),z(n0),x(n4),y(n4),z(n4),xp,yp,zp ) < 0.0D+00 ) then
      n3 = n2
      n2 = n4
      n1s = n1
      if ( n2 /= n2s .and. n2 /= n0 ) then
        go to 8
      end if
    else
      n3 = n1
      n1 = n4
      n2s = n2
      if ( n1 /= n1s .and. n1 /= n0 ) then
        go to 8
      end if
    end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
    n0 = jrand ( n, ix, iy, iz )
    go to 2

  end if
!
!  P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear
!  or P is close to -N0.
!
  if ( b3 >= eps ) then
!
!  B3 /= 0.
!
    b1 = det(x(n2),y(n2),z(n2),x(n3),y(n3),z(n3),xp,yp,zp)
    b2 = det(x(n3),y(n3),z(n3),x(n1),y(n1),z(n1),xp,yp,zp)
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  else
!
!  B3 = 0 and thus P lies on N1->N2. Compute
!  B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
!
    b3 = 0.0D+00
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
    ptn1 = xp * x(n1) + yp * y(n1) + zp * z(n1)
    ptn2 = xp * x(n2) + yp * y(n2) + zp * z(n2)
    b1 = ptn1 - s12 * ptn2
    b2 = ptn2 - s12 * ptn1
!
!  Restart with N0 randomly selected.
!
    if ( b1 < -tol .or. b2 < -tol ) then
      n0 = jrand ( n, ix, iy, iz )
      go to 2
    end if

  end if
!
!  P is in (N1,N2,N3).
!
  i1 = n1
  i2 = n2
  i3 = n3
  b1 = max ( b1, 0.0D+00 )
  b2 = max ( b2, 0.0D+00 )
  return
!
!  P Right N1->N2, where N1->N2 is a boundary edge.
!  Save N1 and N2, and set NL = 0 to indicate that
!  NL has not yet been found.
!
9 continue

  n1s = n1
  n2s = n2
  nl = 0
!
!  Counterclockwise Boundary Traversal:
!
10 continue

  lp = lend(n2)
  lp = lptr(lp)
  next = list(lp)

  if ( det(x(n2),y(n2),z(n2),x(next),y(next),z(next),xp,yp,zp) >= 0.0D+00 ) then
!
!  N2 is the rightmost visible node if P Forward N2->N1
!  or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2.
!
    s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)

    q(1) = x(n1) - s12 * x(n2)
    q(2) = y(n1) - s12 * y(n2)
    q(3) = z(n1) - s12 * z(n2)

    if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0D+00 ) then
      go to 11
    end if

    if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0D+00 ) then
      go to 11
    end if
!
!  N1, N2, NEXT, and P are nearly collinear, and N2 is
!  the leftmost visible node.
!
    nl = n2
  end if
!
!  Bottom of counterclockwise loop:
!
  n1 = n2
  n2 = next

  if ( n2 /= n1s ) then
    go to 10
  end if
!
!  All boundary nodes are visible from P.
!
  i1 = n1s
  i2 = n1s
  i3 = 0
  return
!
!  N2 is the rightmost visible node.
!
11 continue

  nf = n2

  if ( nl == 0 ) then
!
!  Restore initial values of N1 and N2, and begin the search
!  for the leftmost visible node.
!
    n2 = n2s
    n1 = n1s
!
!  Clockwise Boundary Traversal:
!
12  continue

    lp = lend(n1)
    next = -list(lp)

    if ( 0.0D+00 <= &
      det ( x(next), y(next), z(next), x(n1), y(n1), z(n1), xp, yp, zp )  ) then
!
!  N1 is the leftmost visible node if P or NEXT is
!  forward of N1->N2.  Compute Q = N1 X (N2 X N1).
!
      s12 = x(n1) * x(n2) + y(n1) * y(n2) + z(n1) * z(n2)
      q(1) = x(n2) - s12 * x(n1)
      q(2) = y(n2) - s12 * y(n1)
      q(3) = z(n2) - s12 * z(n1)

      if ( xp * q(1) + yp * q(2) + zp * q(3) >= 0.0D+00 ) then
        go to 13
      end if

      if ( x(next) * q(1) + y(next) * q(2) + z(next) * q(3) >= 0.0D+00 ) then
        go to 13
      end if
!
!  P, NEXT, N1, and N2 are nearly collinear and N1 is the rightmost 
!  visible node.
!
      nf = n1
    end if
!
!  Bottom of clockwise loop:
!
    n2 = n1
    n1 = next

    if ( n1 /= n1s ) then
      go to 12
    end if
!
!  All boundary nodes are visible from P.
!
    i1 = n1
    i2 = n1
    i3 = 0
    return
!
!  N1 is the leftmost visible node.
!
13   continue

    nl = n1

  end if
!
!  NF and NL have been found.
!
  i1 = nf
  i2 = nl
  i3 = 0

  return
end
subroutine trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )
