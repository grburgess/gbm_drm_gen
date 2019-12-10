! This code is an adaptation of the code written by Colleen Wilson-Hodge, Geoff Pendalton
! and others using the gbmrsp_release_v2p0 routines. It has been modified for speed and 
! readability. TRFIND is a modified version of R.J. Renka's original source
! ( and redistributed with an ACM liscense). 

SUBROUTINE TRFIND (NST,P,N,X,Y,Z,LIST,LPTR,LEND, B1,B2,B3,I1,I2,I3)
  INTEGER NST, N, LIST(*), LPTR(*), LEND(*)
  REAL    P(3), X(*), Y(*), Z(*)
  INTEGER, INTENT(OUT) :: I1, I2, I3
  REAL, INTENT(OUT) :: B1, B2, B3
  
  INTEGER JRAND, LSTPTR
  INTEGER IX, IY, IZ, LP, N0, N1, N1S, N2, N2S, N3, N4,        NEXT, NF, NL
  REAL    STORE
  REAL    DET, EPS, PTN1, PTN2, Q(3), S12, TOL, XP, YP,        ZP
  REAL    X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2
  
  SAVE    IX, IY, IZ
  DATA    IX/1/, IY/2/, IZ/3/
  
  DET (X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0) = X0*(Y1*Z2-Y2*Z1)     - Y0*(X1*Z2-X2*Z1) + Z0*(X1*Y2-X2*Y1)
  
  XP = P(1)
  YP = P(2)
  ZP = P(3)
  N0 = NST
  IF (N0 .LT. 1  .OR.  N0 .GT. N)  N0 = JRAND(N, IX,IY,IZ )
  
  EPS = 1.E0
1 EPS = EPS/2.E0
  IF (STORE(EPS+1.E0) .GT. 1.E0) GO TO 1
  EPS = 2.E0*EPS
  TOL = 100.E0*EPS
  
2 LP = LEND(N0)
  NL = LIST(LP)
  LP = LPTR(LP)
  NF = LIST(LP)
  N1 = NF
  
  IF (NL .GT. 0) THEN
     
3    IF ( DET(X(N0),Y(N0),Z(N0),X(N1),Y(N1),Z(N1),          XP,YP,ZP) .LT. 0. ) THEN
        LP = LPTR(LP)
        N1 = LIST(LP)
        IF (N1 .EQ. NL) GO TO 6
        GO TO 3
     ENDIF
  ELSE
     
     NL = -NL
     IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF),          XP,YP,ZP) .LT. 0. ) THEN
        
        N1 = N0
        N2 = NF
        GO TO 9
     ENDIF
     IF ( DET(X(NL),Y(NL),Z(NL),X(N0),Y(N0),Z(N0),           XP,YP,ZP) .LT. 0. ) THEN
        
        N1 = NL
        N2 = N0
        GO TO 9
     ENDIF
  ENDIF
  
4 LP = LPTR(LP)
  N2 = ABS(LIST(LP))
  IF ( DET(X(N0),Y(N0),Z(N0),X(N2),Y(N2),Z(N2),        XP,YP,ZP) .LT. 0. ) GO TO 7
  N1 = N2
  IF (N1 .NE. NL) GO TO 4
  IF ( DET(X(N0),Y(N0),Z(N0),X(NF),Y(NF),Z(NF),         XP,YP,ZP) .LT. 0. ) GO TO 6
  
  IF (STORE(ABS(X(N0)*XP + Y(N0)*YP + Z(N0)*ZP))   .LT. 1.0-4.0*EPS) THEN
5    IF ( DET(X(N1),Y(N1),Z(N1),X(N0),Y(N0),Z(N0),           XP,YP,ZP) .GE. 0. ) THEN
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
        IF (N1 .EQ. NL) GO TO 14
        GO TO 5
     ENDIF
  ENDIF
  
  N0 = N1
  GO TO 2
  
6 N2 = NF
  
7 N3 = N0
  N1S = N1
  N2S = N2
  
8 B3 = DET(X(N1),Y(N1),Z(N1),X(N2),Y(N2),Z(N2),XP,YP,ZP)
  IF (B3 .LT. 0.) THEN
     
     LP = LSTPTR(LEND(N2),N1,LIST,LPTR)
     IF (LIST(LP) .LT. 0) GO TO 9
     LP = LPTR(LP)
     N4 = ABS(LIST(LP))
     
     IF ( DET(X(N0),Y(N0),Z(N0),X(N4),Y(N4),Z(N4),           XP,YP,ZP) .LT. 0. ) THEN
        N3 = N2
        N2 = N4
        N1S = N1
        IF (N2 .NE. N2S  .AND.  N2 .NE. N0) GO TO 8
     ELSE
        N3 = N1
        N1 = N4
        N2S = N2
        IF (N1 .NE. N1S  .AND.  N1 .NE. N0) GO TO 8
     ENDIF
     
     N0 = JRAND(N, IX,IY,IZ )
     GO TO 2
  ENDIF
  
  IF (B3 .GE. EPS) THEN
     
     B1 = DET(X(N2),Y(N2),Z(N2),X(N3),Y(N3),Z(N3),           XP,YP,ZP)
     B2 = DET(X(N3),Y(N3),Z(N3),X(N1),Y(N1),Z(N1),           XP,YP,ZP)
     IF (B1 .LT. -TOL  .OR.  B2 .LT. -TOL) THEN
        
        N0 = JRAND(N, IX,IY,IZ )
        GO TO 2
     ENDIF
  ELSE
     
     B3 = 0.
     S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
     PTN1 = XP*X(N1) + YP*Y(N1) + ZP*Z(N1)
     PTN2 = XP*X(N2) + YP*Y(N2) + ZP*Z(N2)
     B1 = PTN1 - S12*PTN2
     B2 = PTN2 - S12*PTN1
     IF (B1 .LT. -TOL  .OR.  B2 .LT. -TOL) THEN
        
        N0 = JRAND(N, IX,IY,IZ )
        GO TO 2
     ENDIF
  ENDIF
  
  I1 = N1
  I2 = N2
  I3 = N3
  IF (B1 .LT. 0.0) B1 = 0.0
  IF (B2 .LT. 0.0) B2 = 0.0
  
  RETURN
  
9 N1S = N1
  N2S = N2
  NL = 0
  
10 LP = LEND(N2)
  LP = LPTR(LP)
  NEXT = LIST(LP)
  IF ( DET(X(N2),Y(N2),Z(N2),X(NEXT),Y(NEXT),Z(NEXT),         XP,YP,ZP) .GE. 0. ) THEN
     
     S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
     Q(1) = X(N1) - S12*X(N2)
     Q(2) = Y(N1) - S12*Y(N2)
     Q(3) = Z(N1) - S12*Z(N2)
     IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 11
     IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3)     .GE. 0.) GO TO 11
     
     NL = N2
  ENDIF
  
  N1 = N2
  N2 = NEXT
  IF (N2 .NE. N1S) GO TO 10
  
  I1 = N1S
  I2 = N1S
  I3 = 0
  
  RETURN
  
11 NF = N2
  IF (NL .EQ. 0) THEN
     
     N2 = N2S
     N1 = N1S
     
12   LP = LEND(N1)
     NEXT = -LIST(LP)
     IF ( DET(X(NEXT),Y(NEXT),Z(NEXT),X(N1),Y(N1),Z(N1),           XP,YP,ZP) .GE. 0. ) THEN
        
        S12 = X(N1)*X(N2) + Y(N1)*Y(N2) + Z(N1)*Z(N2)
        Q(1) = X(N2) - S12*X(N1)
        Q(2) = Y(N2) - S12*Y(N1)
        Q(3) = Z(N2) - S12*Z(N1)
        IF (XP*Q(1) + YP*Q(2) + ZP*Q(3) .GE. 0.) GO TO 13
        IF (X(NEXT)*Q(1) + Y(NEXT)*Q(2) + Z(NEXT)*Q(3)      .GE. 0.) GO TO 13
        
        NF = N1
     ENDIF
     
     N2 = N1
     N1 = NEXT
     IF (N1 .NE. N1S) GO TO 12
     
     I1 = N1
     I2 = N1
     I3 = 0
     
     RETURN
     
13   NL = N1
  ENDIF
  
  I1 = NF
  I2 = NL
  I3 = 0
  
  RETURN
  
14 I1 = 0
  I2 = 0
  I3 = 0
  
  RETURN
END SUBROUTINE TRFIND



INTEGER FUNCTION JRAND (N, IX,IY,IZ )
  INTEGER N, IX, IY, IZ
  REAL U,X
  IX = MOD(171*IX,30269)
  IY = MOD(172*IY,30307)
  IZ = MOD(170*IZ,30323)
  X = (REAL(IX)/30269.) + (REAL(IY)/30307.) +    (REAL(IZ)/30323.)
  U = X - INT(X)
  JRAND = REAL(N)*U + 1.
  RETURN
END FUNCTION JRAND

INTEGER FUNCTION LSTPTR (LPL,NB,LIST,LPTR)
  INTEGER LPL, NB, LIST(*), LPTR(*)
  
  INTEGER LP, ND
  
  LP = LPTR(LPL)
1 ND = LIST(LP)
  IF (ND .EQ. NB) GO TO 2
  LP = LPTR(LP)
  IF (LP .NE. LPL) GO TO 1
  
2 LSTPTR = LP
  RETURN
END FUNCTION LSTPTR

REAL FUNCTION STORE (X)
  REAL X
  
  REAL Y
  COMMON/STCOM/Y
  Y = X
  STORE = Y
  RETURN
END FUNCTION STORE


subroutine geocords(theta_geo, phi_geo, theta_source, phi_source, gx, gy, gz, sl)
  
  implicit none
  real, intent(in)  :: theta_geo, phi_geo, theta_source, phi_source
  real, intent(out) :: gx(3), gy(3), gz(3), sl(3)
  real gyr,gxr,gzr,slr
  real*4 dtr
  dtr=acos(-1.0)/180.0
  gz(1)=sin(theta_geo*dtr)*cos(phi_geo*dtr)
  gz(2)=sin(theta_geo*dtr)*sin(phi_geo*dtr)
  gz(3)=cos(theta_geo*dtr)
  
  gzr=sqrt(gz(1)*gz(1)+gz(2)*gz(2)+gz(3)*gz(3))
  gz(1)=gz(1)/gzr
  gz(2)=gz(2)/gzr
  gz(3)=gz(3)/gzr
  
  
  sl(1)=sin(theta_source*dtr)*cos(phi_source*dtr)
  sl(2)=sin(theta_source*dtr)*sin(phi_source*dtr)
  sl(3)=cos(theta_source*dtr)
  
  slr=sqrt(sl(1)*sl(1)+sl(2)*sl(2)+sl(3)*sl(3))
  sl(1)=sl(1)/slr
  sl(2)=sl(2)/slr
  sl(3)=sl(3)/slr
  
  
  gy(1)=gz(2)*sl(3)-gz(3)*sl(2)
  gy(2)=gz(3)*sl(1)-gz(1)*sl(3)
  gy(3)=gz(1)*sl(2)-gz(2)*sl(1)
  
  gyr=sqrt(gy(1)*gy(1)+gy(2)*gy(2)+gy(3)*gy(3))
  
  gy(1)=gy(1)/gyr
  gy(2)=gy(2)/gyr
  gy(3)=gy(3)/gyr
  
  
  gx(1)=gy(2)*gz(3)-gy(3)*gz(2)
  gx(2)=gy(3)*gz(1)-gy(1)*gz(3)
  gx(3)=gy(1)*gz(2)-gy(2)*gz(1)
  
  gxr=sqrt(gx(1)*gx(1)+gx(2)*gx(2)+gx(3)*gx(3))
  gx(1)=gx(1)/gxr
  gx(2)=gx(2)/gxr
  gx(3)=gx(3)/gxr
  
  
end subroutine geocords




subroutine geo_to_space(theta_u, phi_u, gx, gy, gz,      dirx, diry, dirz, az, el)
  
  implicit none
  real, intent(in)  :: theta_u, phi_u, gx(3), gy(3), gz(3)
  real, intent(out) :: dirx, diry, dirz, az, el
  real r, xg, yg, zg
  real*8 dtr
  dtr=acos(-1.d0)/180.d0
  xg=sin(theta_u*dtr)*cos(phi_u*dtr)
  yg=sin(theta_u*dtr)*sin(phi_u*dtr)
  zg=cos(theta_u*dtr)
  
  dirx=xg*gx(1)+yg*gy(1)+zg*gz(1)
  diry=xg*gx(2)+yg*gy(2)+zg*gz(2)
  dirz=xg*gx(3)+yg*gy(3)+zg*gz(3)
  
  r = sqrt(dirx*dirx+diry*diry+dirz*dirz)
  dirx=dirx/r
  diry=diry/r
  dirz=dirz/r
  
  az = atan2(diry, dirx)/dtr
  if (az.lt.0.0) az = az +360.0
  el = 90.0 - acos(dirz)/dtr
end subroutine geo_to_space



subroutine calc_sphere_dist(ra1,dec1,ra2,dec2,dist)
  
  implicit none
  real, intent(in)  ::  ra1,dec1,ra2,dec2
  real, intent(out) ::  dist
  real dtr,x,y
  
  
  dtr=acos(-1.0)/180.0
  y=sqrt(     (cos(dec2*dtr)*sin((ra1-ra2)*dtr))**2+     (cos(dec1*dtr)*sin(dec2*dtr)-sin(dec1*dtr)*cos(dec2*dtr)*     cos((ra1-ra2)*dtr))**2)
  x=sin(dec1*dtr)*sin(dec2*dtr)+cos(dec1*dtr)*cos(dec2*dtr)*     cos((ra1-ra2)*dtr)
  dist=atan2(y,x)/dtr
  return
end subroutine calc_sphere_dist


subroutine sum_at_scat(direct_diff_matrix,     at_scat_data_lo,at_scat_data_hi,l_fract,     ienerg,nobins_out,out_matrix)
  
  
  
  integer, intent(in):: ienerg,nobins_out
  integer:: nth
  
  !      real,dimension(ienerg,nobins_out), intent(in) :: last_matrix
  real, dimension(ienerg,nobins_out), intent(in) ::     direct_diff_matrix
  real, dimension(ienerg,ienerg)    , intent(in) :: at_scat_data_lo
  real, dimension(ienerg,ienerg)    , intent(in) :: at_scat_data_hi
  
  real l_fract
  !      real coslat_corr
  
  integer i,j,k
  
  real,dimension(ienerg,nobins_out), intent(out) :: out_matrix
  
  
  
  
  
  
  
  
  do i=1,ienerg
     do j=1,nobins_out
        do k=1,ienerg
           out_matrix(i,j) =out_matrix(i,j)+(at_scat_data_lo(i,k)*              l_fract+at_scat_data_hi(i,k)*         (1.0-l_fract))*direct_diff_matrix(k,j)
           
        enddo
     enddo
  enddo
  
  
  
end subroutine sum_at_scat


subroutine highres_ephoton_interpolator (ebin_edge_in, nobins_in, ein, nvbins, matrix, edif_edge_lo, edif_edge_hi, nhbins, new_epx_lo, new_epx_hi, diff_matrix)
  implicit none
  real, dimension(nobins_in), intent(in) :: ebin_edge_in
  integer, intent(in) :: nobins_in
  real, dimension(nvbins+1), intent(in) :: ein
  integer, intent(in) :: nvbins
  real, dimension(70,64), intent(in) :: matrix
  real, dimension(70,64), intent(in) :: edif_edge_lo, edif_edge_hi
  integer, intent(in) ::nhbins
  real, dimension(nobins_in,64), intent(out) ::  new_epx_lo, new_epx_hi,diff_matrix
  
  
  
  integer i,j,k, ivfind
  
  real mu


  ivfind=1

  
  do i=1,nobins_in  
     do j=ivfind, nvbins
        if((ebin_edge_in(i).ge.ein(j)).and.(ebin_edge_in(i).lt.ein(j+1))) then
           ivfind=j
           
           mu=(log(ebin_edge_in(i))-log(ein(ivfind)))/(log(ein(ivfind+1))-log(ein(ivfind)))
           if ((mu.lt.0.0).and.(mu.gt.-1e-5)) mu=0.0
           if ((mu.gt.1.0).and.(mu.lt.1.00001)) mu=1.0

           
           
           do k=1,nhbins
              new_epx_lo(i,k)=(edif_edge_lo(ivfind,k)/ein(ivfind)*(1-mu)+edif_edge_lo(ivfind+1,k)/ein(ivfind+1)*mu)*ebin_edge_in(i)
  
              new_epx_hi(i,k)=(edif_edge_hi(ivfind,k)/ein(ivfind)*(1-mu)+edif_edge_hi(ivfind+1,k)/ein(ivfind+1)*mu)*ebin_edge_in(i)
              diff_matrix(i,k)=(matrix(ivfind,k)*(1-mu)+ matrix(ivfind+1,k)*mu)
           enddo
           
  
        endif
     enddo
  enddo
  
  return
end subroutine highres_ephoton_interpolator



subroutine atscat_highres_ephoton_interpolator (ebin_edge_in, nobins_in, e_in, nvbins, matrix, nobins_out, new_matrix)
  
  implicit none
  
  integer, intent(in) :: nobins_in,nvbins,nobins_out
  real, dimension(nobins_in+1), intent(in) :: ebin_edge_in !Desired input bin edges
  real, dimension(nvbins+1), intent(in) :: e_in !original input bin edges
  real, dimension(70,nobins_out), intent(in) :: matrix !original matrix
  
  real, dimension(nobins_in,nobins_out), intent(out) :: new_matrix !interpolated matrix
  
  integer i,j,k, ivfind
  real mu

  do i=1, nobins_in
     do j=1, nobins_out
        new_matrix(i,j) = 0.0
     enddo
  enddo
  ivfind=1
  do i=1,nobins_in  
     do j=ivfind, nvbins
        if((ebin_edge_in(i).ge.e_in(j)).and.(ebin_edge_in(i).lt.e_in(j+1))) then
           ivfind=j
           mu=(log(ebin_edge_in(i))-log(e_in(ivfind)))/(log(e_in(ivfind+1))-log(e_in(ivfind)))
           if ((mu.lt.0.0).and.(mu.gt.-1e-5)) mu=0.0
           if ((mu.gt.1.0).and.(mu.lt.1.00001)) mu=1.0

           
           do k=1,nobins_out
              new_matrix(i,k)=(matrix(ivfind,k)*(1-mu)+ matrix(ivfind+1,k)*mu)
           enddo
           
        endif
     enddo
  enddo
  return
end subroutine atscat_highres_ephoton_interpolator


subroutine echan_integrator (nobins_in, diff_matrix, edif_edge_lo, edif_edge_hi, nhbins,&
     ebin_edge_out, nobins_out, binned_matrix)
  implicit none
  
  
  integer, intent(in) ::  nobins_in, nobins_out, nhbins
  real, dimension(nobins_in,64), intent(in)  :: diff_matrix
  real, dimension(nobins_in,64), intent(in) :: edif_edge_lo, edif_edge_hi
  real, dimension(nobins_out+1), intent(in) :: ebin_edge_out
  
  
  real, dimension(nobins_in,nobins_out), intent(out) ::  binned_matrix
  
  
  
  integer jcdif,ivh,ihlow,ihhigh,ihbin,ihlfind,ihfind,ivhsum,ihover,nhpoints,icbin,icdif,i,j
  real diff_matrix_vec(nhbins),edif_edgeh(nhbins+1), edif_cent(nhbins),row_tot(nobins_out+1),hlow,hhigh,hwide, hchunk,euse,row_entry
  
  real pi, mu, mu2
  
  pi = 3.1415926535897931
  
  
  

  
  

  do i=1,nobins_out+1
     row_tot(i)=0.
  enddo
  do i=1, nobins_in
     do j=1,nobins_out
        binned_matrix(i,j)=0.
     enddo
  enddo
  
  !loop over input (photon) energies   	  
  do jcdif=1,nobins_in 
     !extract row and divide by bin width to make cm^2/kev units.	  
     do ivh=1,nhbins
        diff_matrix_vec(ivh) = diff_matrix(jcdif,ivh)/(edif_edge_hi(jcdif,ivh)-edif_edge_lo(jcdif,ivh))
        edif_edgeh(ivh)=edif_edge_hi(jcdif,ivh)
        edif_cent(ivh)=     (edif_edge_lo(jcdif,ivh)+edif_edge_hi(jcdif,ivh))/2.
        
  
     enddo !do ivh=1,nhbins
     edif_edgeh(nhbins+1)=edif_edge_hi(jcdif,nhbins)+(edif_edge_hi(jcdif,nhbins)-edif_edge_hi(jcdif,nhbins-1))
     
  
     ihlow=0
     ihhigh=0
     
     !ready for  the horizontal summation
     !---------------- horizontal loop -----------------
     
     do ihbin=1,nobins_out ! loop over output bins
        
          hlow=ebin_edge_out(ihbin)
          hhigh=ebin_edge_out(ihbin+1)
          hwide=hhigh-hlow
          if(ihlow.eq.0) then 
             !next we identify the compressed matrix horizontal bin indices
             !hlow and ihhigh. this is the first horizontal bin so search for ihlow
             
             do ihlfind=1,nhbins-1 
                if((hlow.gt.edif_cent(ihlfind)).and.(hlow.le.edif_cent(ihlfind+1))) then 
                ihlow=ihlfind
             goto 150
              endif !if( (hlow.gt.edif_cent(ihlfind)).and.
           enddo ! do ihlfind=1,nhbins-1
150        continue
!             ihlfind =1
!             do while ((ihlow.eq.0).and.ihlfind.lt.nhbins)
!                if (hlow.gt.edif_cent(ihlfind)).and.(hlow.le.edif_cent(ihlfind+1)) then
 !               ihlow=ihlfind
             
 !               ihlfind = ihlfind +1
 !            enddo
                

                
        endif !if(ihlow.eq.0) then   
        
        ! if ihlow is still zero after this tremedous effort to 
        !non zero it then hlow lt edif_cent(1). assume ihhigh on scale
        !next we are going to break the search into three sections for the 
        !horizontal bin calculation. the first section is for binned matrix 
        !elements that overflow the bottom of the compressed horizontal matrix
!data ( delinieatied by (c111111111111111111111111111) commnets).
        !the second section is for overflowing the top end of the horizontal
        !matrix ( delinieatied by (c22222222222222222222222222) comments)
        !the third section is for horizontal rows contained by the comressed matrix energies.
        !( delinieatied by (c333333333333333333333333333333) comments) 
        
        !secton 1 for overflowing the bottom of the hrizontal matrix
        !111111111111111111111111111111111111111111111111111111111111111
        !111111111111111111111111111111111111111111111111111111111111111
        
          if(hlow.le.edif_cent(1)) then 
             if( hhigh.gt.edif_cent(1)) then 
                !locate ihhigh
                do ihfind=1,nhbins-1 
                   if((hhigh.gt.edif_cent(ihfind)).and.(hhigh.le.edif_cent(ihfind+1))) then
                      ihhigh=ihfind
                   endif !if( (hhigh.gt.edif_cent(ihfind)) .and.
                enddo !do ihfind=1,nhbins-1
                !               if(hhigh.gt.edif_cent(nhbins) ) ihhigh=nhbins-1
                nhpoints=ihhigh+2
                hchunk=hwide/float(nhpoints-1)
                do icbin=1,nhpoints
                   euse=hlow+hchunk*float(icbin-1)
                   if(euse.le.edif_cent(1)) then
                      row_entry=diff_matrix_vec(1)*euse/edif_cent(1)
                      !the  *euse/edif_cent(1) results in tapering to 0.0 below lowest
                      !horizontal compressed point
                      goto 1500
                   endif ! if(euse.le.edif_cent(1)) then
                   if(euse.gt.edif_cent(1)) then 
                      do icdif=1,ihhigh
                         if((euse.gt.edif_cent(icdif)).and.(euse.le.edif_cent(icdif+1))) then 
                            row_entry=diff_matrix_vec(icdif) + (diff_matrix_vec(icdif+1)-diff_matrix_vec(icdif))*&
                                 (euse-edif_cent(icdif))/(edif_cent(icdif+1)-edif_cent(icdif))
                            !standard interpolation betweeen two points in horizontal energy
                            goto 1500
                         endif  ! if( (euse.gt.edif_cent(icdif)) .and.
                      enddo !do icdif=1,ihhigh
                   endif ! if(euse.gt.edif_cent(1)) 
1500               continue
                   row_tot(ihbin)=row_entry+row_tot(ihbin) ! sum up horizontal 
                   !contributions each with the same weighting defined by the width hchunk
                   row_entry=0.0 ! reset this each time
                enddo ! end of icbin loop
                row_tot(ihbin)=row_tot(ihbin)*hwide/float(nhpoints)
                !hwide = horizontal bin width 
                !nhpoints = # of samples used
                !convert from counts/(unit energy) to counts/bin
             endif ! hhigh gt edif_cent(1)
             if(hhigh.le.edif_cent(1)) then 
                row_tot(ihbin)=diff_matrix_vec(1)*((hlow+hhigh)/2.)/edif_cent(1)*hwide
                !here just load this one element extrapolated downward from edif_cent(1)
             endif ! hhigh le edif_cent(1)
          endif ! end of low hlow loop
          !end secton 1 for overflowing the bottom of the hrizontal matrix
          !111111111111111111111111111111111111111111111111111111111111111
          !111111111111111111111111111111111111111111111111111111111111111
          !section 2 for overflowwing the top of the horizontal compressed matrix
          !22222222222222222222222222222222222222222222222222222222222222222
          !22222222222222222222222222222222222222222222222222222222222222222
          !ihlow= last lower edif_cent index: below previous edge
          
          if(ihlow.ge.nhbins) then 
             if(hlow.gt.edif_edgeh(nhbins+1)) then 
                row_tot(ihbin)=-1.0 ! not used in sum below (check 700 continue statment)
                ihover=ihbin
                goto 600
             endif ! if(hlow.gt.edif_edgeh(65))
             if(hlow.le.edif_edgeh(nhbins+1)) then 
                !could still be below this edge
                if(hhigh.le.edif_edgeh(nhbins+1)) then 
                   !bothe bin edges greater than edif_cent(64) and
                   !less than edif_edgeh(65)
                   row_tot(ihbin)= diff_matrix_vec(nhbins)*(edif_edgeh(nhbins+1)-(hlow+hhigh)/2.0)/&
                        (edif_edgeh(nhbins+1)-edif_cent(nhbins))*hwide
                   !interpolation to bin center with hwide providing counts/bin correction
                   goto 600
                endif !if(hhigh.le.edif_edgeh(65)) 
                if(hhigh.gt.edif_edgeh(nhbins+1)) then 
                   !one edge less, one edge greater
                   row_tot(ihbin)=((edif_edgeh(nhbins+1)-hlow)**2 )*diff_matrix_vec(nhbins)/&
                        (2.*(edif_edgeh(nhbins+1)-edif_cent(nhbins)))
                   !interpolation with truncated bin width since the bin exstends outside the
                   !horizontal range
                   goto 600
                endif ! if(hhigh.gt.edif_edgeh(65))
             endif ! if(hhigh.gt.edif_edgeh(65))
          endif  ! if(ihlow.ge.64) then
          
          !end section 2 for overflowwing the top of the horizontal compressed matrix
          !22222222222222222222222222222222222222222222222222222222222222222
          !22222222222222222222222222222222222222222222222222222222222222222
          
          !section 3 for horizontal elements contained by the compressed matrix
          !33333333333333333333333333333333333333333333333333333333333333333
          !33333333333333333333333333333333333333333333333333333333333333333
          
          
          if((ihlow.lt.nhbins).and.(ihlow.ge.1)) then 
             if(hhigh.gt.edif_edgeh(nhbins+1)) then 
                hwide=edif_edgeh(nhbins+1)-hlow ! total width adjusted for active response range
                nhpoints=nhbins-ihlow+2
                hchunk=hwide/float(nhpoints-1)
                do icbin=1,nhpoints 
                   euse=hlow+hchunk*float(icbin-1)
                   do icdif=ihlow,nhbins-1 
                      if((euse.gt.edif_cent(icdif)).and.(euse.le.edif_cent(icdif+1))) then 
                         mu = (euse-edif_cent(icdif))/(edif_cent(icdif+1)-edif_cent(icdif))
                         mu2 = (1-cos(mu*pi))/2.0
                         row_entry = diff_matrix_vec(icdif)+(diff_matrix_vec(icdif+1)-diff_matrix_vec(icdif))*mu2
                         goto 400
                      endif ! if( (euse.gt.edif_cent(icdif)) .and.
                   enddo !do icdif=ihlow,63   
                   
400                row_tot(ihbin)=row_entry+row_tot(ihbin)     
                   row_entry=0.0
                   
                enddo ! do icbin=1,nhpoints
                row_tot(ihbin)=row_tot(ihbin)*hwide/float(nhpoints) !cawh 12-16-2005
                ihlow=nhbins
             endif ! if(hhigh.gt.edif_edgeh(65)  
             
             if(hhigh.le.edif_edgeh(nhbins+1)) then 
                !here we are still inside the compressed row
                !so this is the basi!sumation and hwide is left 
                !unviolated
                
                do ihfind=ihlow,nhbins-1 
                   if((hhigh.gt.edif_cent(ihfind)).and.(hhigh.le.edif_cent(ihfind+1)) ) then
                      !search test
                      ihhigh=ihfind
                   endif ! if( (hhigh.gt.edif_cent(ihfind))                 
                enddo !do ihfind=ihlow,63
                if( hhigh.gt.edif_cent(nhbins)) ihhigh=nhbins
                nhpoints=ihhigh-ihlow+2
                if(nhpoints.lt.9) nhpoints=9
                hchunk=hwide/float(nhpoints-1)
                do icbin=1,nhpoints
                   euse=hlow+hchunk*float(icbin-1)             
                   do icdif=ihlow,ihhigh 
                      if(icdif.le.nhbins-1) then 
                         if((euse.gt.edif_cent(icdif)).and.(euse.le.edif_cent(icdif+1))) then 
                            !note: this samples to edif_cent(ihhigh+!)
                            
                            row_entry=diff_matrix_vec(icdif)+(diff_matrix_vec(icdif+1)-diff_matrix_vec(icdif))*&
                                 (euse-edif_cent(icdif))/(edif_cent(icdif+1)-edif_cent(icdif))
                            
                            !standard interpolation
                            goto 500
                         endif !if( (euse.gt.edif_cent(icdif)) .and. 
                      endif ! if(icdif.le.63) 
                      if(icdif.gt.nhbins-1) then 
                         row_entry=diff_matrix_vec(icdif)*(hhigh-edif_cent(nhbins))/&
                              (edif_cent(nhbins)-edif_cent(nhbins-1))
                         !extrapolation
                         goto 500
                      endif ! if(icdif.gt.63)                    
                   enddo ! do icdif=ihlow,ihhigh
500                row_tot(ihbin)=row_entry+row_tot(ihbin)
                   row_entry=0.0
                enddo ! do icbin=1,nhpoints
                row_tot(ihbin)=row_tot(ihbin)*hwide/float(nhpoints) !cawh 12-16-2005
                ihlow=ihhigh
             endif ! if(hhigh.le.edif_edgeh(65)) 
          endif  ! if((ihlow.lt.64).and.(ihlow.ge.1)) 
          !end section 3 for horizontal elements contained by the compressed matrix
          !33333333333333333333333333333333333333333333333333333333333333333
          !33333333333333333333333333333333333333333333333333333333333333333
600       continue
          if (row_tot(ihbin).eq.-1.0) goto 700
          if(ihbin.eq.nobins_out)  ihover=nobins_out+1
       enddo !do ihbin=1,nobins_out loop over output bins
       
700    continue
       !now in vertical bit only goto ihover which can be lower than nobins_out+1
       !remember to reset row_tot
       do ivhsum=1,ihover-1
          binned_matrix(jcdif,ivhsum)=binned_matrix(jcdif,ivhsum)+row_tot(ivhsum)
          row_tot(ivhsum)=0.0 ! empty for next round
       enddo
       
       row_tot(ihover)=0.0
    enddo ! end of loop over photon energies
    
  end subroutine echan_integrator


	
	
   
