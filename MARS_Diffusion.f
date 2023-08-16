      program mars_diff
c  This code is built from routines from Paul Makar's canopy model code
c  See Makar et al. 1999; Stroud et al. 2005 for details

c Set model array sizes: 
c  n:  number of vertical levels in the model
c  nc:  level indicator for near surface data (assumed constant above this level.)
      integer rea,spe,n,nsb,nsteps
      parameter (n=3001,nc=120)
c  Set up model arrays.
c  Species concentration (density)
      real uu(1,1,n)
c  Species concentrations (ppbv, model level)
      real cm(n,1)
c  Model level heights, in m:
      real z(n)
c  Values of K, Tl, sigma-w, kn as a function of height:
      real sigmaw, ustar, nf, ofactor
      real phi, TKE, OL, XX, AX, BX
      real KMethod
      real sigw(n), tl(n), km(n)
c  Values of K only as a function of height: 
      real k(n) !MG: (Tl, sigma-w, and kn are not used)
c  Initial CH4 value (near surface):
      real ch4
c  Work arrays for tridiagonal matrix inversion used in diffusion routine:
      real aa(n),bb(n),cc(n),u(n,1),r(1,n),wrk(n)
      real te1(n),te2(n),rb(n,1),g(1,n)
c Crank-Nicholson Arrays:
      real kkt(1,1,n,2),kkm(1,1,n,2),zzm(1,1,n),zzt(1,1,n)
      real ddtt,flxt(1,1),flxb(1,1),vvd(1,1)
c  Temperature (height):jupyter nto
      real temp(n)
c  Incoming radiation:
      real par
c  csza (degrees):
      real csza
c  Emission rate of CH4 at surface (ppbv/min):
      real ech4
c  Solar zenith angle in radians:S
      real theta
c  Latitude, longitude, time zone, of the site:
      real zla,zlo,tz
c  Time factor tau, and scaling ratio R (See Makar et al., 1998):
      real tr(n),rr(n)
c  Integer year, month, day
      integer iy,im,id
c  Double precision values of solar zenith angle calculations:
      double precision lat,long,zen,ttz
c  Current hour, minute, seconds
      integer hh,mm,ss
c  Pressure value:
      real press
c  Variables for tracking cpu clock time:
      character dte*10,tme*10,zne*10
      integer values(8)
c  sum: 
      real sum 
c-------------------------------------------------------------------------
c
      pi=acos(-1.0)
c
c  Open the data, output files:
      open(unit=20,file='input/input.dat',status='old')
      open(unit=21,file='output/output.out',status='unknown')
      open(unit=22,file='output/vdiff',status='unknown')
      open(unit=37,file='output/errors',status='unknown')
      open(unit=38,file='output/warnings',status='unknown')
      open(unit=91,file='output/massdiff',status='unknown')

c      open(unit=81,file='output/icount',status='unknown')
c      open(unit=82,file='output/timechk',status='unknown')
c      open(unit=83,file='output/solar',status='unknown')
c      open(unit=84,file='output/fluxchk',status='unknown')
c      open(unit=85,file='output/TP.txt',status='unknown')

      rewind(20)
      rewind(22)
      rewind(21)
      rewind(37)
      rewind(38)
      write(6,*) 'Files opened'
c 
c  Read in various parameters for the model startup (vertical levels,
c  initial concentrations, concentrations for constant voc's):
      call initread(n,z,zla,zlo,iy,im,id,tz,ch4)
c  Place the initial concentrations in the main model concentration array.
      do ii=1,nc
        cm(ii,1)=ch4
      end do
      do ii=nc+1,n
        cm(ii,1)=0
      end do
      do ii=1,n 
        cm(ii,1)=0.1  ! BG mixing ratio 
      end do
c 
      nsteps =  48 * 5 ! Set the number of 30-min steps to do here ********
      
c  MAIN LOOP incrementing in TIME STEPS
      write(*,*) 'total nsteps = ',nsteps
      do i=1,nsteps
 
c  Read in measurement data for given time step:
        call readstep(dtime,press,par,n,nc,temp,k, 
     &  sigmaw,ustar,TKE,OL,csza)
c
c  First, reset time to start of time step:
        dtime=dtime-0.5
        if(dtime.ge.24.) then
          idd=id+int(dtime/24.)
          time=24.*(dtime/24.-float(int(dtime/24.)))  !dec. hours
          tt2=time
c  convert time in minutes to hhmm format
          time=float(int(time))*100.+
     #         60.*(time-float(int(time)))
        else
          idd=id
          time=dtime
          tt2=time
          time=float(int(time))*100.+
     #         60.*(time-float(int(time)))
        endif
c  Calculate solar zenith angle for PAR and jvalue
c  attenuation in canopy.
c        long=dble(zlo)
c        lat=dble(zla)
c        ttz=-dble(tz) !new solar routine: time zone <0 if west of Greenwich
c        hh=int(tt2)
c        mm=int( (tt2-float(hh))*60.)
c        ss=0
c        call solzen(iy,im,idd,hh,mm,ss,ttz,lat,long,zen)
c        theta=real(zen)
c  Convert to radians
c        theta = theta / 180. * pi
         theta = acos(csza)
        call date_and_time(dte,tme,zne,values)
c  Call routine to calculate CH4 emissions at surface
        call emit(par,theta,ech4,press,temp)
c***********************************************************
c  operator splitting:  repeat the operators for 30 steps,
c  each one of duration 1 minute.  CH4 deposition/emission is
c  therefore 1 minute, while flanking diffusion steps are
c  0.5 minutes = 30 seconds.
        do ichems=1,32
c          write(6,*) "Hour =", dtime + ichems/60., ", hhmm = ", time
c
c  assign K values to momentum and thermal layers
          do kk=1,n                                                      
            kkm(1,1,kk,1)=k(kk)                                          
            zzm(1,1,kk) = float(kk)
            zzt(1,1,kk) = zzm(1,1,kk) - 0.5 
          end do
c  fluxes     
          flxt(1,1)=1.0E-20  ! No flux through the top layer
          flxb(1,1)= 1*press*6.022E23/1.E6/1013./0.08206/temp(1)/1800 ! flux at bottom layer passed to diffusion solver and used as boundary condition (molec/m2/s)
          vvd(1,1) = 0 ! No deposition to ground
          ddtt=1. ! 1 second time step
c
c  Set top layer to 0        
          cm(n,1)=0
c          
          do nn=1,n
            uu(1,1,nn)=cm(nn,1)*  ! Convert from mixing ratio (ppb) to mass density (molecules/m3)
     #              press*6.022E23/1.E6/1013./0.08206/temp(nn)
          end do
c
          do kk=1,30
            call vdif_b(uu,kkt,kkm,zzt,zzm,vvd,flxt,flxb,1,1,n,ddtt) 
          end do
c
          do nn=1,n
            cm(nn,1)=uu(1,1,nn)
     #              /(press*6.022E23/1.E6/1013./0.08206/temp(nn))
c            cm(n,1) = 0
          end do
c
c  Add surface CH4 emissions to lowest layer
          cm(2,1)=cm(2,1)+ech4 !cm is mixing ratio array (ppbv) starting from second level and adding emission rate at surface
          if(cm(1,1).lt.0.) then
            cm(1,1)=0.
          end if
c  Set top layer to 0        
          cm(n,1)=0
c
          do nn=1,n
            uu(1,1,nn)=cm(nn,1)*  ! Convert from mixing ratio to mass density (ppb to molecules)
     #              press*6.022E23/1.E6/1013./0.08206/temp(nn)
          end do
c          
          do kk=1,30
            call vdif_b(uu,kkt,kkm,zzt,zzm,vvd,flxt,flxb,1,1,n,ddtt)
          end do
          do nn=1,n
            cm(nn,1)=uu(1,1,nn)
c            cm(n,1) = 0
     #              /(press*6.022E23/1.E6/1013./0.08206/temp(nn))
          end do
        end do
        sum = 0
        do ii = 1,n 
          sum = sum + cm(ii,1)
        end do 
c********************************************************
c  End of one operator split model step (30 min).  Write out
c  model results to output file:
        ttot=ttot+dt
c  Output results.
c        write(21,40) float(i)
        write(21,40) (cm(ii,1), ii=1,n)
        write(22,40) (cm(ii,1), ii=1,n)
        write(91,*) sum 
c        write(21,41) cm(1,1),cm(10,1),cm(20,1),cm(60,1),cm(70,1)
      end do
c  End of MAIN LOOP

40    format(3001(1x,1pe10.3))
41    format(5(1pe10.3))
      close(20) 
      close(21)
      close(37)
      close(38)
c      close(81)
c      close(82)
c      close(83)
c      close(84)
c      close(85)

      end
c
c
c     !***********************!
c     !**      vdif_b       **!
c     !***********************!
      SUBROUTINE vdif_b(u,kt,km,zt,zm,vd,flxt,flxb,n1,n2,n,dt)
!
!***************************************************************************
!***************************************************************************
!
!  MODIFICATION 4:  April 2004:
!   Code modified for use by Craig Stroud:  old code option ; now just 
!   Crank-Nicholson. 
!
!  MODIFICATION 3:
!  Author: P.A. Makar, February 2004.
!  Offline tests have revealed that the diffusion time step should
!  be on the order of 10 seconds or less in order to achieve a monotonic
!  mass conservative field with the Crank-Nicholson approach.  A time
!  step of one second will therefore be used here.  In addition,
!  the original diffusion code (probably based on the Laasonen approach)
!  has poor mass conservation tendencies, particularly in the case where
!  the lowest model layer has a "spike" in mass (this is regardless of
!  the time step chosen).  The Crank-Nicholson approach conserves mass
!  regardless of the timestep chosen and the initial conditions, though
!  large time steps (> 10 seconds)  result in a loss of monotonicity of
!  the field.  The Crank-Nicholson approach is therefore
!  recommended ( cr == .true. should be set in the calling code).
!  MODIFICATION 2:  
!  Author: P.A. Makar, December 2003.
!  Original code extensively modified in order to allow calculations
!  with either the original vdif algorithm or the Crank-Nicholson
!  approach.  Part of the modification is the inclusion of vertical
!  diffusion constants for the momentum levels (as the values at the
!  interfaces between layers centered on the thermodynamic levels).
!
!  MODIFICATION 1:
!  Author: P.A. Makar, Jan 18, 1999
!  Local workspace arrays are now declared internally using
!  dynamic allocation, rather than making use of the "stor" 
!  storage array.
!  Original vdif code:
!     Authors:  Janusz Pudykiewicz (AES) Canada
!     Creation date:    FY 1989-90

!***************************************************************************
!***************************************************************************
!  Input variables:
!
       implicit none

       INTEGER  :: n1,n2,n,i1,i2,k,j,nm1,kk

       REAL  :: u(n1,n2,n),kt(n1,n2,n,2),zt(n1,n2,n)
       REAL  :: vd(n1,n2),flxt(n1,n2),flxb(n1,n2),dt
!  top value (z = n) of the following two arrays are not 
!  actually used; just included here for dimension consistency
       REAL  :: km(n1,n2,n,2),zm(n1,n2,n)  
       integer :: istop
       real :: vkm_s, vkp_s, alpha_s, beta_s, gamma_s
       real :: kappa_s, epsil_s
 
!  Local workspace arrays from include file:

! (generated in preprocessor)
       INCLUDE 'vdif_b.cdk'
!             km  -  K at momentum levels (layer interfaces)
!             kt  -  K at thermodynamic levels (layer midpoints)
!             vd  -  deposition velocity 
!             flxt-  flux through upper boundary 
!             flxb-  flux at the surface 
!             zt  -  thermodynamic levels (location of layer midpoints)
!             zm  -  momentum levels (location of layer interfaces)
!             n1  -  number of gridpoints in the x direction
!             n2  -  number of gridpoints in the y direction
!-------------------------------------------------------------------
!
c        do kk=1,n
c          write(90,*) kk, u(1,1,kk), km(1,1,kk,1)
c        end do
c
c        write(90,*) 'vd(1,1) = ', vd(1,1)
c        write(90,*) 'flxt(1,1) = ', flxt(1,1)
c        write(90,*) 'flxb(1,1) = ', flxb(1,1)
c        write(90,*) 'n1 = ', n1
c        write(90,*) 'n2 = ', n2
c        write(90,*) 'n = ', n
c        write(90,*) 'dt = ', dt
c
!  Use Crank-Nicholson method:
!
!__evaluate coefficients of tridiagonal system:
!__interior points (all levels except the top two and bottom two):
!      
        DO k  = 3,n-2
        DO i2 = 1,n2
        DO i1 = 1,n1
    
          vkm_s   = km(i1,i2,k-1,1) / (zt(i1,i2,k) - zt(i1,i2,k-1)) 
          vkp_s   = km(i1,i2,k,1) / (zt(i1,i2,k+1) - zt(i1,i2,k))
          alpha_s = zm(i1,i2,k)-zm(i1,i2,k-1)
          beta_s  = vkm_s + vkp_s
          gamma_s = beta_s / alpha_s
          kappa_s = vkm_s  / alpha_s
          epsil_s = vkp_s  / alpha_s

          a(i1,i2,k) = - 0.5 * dt * kappa_s
   
          b(i1,i2,k) = 1.0 + dt * 0.5 * gamma_s
   
          c(i1,i2,k) = - 0.5 * dt * epsil_s

          r(i1,i2,k) = (1.0 - dt * 0.5 * gamma_s) * u(i1,i2,k) 
     &              +(dt * 0.5 * epsil_s) * u(i1,i2,k+1) 
     &              +(dt * 0.5 * kappa_s) * u(i1,i2,k-1)
        END DO
        END DO
        END DO
!
!__ Lower boundary:  k=1
!
        DO i2 = 1,n2
        DO i1 = 1,n1
          alpha_s = (zm(i1,i2,1)-zt(i1,i2,1))
          beta_s  = vd(i1,i2)*dt/alpha_s
          gamma_s = zt(i1,i2,2)-0.5*(zm(i1,i2,1) + zt(i1,i2,1))

          vkm_s = km(i1,i2,1,1) / gamma_s
       
          a(i1,i2,1) = 0.0
        
          b(i1,i2,1) = 1.0 + 0.5 * dt * vkm_s / alpha_s  
     &                    + 0.5 * dt * vd(i1,i2) / alpha_s
 
          c(i1,i2,1) = - 0.5 * dt * vkm_s / alpha_s

          r(i1,i2,1)=( 1.0 - 0.5*dt* (vkm_s + vd(i1,i2))/alpha_s) 
     &              * u(i1,i2,1) 
     &               + (0.5*dt* vkm_s / alpha_s) * u(i1,i2,2)
     &               + dt * flxb(i1,i2) / alpha_s

        END DO
        END DO
!
!__Next layer up:  k=2  
        DO i2 = 1,n2
        DO i1 = 1,n1

          alpha_s = (zm(i1,i2,2)-zm(i1,i2,1)) 
          beta_s = km(i1,i2,2,1) / (zt(i1,i2,3)-zt(i1,i2,2))
          gamma_s = zt(i1,i2,2)-0.5*(zm(i1,i2,1) + zt(i1,i2,1))

          vkm_s = km(i1,i2,1,1) / gamma_s

          a(i1,i2,2) = - 0.5 * dt * vkm_s / alpha_s

          b(i1,i2,2) = 1.0 + dt * 0.5 * (beta_s+vkm_s)/alpha_s

          c(i1,i2,2) = - 0.5 * dt * beta_s/alpha_s

          r(i1,i2,2) = ( 1.0 - 0.5 * dt *   
     &                (beta_s+vkm_s) / alpha_s ) * u(i1,i2,2)
     &                 + (0.5 * dt * beta_s/alpha_s)*u(i1,i2,3)
     &                 + (0.5 * dt * vkm_s/alpha_s)*u(i1,i2,1)
        END DO
        END DO
!__Upper boundary:  k=N
        DO i2 = 1,n2
        DO i1 = 1,n1

           alpha_s = (zt(i1,i2,n)-zm(i1,i2,n-1))
           gamma_s = 0.5 * (zt(i1,i2,n)+zm(i1,i2,n-1))-zt(i1,i2,n-1)

           vkp_s = km(i1,i2,n-1,1)/gamma_s
 
           a(i1,i2,n) = - 0.5 * dt * vkp_s/alpha_s
 
           b(i1,i2,n) = 1.0 + 0.5 * dt * vkp_s/alpha_s
 
           c(i1,i2,n) = 0.0

           r(i1,i2,n) = (1.0 - 0.5 * dt * vkp_s/alpha_s)* u(i1,i2,n) 
     &                + (0.5 * dt * vkp_s/alpha_s)* u(i1,i2,n-1)
     &                 - dt * flxt(i1,i2)/ alpha_s
        END DO
        END DO
!
!__Next layer down:  k=N-1; some variables from k=N are used here, too:
        DO i2 = 1,n2
        DO i1 = 1,n1
           alpha_s = (zm(i1,i2,n-1)-zm(i1,i2,n-2))
           beta_s  =  km(i1,i2,n-2,1)/(zt(i1,i2,n-1)-zt(i1,i2,n-2))
           gamma_s = 0.5 * (zt(i1,i2,n)+zm(i1,i2,n-1))-zt(i1,i2,n-1)

           vkp_s = km(i1,i2,n-1,1)/gamma_s
           
           a(i1,i2,n-1) = -0.5 * dt * beta_s/alpha_s
     
           b(i1,i2,n-1) = 1.0 + 0.5*dt*(beta_s+vkp_s)/alpha_s          
           c(i1,i2,n-1) = - 0.5 * dt * vkp_s/alpha_s
          
           r(i1,i2,n-1) =  (1.0 - 0.5 * dt 
     &               *   (beta_s+vkp_s)/alpha_s )*u(i1,i2,n-1) 
     &               + (0.5 * dt * vkp_s/alpha_s )*u(i1,i2,n) 
     &               + (0.5 * dt * beta_s/alpha_s )*u(i1,i2,n-2)


        END DO
        END DO
!
!__solve tridiagonal system of equations

        CALL tridiag(a,b,c,r,bet,gam,u,n1,n2,n)
!
      RETURN 
      END
      
      
       SUBROUTINE tridiag (a,b,c,r,bet,gam,u,n1,n2,n)

!
!***************************************************************************
!***************************************************************************
!
!     Authors:        Janusz Pudykiewicz (AES) Canada 
!     Creation date:        FY 1989-90 
!     Current Version:        March 1995
!
!***************************************************************************
!***************************************************************************
!
!
!__solution of the tridiagonal system of equations 
!  Comment by P.A. Makar:  (...set up in vertical
!  diffusion routines 'difvv" and "vdif".
!  A vectorized version of the code found in "Numerical Recipes".
      INTEGER  :: n1,n2,n
      INTEGER  :: j
      REAL  :: a(n1,n2,n),b(n1,n2,n),c(n1,n2,n),r(n1,n2,n)
      REAL  :: u(n1,n2,n)
      REAL  :: bet(n1,n2),gam(n1,n2,n)
      integer :: thdid,thdnum
      integer omp_get_num_threads, omp_get_thread_num
!
!--------------------------------------------------------------------------
!
!  solves for a vector u of length n the tridiagonal linear set:
!
!       b1  c1                        u1      r1
!       a2  b2  c2                 =  u2      r2
!
!
!                           aN  bN    uN      rN
!
!  a,b,c and r are input vectors and are not modified
!
!--------------------------------------------------------------------------
!
!        write(6,*) 'tridiag first loop'
!        call flush(6)
        DO i2 = 1,n2
        DO i1 = 1,n1
       bet(i1,i2) = b(i1,i2,1)
       u(i1,i2,1) = r(i1,i2,1) /  b(i1,i2,1)

        END DO
        END DO
!
        DO j = 2,n
        DO i2 = 1,n2
        DO i1 = 1,n1

        gam(i1,i2,j) = c(i1,i2,j - 1) / bet(i1,i2)

        bet(i1,i2) = b(i1,i2,j) - a(i1,i2,j) * gam(i1,i2,j)

        u(i1,i2,j)=(r(i1,i2,j) - a(i1,i2,j) * u(i1,i2,j - 1))/bet(i1,i2)

        END DO
        END DO
        END DO
!
!__backsubstitution
!
!
        DO j = n - 1,1,-1
        DO i2 = 1,n2
        DO i1 = 1,n1

            u(i1,i2,j) = u(i1,i2,j) - gam(i1,i2,j + 1) * u(i1,i2,j + 1)

        END DO
        END DO
        END DO
      RETURN 
      END
c
c
c     *************************
c     **       solzen        **
c     *************************
      subroutine solzen(iy,im,id,hh,mm,ss,tz,lat,long,zen)
c  A subroutine to calculate the solar zenith angle to an accuracy of
c  about 0.01 degree (ie. timing accuracy of 2.4 seconds)
c       Reference:  J. Meeurs, "Astronomical Algorithms", 1991, 
c       Willman-Bell Inc., Richmond, Va.  429 pp
c
c  Note: Bracketed numbers in the comment lines (eg.  (12.1) )
c  refer to equations in the chapters of the above reference.
c
c    Coding by Paul A. Makar, ARQI, Atmospheric Environment Service,
c    January 30, 1998. 4905 Dufferin Street, Downsview, Ontario, 
c    Canada, M3H 5T4,  email: paul.makar@ec.gc.ca, ph:416-739-4692
c
c*********************************************************************
      implicit none
c  i/o variables
      integer iy,im,id,hh,mm,ss
      double precision lat,long,zen,tz
c  internal variables
      integer i
      double precision zhr,jd,d
      integer m,a,b,iiy
      double precision theta,dt
      double precision t,l0,c,eps0,omega,lp,deleps,eps,ms
      double precision alpha,delta
      double precision rlat, rlong, thet0
      double precision delpsi,corr,houran,altit
      double precision jd0,t0
c
      double precision gmts
      integer mday(12)
      data mday/31,29,31,30,31,30,31,31,30,31,30,31/
c  Leap year conversion:  if mod(iy,4)=0, its a leap year
c  and the number of months in February should be 28:
      if(mod(iy,4).lt.1) mday(2)=28
c  Convert latitude and longitude to radians:
        rlong=long/180.d+00*dacos(-1.d+00)
        rlat=lat/180.d+00*dacos(-1.d+00)
        print*, iy, im, id, hh, mm, ss
c
c  Convert local time to GMT
c
c  (a)  hh mm ss to decimal hours
c     
      zhr=dble(hh)+dble(mm)/60.+dble(ss)/3600.
c
c  (b)  add time zone to get GMT
c
      zhr=zhr-tz
c
c  (c)  correct local day month and year, if necessary
c
      if(zhr.lt.0.0d+00) then
        zhr=24.d+00+zhr
        id=id-1
        if(id.lt.1) then
          if(im.gt.1) then
            im=im-1
            id=mday(im)
          else
            iy=iy-1
            im=12
            id=mday(im)
          end if
        end if
      end if
      if(zhr.ge.24.0d+00) then
        zhr=zhr-24.d+00
        id=id+1
        if(id.gt.mday(im)) then
           im=im+1
           if(im.gt.12) then
             im=1
             iy=iy+1
             id=1
           else
             id=1
           end if
        end if
      end if 
c
c  Correct universal time to dynamical time for years prior to 
c  1989, from Schmadel and Zech formula (pg 74, Meeurs)
c    
      if(iy.ge.1900.and.iy.le.1987) then
c  Calculate time elapsed since 1900.0 in Julian Centuries
       theta=0.d+00
       do i=1900,iy-1
         if(mod(i,4).lt.1) then
            theta=theta+365.d+00
         else
            theta=theta+366.d+00
         end if
       end do
       do i=1,im-1
         theta=theta+dble(mday(i))
       end do
       d=zhr/24.+dble(id)
       theta=theta+d
c convert total days since 1900. to julian centuries:
       theta=theta/365.25d+00/100.d+00
c Use Schmadel and Zech' formula:
       dt=-2.0d-05+2.97d-04*theta+2.5184d-02*theta**2.
     #    -1.81133d-01*theta**3.
     #    +5.553040d-01*theta**4.
     #    -8.61938d-01*theta**5.
     #    +6.77066d-01*theta**6.
     #    -2.12591d-01*theta**7.
c  Dynamical time correction in seconds:
       dt=dt*24.*60.*60.
       end if
       if(iy.eq.1988) dt=55.8d+00
       if(iy.eq.1989) dt=56.4d+00
       if(iy.eq.1990) dt=56.9d+00
       if(iy.gt.1990) dt=0.d+00
c  Convert dynamical time correction to hours
       dt=dt/60./60.
c  Dynamical time = GMT + dt:
       zhr=zhr+dt
c
c  (c)  correct local day month and year, if necessary
c
      if(zhr.lt.0.0d+00) then
        zhr=24.d+00+zhr
        id=id-1
        if(id.lt.1) then
          if(im.gt.1) then
            im=im-1
            id=mday(im)
          else
            iy=iy-1
            im=12
            id=mday(im)
          end if
        end if
      end if
      if(zhr.ge.24.0d+00) then
        zhr=zhr-24.d+00
        id=id+1
        if(id.gt.mday(im)) then
           im=im+1
           if(im.gt.12) then
             im=1
             iy=iy+1
             id=1
           else
             id=1
           end if
        end if
      end if
c
c
c  Calculate the julian ephemiris day of the year (GMT Gregorian
c  calendar assumed as the starting point)  (7.1)
c     
      d=zhr/24.+dble(id) 
      if(im.gt.2) then
         m=im
         iiy=iy
      else
         m=im+12
         iiy=iy-1
      end if
      a=int(dble(iiy)/100.)
      b=2-a+int(dble(a/4.0))
      jd=int(365.25*(dble(iiy)+4716.0))
     #  +int(30.6001*(dble(m)+1.0))
     #  + d + dble(b) -1524.5
c  Calculate the julian ephemiris day of the year (GMT Gregorian
c  calendar assumed as the starting point) for 0 hr GMT for that 
c  date.
       a=int(dble(iiy)/100.)
       b=2-a+int(dble(a/4.0))
       jd0=int(365.25*(dble(iiy)+4716.0))
     #  +int(30.6001*(dble(m)+1.0))
     #  + dble(id) + dble(b) -1524.5
c  
c  Calculate the Time in centuries since the Julian epoch J2000.0 (24.1)
c
      t=(jd-2451545.0)/36525.0
c 
c  Calculate the Time in centuries since the Julian epoch J2000.0 (24.1)
c
      t0=(jd0-2451545.0)/36525.0
c
c  Calculate mean sidereal time at Greenwhich at 0h UT (11.2)
c  in degrees and decimals:
c
      thet0=1.0046061837d+02 + 3.6000770053608d+04*t0
     #+3.87933d-04*t0*t0 - t0*t0*t0/3.871d+07
c  
c  Convert to decimal seconds of arc:
c
      thet0=thet0/360.*24.*60.*60.
c
c  Convert current GMT to seconds:
c   
      gmts=zhr*60.*60.
c  Multiply by constant for sidereal time:
      gmts=gmts*1.00273790935d+00
c
c  Add to thet0 to get sidereal time in seconds of arc:
c
      thet0=thet0+gmts
c  
c  Convert back to degrees and decimals:
c
      thet0=thet0*360./(24.*60.*60.)
c  Convert to radians:
      thet0=thet0/180.d+00*dacos(-1.d+00)
c
c  Calculate the gemetric mean longitude of the sun (24.2)
c
      l0=280.46645 + 36000.76983*t + 3.032d-04*t**2.
      l0=l0/180.d+00*dacos(-1.d+00)
c
c  Calculate mean anomaly of the sun (24.3)
c
      ms=357.52910 + 35999.05030*t - 1.559d-04*t**2. - 4.8d-07*t**3.
c  convert m to radians:
      ms=ms/180.d+00*dacos(-1.d+00)
c
c  Sun's equation of center (pg 152):
c
      c=(1.9146d+00-4.817d-03*t-1.4d-5*t**2.)*dsin(ms)
     # +(1.9993d-02 - 1.01d-04*t)*dsin(2.*ms)
     # + 2.90d-04*dsin(3.*ms)
c  Convert c to radians
      c=c/180.d+00*dacos(-1.d+00)
c
c  Add c and l0 to get true longitude of the sun:
c
      theta=c+l0
c
c  Calculate obliquity of the ecliptic (21.2) :
c
       eps0=21.448d+00-4.68150d+01*t-5.9d-04*t**2.+1.813d-03*t**3.
       eps0=eps0/60.d+00+26.d+00
       eps0=eps0/60.d+00+23.d+00
c  Convert eps0 to radians:
       eps0=eps0/180.d+00*dacos(-1.d+00)
c  
c  Calculate the nutation in the ecliptic:
c
c     Longitude of the ascending node of the Moon's mean orbit
c      (pg 132)
         omega=1.2504452d+02 - 1.934136261d+03*t 
     #         +2.0708d-03*t**2. +t**3./4.5d+05
         omega=omega/180.*dacos(-1.d+00)
c     Mean longitude of the moon:
         lp=(218.3165d+00+4.812678813d+05*t)
         lp=lp/180.d+00*dacos(-1.d+00)
c     Nutation in the ecliptic (seconds of arc):
         deleps=9.20d+00*dcos(omega)+5.7d-01*dcos(2.d+00*l0)
     # +1.d-01*dcos(2.d+00*lp) - 9.d-02*dcos(2.d+00*omega)
c     Convert nutation to radians
         deleps=deleps/60./60.
         deleps=deleps/180.d+00*dacos(-1.d+00)
c     Nutation in the longitude (seconds of arc):
         delpsi=-1.720d+01*dsin(omega) -1.32d+00*dsin(2.d+00*l0)
     # -2.3d-01*dsin(2.d+00*lp)+2.1d-01*dsin(2.d+00*omega)
c     Convert nutation in longitude to radians:
         delpsi=delpsi/60./60.
         delpsi=delpsi/180.d+00*dacos(-1.d+00)
c     Add to mean obliquity to get true obliquity:
         eps=eps0+deleps
c
c     Calculate correction factor for sidereal time (pg. 84)
c     (seconds of arc)
         corr=delpsi*dcos(eps)/15.
c
c     Convert to radians and add to mean sidereal time:
c
         corr=corr/60./60./180.d+00*dacos(-1.d+00)
         thet0=thet0+corr
c  
c     Calculate right ascension of the sun:
c  
       alpha=datan2((dcos(eps)*dsin(theta)),dcos(theta))
c
c     Calculate declination of the sun:
c
       delta=dasin(dsin(eps)*dsin(theta))
c
c     Calculate local hour angle, measured westwards from the South (pg 88)
c
       houran=thet0-rlong-alpha
c
c     Calculate solar altitude angle (12.6):
c
      altit=dasin (dsin(rlat)*dsin(delta) 
     #    + dcos(rlat)*dcos(delta)*dcos(houran) )
c  
c     Convert to degrees and to solar zenith angle:
c
      zen=90.-altit*180.d+00/dacos(-1.d+00)
c      print*, zen
      return
      end
c
c     *************************
c     **       initread      **
c     *************************
      subroutine initread(n,z,zla,zlo,iy,im,id,tz,ch4)
c  Reads in the model heights and initial concentrations
c  for the canopy model
c  Paul Makar, Atmospheric Environment Service, 4905 Dufferin
c  Street, Downsview, Ontario Canada, M3H 5T4.  email:
c  paul.makar@ec.gc.ca, ph: (416)-739-4692
      implicit none
      integer n,i,j,iy,im,id
      real z(n),zla,zlo,tz,ch4
      write(6,*) 'reading initial conditions'
      do i=1,6
        read(20,*)
      end do
      read(20,*) zla,zlo,iy,im,id,tz
      read(20,*)
c      do i=1,3001
c       z(i) = i
c      end do
      read(20,*) (z(j), j=1,n)
c  time varying species:
      write(6,*) 'heights from ', z(1), ' to ', z(n)
      read(20,*)
      read(20,*) ch4
      return
      end
c
c     *************************
c     **       readstep      **
c     *************************
      subroutine readstep(dtime,press,par,n,nc,temp,k,
     #                  sigmaw,ustar,TKE,OL,csza)
c  Reads in canopy model data for one time step:
c  Paul Makar, Atmospheric Environment Service, 4905 Dufferin
c  Street, Downsview, Ontario Canada, M3H 5T4.  email:
c  paul.makar@ec.gc.ca, ph: (416)-739-4692
      implicit none
      integer n,ii,nc
      real dtime,press,par,sigmaw,ustar,TKE,OL,csza
      real temp(n),rh(n),k(n)
c  Read in data time, pressure, par
c      write(6,*) 'reading step values'
      read(20,*)
      read(20,*) dtime,press,par,csza
c      write(6,*) dtime,press,par,csza
c  Read in temperatures for lower nc levels
      read(20,*)
      read(20,*) (temp(ii), ii=1,nc)
c      write(6,*) 'T = ', temp(1)
c  Remaining upper levels constant
      do ii=nc+1,n
        temp(ii)=temp(nc)
      end do
c      write(6,*) 'beginning diffusion'
c  Read in values of diffusion constant for each level:
      read(20,*) 
      read(20,*) (k(ii), ii=1,n)
c      do ii = 1,n
c            read(20,*) k(ii)
c            write(6,*) 'K = ', k(ii), n 
c      end do
      
      write(6,*) 'K = ', k(1)
      return
      end
c
c
c     *************************
c     **        emit         **
c     *************************
      subroutine emit(par,theta,ech4,p,temp)
   
c  calculates CH4 surface emissions
c     
      implicit none
      integer i,j
      real theta,par,p,temp(60)
      real ech4

c  Add code here to calculate ch4 emissions in ppbv/min
      ech4 = 0
      write(90,*) ech4
      return
      end
c