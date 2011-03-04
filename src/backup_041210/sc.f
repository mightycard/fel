!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!     033010 - This program is for calculating electromagnetic      !
!     potentials and fields of the bunch inside of the undulator.   !
!     (C.S. Park)                                                   !          
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Program Fields !(qn,Esc,Bsc)

      Integer*8 Mmax, Nmax, Nx, Ny, Nz, Nxp, Nyp
      Integer*8 m, n, j1, j2, k1, k2, k3, k4, k5
      Parameter (Mmax = 50, Nmax = 50)
      Parameter (Nx = 500, Ny = 500, Nz = 100)
      Parameter (Nxp = 1D3, Nyp = 1D3)
      Real*8 pi, c, mk, me, eps, mu
      Real*8 Eu, g, bg, bex, bey, emx, emy, sigx, sigy
      Real*8 lb, lu, Bmax, Q, Ql, k, w, x0, k02, kmn2, kp2
      Real*8 a, b, atil, btil, dxn, dyn, dxnp, dynp
      Real*8 erfa, erfb, e0
      Real*8 ma(Mmax), ma2(Mmax), nb(Nmax), nb2(Nmax)
      Real*8 tau, zn, dzn, ztmp, xn, xnp, cx, yn, ynp, cy 
      Real*8 rhox, rhoy, sxm, cxm, syn, cyn, sxpm, cxpm, sypn
      Real*8 sxsy, cxsy, sxcy, cxcy, bx0, bxe, by0, bye
      Real*8 I1(Mmax,Nmax), I2(Mmax,Nmax)
      Real*8 Ix1(Mmax), Ix2(Mmax), Iy1(Nmax)
      Real*8 Jx1(Mmax), Jx2(Mmax), Jy1(Nmax)
      Real*8 Ix1tmp(Nxp+1), Ix2tmp(Nxp+1), Iy1tmp(Nyp+1)
      Real*8 ex, ey, ez, by, bz, ex0, ex1, ey0, ey1
      Real*8 e0x(Nx,Ny), e1x(Nx,Ny)
      Real*8 e0y(Nx,NY), e1y(Nx,Ny), e1z(Nx,Ny)
      Real*8 b1y(Nx,Ny), b1z(Nx,NY)
      Parameter (tau = 0.01)

      Open(unit=1, file='data/tau_0.01.dat')

 100  Format (12(F17.8, 3X))

!     General Physical Parameters

      pi = 3.14159265359
      c = 2.99792458D8                             !speed of light
      mk = 9.1094D-31                               !mass of electron
      me = 0.511D6                                  !mass of electron (MeV/c^2)
      eps = 8.8542D-12                             !epsilon
      mu = pi*4.0D-7                               !mu

!     Undulator Parameters

      Eu = 17.5D9
      bex = 32.0
      bey = 32.0
      emx = 1.4D-6
      emy = 1.4D-6
      lb = 25.0D-6
      lu = 35.6D-3
      Bmax = 1.0
      Q = 1.0D-9

      g = Eu/me
      bg = sqrt(g**2.0-1)
      sigx = sqrt(bex*emx)
      sigy = sqrt(bey*emy)
      Ql = Q/(g*lb)
      k = 2.0*pi*g/lu
      w = 2.0*pi*bg*c/lu
      x0 = Q*Bmax/(g**2.0*mk*k**2.0*bg*c)

      a = 50.0*sigx
      b = 50.0*sigy
      atil = a/pi
      btil = b/pi
      k02 = (1.0 - (bg/g)**2.0)*k**2.0
      cx = (atil/sigx)**2.0/2.0
      cy = (btil/sigy)**2.0/2.0
      c2 = c**2.0
      wc2 = w**2.0/c2

      erfa = derf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      erfb = derf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      e0 = 2.0*Ql/(pi**3.0*eps*sigx*sigy)/(erfa*erfb)

      dxn = pi/dfloat(Nx)
      dyn = pi/dfloat(Ny)
      dzn = 2.0*pi/dfloat(Nz)
      dxnp = pi/dfloat(Nxp)
      dynp = pi/dfloat(Nyp)
      bx0 = 0.0
      bxe = pi
      by0 = 0.0
      bye = pi

      print*, g, sigx, sigy
      print*, k, w, x0
      print*, e0, Ql, k02
      print*, x0, x0c
      print*, erfa, erfb, derf(1.0)
      print*, ""

!      stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Integration Starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Integration Part1 Starts


      Do 11 m = 1, Mmax

         ma(m) = m/atil
         ma2(m) = ma(m)**2.0

         Jx1(m) = 0.0
         Jx2(m) = 0.0

         Do 21 j1 = 1, Nxp+1

            xnp = dfloat(j1-1)*dxnp

            rhox = exp(-cx*(xnp-pi/2.0)**2.0)
            sxpm = sin(m*xnp)
            cxpm = cos(m*xnp)

            Jx1(m) = Jx1(m) + sxpm*rhox
            Jx2(m) = Jx2(m) + cxpm*rhox

            Ix1tmp(j1) = sxpm*rhox
            Ix2tmp(j1) = cxpm*rhox

 21      Continue

         Jx1(m) = Jx1(m)*dxnp
         Jx2(m) = Jx2(m)*dxnp

         Ix1(m) = DSIMPS(Ix1tmp,bx0,bxe,Nxp)
         Ix2(m) = DSIMPS(Ix2tmp,bx0,bxe,Nxp)

 11   Continue

!      stop

!     Integration Part1 Ends
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Integration Part2 Starts

      Do 12 n = 1, Nmax

         nb(n) = n/btil
         nb2(n) = nb(n)**2.0

         Jy1(n) = 0.0

         Do 22 j2 = 1, Nyp+1

            ynp = dfloat(j2-1)*dynp

            rhoy = exp(-cy*(ynp-pi/2.0)**2.0)
            sypn = sin(n*ynp)

            Jy1(n) = Jy1(n) + sypn*rhoy

            Iy1tmp(j2) = sypn*rhoy

 22      Continue

         Jy1(n) = Jy1(n)*dynp

         Iy1(n) = DSIMPS(Iy1tmp,by0,bye,Nyp)

 12   Continue

!     Integration Part2 Ends
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Integration Ends
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      Do 13 m = 1, Mmax

         Do 23 n = 1, Nmax

         kmn2 = ma2(m) + nb2(n)
         kp2 = k02 + kmn2

         I1(m,n) = e0*Ix1(m)*Iy1(n)/kmn2
         I2(m,n) = e0*Ix2(m)*Iy1(n)/kp2

 23      Continue

 13   Continue

      print*, "1111111111"

      Do 10 k1 = 1, Nx+1

         xn = dfloat(k1-1)*dxn

!         Do 20 k2 = 1, Ny+1

!            yn = dfloat(k2-1)*dyn

            k2 = 1
            yn = pi/2.0

            e0x(k1,k2) = 0.0
            e1x(k1,k2) = 0.0
            e0y(k1,k2) = 0.0
            e1y(k1,k2) = 0.0
            e1z(k1,k2) = 0.0
            b1y(k1,k2) = 0.0
            b1z(k1,k2) = 0.0
           
            Do m = 1, Mmax

               cxm = cos(m*xn)
               sxm = sin(m*xn)

               Do n = 1, Nmax

                  cyn = cos(n*yn)
                  syn = sin(n*yn)

                  sxsy = sxm*syn
                  cxsy = cxm*syn
                  sxcy = sxm*cyn
                  cxcy = cxm*cyn

                  e0x(k1,k2) = e0x(k1,k2) - cxsy*I1(m,n)*ma(m)
                  e1x(k1,k2) = e1x(k1,k2) - cxsy*I2(m,n)*(ma2(m)-wc2)

                  e0y(k1,k2) = e0y(k1,k2) - sxcy*I1(m,n)*nb(n)
                  e1y(k1,k2) = e1y(k1,k2) - sxcy*I2(m,n)*ma(m)*nb(n)

                  e1z(k1,k2) = e1z(k1,k2) + sxsy*I2(m,n)*ma(m)

                  b1y(k1,k2) = b1y(k1,k2) - cxsy*I2(m,n)

                  b1z(k1,k2) = b1z(k1,k2) + cxcy*I2(m,n)*nb(n)

               End Do

            End Do

!         print*, "2222222222"

! 20      Continue
      
 10   Continue

      Do 30 k3 = 1, Nz+1

         zn = dfloat(k3-1)*dzn

         ztmp = zn + tau
         x0c = x0*cos(ztmp)
         x0s = x0*sin(ztmp)

         Do 40 k4 = 1, Nx+1

            xn = dfloat(k4-1)*dxn

!            Do 50 k5 = 1, Ny+1

!            yn = dfloat(k5-1)*dyn

            k5 = 1
            yn = pi/2.0

            ex = e0x(k4,k5) + e1x(k4,k5)*x0c
            ey = e0y(k4,k5) + e1y(k4,k5)*x0c
            ez =              e1z(k4,k5)*x0s*k

            by =              b1y(k4,k5)*x0c*w*k/c2
            bz =              b1z(k4,k5)*x0s*w/c2

            ex0 = e0x(k4,k5)
            ex1 = e1x(k4,k5)*x0c
            ey0 = e0y(k4,k5)
            ey1 = e1y(k4,k5)*x0c

            Write(1,100) xn, yn, zn, ex, ey, ez, by, bz,
     &                   ex0, ex1, ey0, ey1
            print*, xn, zn, ex

! 50         Continue

 40      Continue

 30   Continue

      Stop
      End
