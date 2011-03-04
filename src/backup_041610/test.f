      Program Main

      Integer j1, id, it, n1, n2
      Integer k1, k2, k3, k4
      Integer NT, NP
      Integer Mmax, Nmax, Nxp, Nyp
      Parameter (NT = 10, NP = 100000)
      Parameter (Mmax = 50, Nmax = 50)
      Parameter (Nxp = 1e3, Nyp = 1e3)
      Real pi, c, e, mk, me, eps, mu
      Real Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real xvec(NP), yvec(NP), zvec(NP)
      Real h, tau, taumax
      Real qn(6,NP), qnp(6), dqn(6)
      Real k02, kmn2, kp2
      Real erfa, erfb, e0
      Real ma(Mmax), ma2(Mmax), nb(Nmax), nb2(Nmax)
      Real xnp, cx, ynp, cy, dxnp, dynp
      Real rhox, rhoy, sxpm, cxpm, sypn
      Real bx0, bxe, by0, bye
      Real I1(Mmax,Nmax), I2(Mmax,Nmax)
      Real Ix1(Mmax), Ix2(Mmax), Iy1(Nmax)
      Real Ix1tmp(Nxp+1), Ix2tmp(Nxp+1), Iy1tmp(Nyp+1)
      Real Esc(3), Bsc(3), Bext(3)
      Real qmw, awp, bwp, bc


!      Open(unit=1, file='input/run2.inp', status='old')
      Open(unit=2, file='data/pardyn.txt')

! 100  Format (F17.8)
! 200  Format (2(I10, 3X), 10(F25.12, 3X))

!     General Physical Parameters

      pi = 4.0*atan(1.0)
      c = 2.99792458e8                             !speed of light
      e = 1.6022e-19                               !electron charge
      mk = 9.1094e-31                              !mass of electron
      me = 0.511e6                                 !mass of electron (MeV/c^2)
      eps = 8.8542e-12                             !epsilon
      mu = pi*4.0e-7                               !mu

!     Undulator Parameters

      Eu = 17.5e9                                  !beam energy
      bex = 32.0                                   !
      bey = 32.0                                   !
      emx = 1.4e-6                                 !
      emy = 1.4e-6                                 !
      lb = 25.0e-6                                 !bunch charge
      lu = 35.6e-3                                 !undulator period
      Bmax = 1.0                                   !undulator peak field
      Q = 1.0e-9                                   !beam charge

      g = Eu/me                                    !gamma
      bg = sqrt(g**2.0-1)                          !beta*gamma
      sigx = sqrt(bex*emx)
      sigy = sqrt(bey*emy)
      Ql = Q/(g*lb)                                !Q_{\lambda}
      k = 2.0*pi*g/lu
      w = 2.0*pi*bg*c/lu
      x0 = Q*Bmax/(g**2.0*mk*k**2.0*bg*c)

      a = 25.0*sigx
      b = 25.0*sigy
      atil = a/pi
      btil = b/pi
     
      taumax = 1.0*2.0*pi
      h = taumax/float(NT)

      qmw = e/mk/w
      awp = a*w/pi
      bwp = b*w/pi
      bc = bg/g*c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      print*, "Test Particle Generation"

!     Transverse distribution of test particles

      Call RNORML(xvec,NP)
      Call RNORML(yvec,NP)

!     Longitudinal distribution of test particles

      n1 = 0
      n2 = 0

      Call DATIME(id,it)
      Call RM48IN(id-it,n1,n2)
      Call RM48(zvec,NP)

      Do j1 = 1, NP

         xvec(j1) = xvec(j1)*sigx/atil + pi/2.0
         yvec(j1) = yvec(j1)*sigy/btil + pi/2.0
         zvec(j1) = zvec(j1)*2.0*pi
c         print*, j1, xvec(j1)
c         write(3,*) j1, xvec(j1), yvec(j1)

      End Do

      print*, "Complete Particle Generation"

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!     Prepare Transverse Integrations

      k02 = (1.0 - (bg/g)**2.0)*k**2.0
      cx = (atil/sigx)**2.0/2.0
      cy = (btil/sigy)**2.0/2.0

      erfa = erf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      erfb = erf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      e0 = 2.0*Ql/(pi**3.0*eps*sigx*sigy)/(erfa*erfb)

      dxnp = pi/float(Nxp)
      dynp = pi/float(Nyp)
      bx0 = 0.0
      bxe = pi
      by0 = 0.0
      bye = pi

      Do 11 m = 1, Mmax

         ma(m) = m/atil
         ma2(m) = ma(m)**2.0

         Do 21 j1 = 1, Nxp+1

            xnp = float(j1-1)*dxnp

            rhox = exp(-cx*(xnp-pi/2.0)**2.0)
            sxpm = sin(m*xnp)
            cxpm = cos(m*xnp)

            Ix1tmp(j1) = sxpm*rhox
            Ix2tmp(j1) = cxpm*rhox

 21      Continue

         Ix1(m) = SIMPS(Ix1tmp,bx0,bxe,Nxp)
         Ix2(m) = SIMPS(Ix2tmp,bx0,bxe,Nxp)

 11   Continue

      Do 12 n = 1, Nmax

         nb(n) = n/btil
         nb2(n) = nb(n)**2.0

         Do 22 j2 = 1, Nyp+1

            ynp = float(j2-1)*dynp

            rhoy = exp(-cy*(ynp-pi/2.0)**2.0)
            sypn = sin(n*ynp)

            Iy1tmp(j2) = sypn*rhoy

 22      Continue

         Iy1(n) = SIMPS(Iy1tmp,by0,bye,Nyp)

 12   Continue

      Do 13 m = 1, Mmax

         Do 23 n = 1, Nmax

         kmn2 = ma2(m) + nb2(n)
         kp2 = k02 + kmn2

         I1(m,n) = e0*Ix1(m)*Iy1(n)/kmn2
         I2(m,n) = e0*Ix2(m)*Iy1(n)/kp2

 23      Continue

 13   Continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     Track the test particles

      print*, "Start tracking test particles"

!     Initialization of particle coordinates

      tau = 0.0
      Do j1 = 1, NP

         qn(1,j1) = xvec(j1)
         qn(2,j1) = yvec(j1)
         qn(3,j1) = zvec(j1)
         qn(4,j1) = 0.0
         qn(5,j1) = 0.0
         qn(6,j1) = 0.0

!     Store initial beam distribution coordiantes
         Write(2,*) tau, qn(1,j1), qn(2,j1)!, qn(3,j1), 
c     &                   qn(4,j1), qn(5,j1), qn(6,j1)

      End Do

      Do 10 k1 = 1, NT+1

         tau = h*float(k1-1)

         Do 20 k2 = 1, NP

            Do k3 = 1, 6
               qnp(k3) = qn(k3,k2)
            End Do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     Calculate Space Charge Fields

c            Call Fields(tau,qnp,Esc,Bsc)

!     Caclulate External Magnetic Fields

            Bext(1) = 0.0
            Bext(2) = 0.0
            Bext(3) = Bmax*cos(qnp(3)+tau)

!     Differential Equations from Lorentz Force Law

!     dx/dt
            dqn(1) = qnp(4)

!     dy/dt
            dqn(2) = qnp(5)

!     dz/dt
            dqn(3) = qnp(6)

!     dvx/dt
            dqn(4) = (Esc(1) + bwp*qnp(5)*Bsc(3)
     &                       - bc*qnp(6)*(Bsc(2) + Bext(2)))
     &               *qmw/awp

!     dvy/dt
            dqn(5) = (Esc(2) - awp*qnp(4)*Bsc(3))*qmw/awp

!     dvz/dt
            dqn(6) = (Esc(3) + awp*qnp(4)*(Bsc(2)+Bext(2)))*qmw/bc

            Do k4 = 1, 6
               qn(k4,k2) = qnp(k4) + dqn(k4)
            End Do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            xvec(k2) = qn(1,k2)
            yvec(k2) = qn(2,k2)

c            print*, k1, k2, qnp(1), qn(1,k2)

!     Store particle's coordinates
            Write(2,*) tau, qn(1,k2), qn(2,k2)!, qn(3,k2),
!     &                      qn(4,k2), qn(5,k2), qn(6,k2)
c            write(2,*) xvec(k2), yvec(k2)

 20      Continue

         print*, k1, tau
c         Call Beamdist(NP,xvec,vec)

 10   Continue

      Stop
      End


