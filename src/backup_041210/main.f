      Program Main

      Integer*8 j1
      Integer*8 k1, k2, k3, k4
      Integer*8 NT, NP
      Integer*8 Mmax, Nmax
      Parameter (NT = 1000, NP = 100000)
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real*8 xvec(NP), yvec(NP), zvec(NP)
      Real*8 h, tau, taumax
      Real*8 qn(6,NP), qnp(6), dqn(6)
c      Real*8 work(18)
      Real*8 I1(100,100), I2(100,100)
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /eigen/ Mmax, Nmax
      Common /intg/ I1, I2
      External dqdt

!      Open(unit=1, file='input/run2.inp', status='old')
      Open(unit=2, file='data/pardyn.txt')

! 100  Format (F17.8)
! 200  Format (2(I10, 3X), 10(F25.12, 3X))

!     General Physical Parameters

      pi = 4.0*atan(1.0)
      c = 2.99792458D8                             !speed of light
      e = 1.6022D-19                               !electron charge
      mk = 9.1094D-31                              !mass of electron
      me = 0.511D6                                 !mass of electron (MeV/c^2)
      eps = 8.8542D-12                             !epsilon
      mu = pi*4.0D-7                               !mu

!     Undulator Parameters

      Eu = 17.5D9                                  !beam energy
      bex = 32.0                                   !
      bey = 32.0                                   !
      emx = 1.4D-6                                 !
      emy = 1.4D-6                                 !
      lb = 25.0D-6                                 !bunch charge
      lu = 35.6D-3                                 !undulator period
      Bmax = 1.0                                   !undulator peak field
      Q = 1.0D-9                                   !beam charge

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
      h = taumax/dfloat(NT)

      Mmax = 50
      Nmax = 50

!     Call the Beam Distribution subroutine
      Call Beamgen(NP,xvec,yvec,zvec)
      Call Beamdist(NP,xvec,yvec)

      print*, "Complete Particle Generation"

c      Call intd2r(Mmax,Nmax,I1,I2)

!     Track the test particles

      print*, "Start tracking test particles"

!     Initialization of particle coordinates

      Do j1 = 1, NP

         qn(1,j1) = xvec(j1)
         qn(2,j1) = yvec(j1)
         qn(3,j1) = zvec(j1)
         qn(4,j1) = 0.0
         qn(5,j1) = 0.0
         qn(6,j1) = 0.0
         tau = 0.0

!     Store initial beam distribution coordiantes
         write(2,*) tau, qn(1,j1), qn(2,j1), xvec(j1), yvec(j1)
c         Write(2,*) tau, qn(1,j1), qn(2,j1), qn(3,j1), 
c     &                   qn(4,j1), qn(5,j1), qn(6,j1)

      End Do

      Stop

      Do 10 k1 = 1, 1 !NT+1

         tau = h*dfloat(k1-1)

         Do 20 k2 = 1, NP

!            tau = h*dfloat(k1-1)
!     RK4 subroutine needs n, h, and tau, but also returns tau and qn0(n)
!     Thus, caution for tau values

            Do k3 = 1, 6
               qnp(k3) = qn(k3,k2)
            End Do

c            print*, "before", k1, k2
!            Call DRKSTP(n, h, tau, qnp, dqdt, work)
c            print*, "after",  k1, k2

c            Write(3,*) tau, qn(1,k2), qn(2,k2) !, qn(3,k2),
c     &                      qn(4,k2), qn(5,k2), qn(6,k2)


            Call dqdt(tau,qnp,dqn)

            Do k4 = 1, 6
               qn(k4,k2) = qnp(k4) + dqn(k4)
            End Do

            xvec(k2) = qn(1,k2)
            yvec(k2) = qn(1,k2)

            print*, k1, k2, qnp(1), qn(1,k2)
!     Store particle's coordinates
c            Write(2,*) tau, qn(1,k2), qn(2,k2) !, qn(3,k2),
c     &                      qn(4,k2), qn(5,k2), qn(6,k2)

 20      Continue

c         Call Beamdist(NP,xvec,vec)

 10   Continue

      Stop
      End


