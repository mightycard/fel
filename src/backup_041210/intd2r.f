      Subroutine intd2r(Mmax,Nmax,I1,I2)

      Integer*8 Mmax, Nmax, Nxp, Nyp
      Integer*8 m, n, j1, j2
      Parameter (Nxp = 1D3, Nyp = 1D3)
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real*8 k02, kmn2, kp2
      Real*8 erfa, erfb, e0
      Real*8 ma(Mmax), ma2(Mmax), nb(Nmax), nb2(Nmax)
      Real*8 xnp, cx, ynp, cy, dxnp, dynp
      Real*8 rhox, rhoy, sxpm, cxpm, sypn
      Real*8 bx0, bxe, by0, bye
      Real*8 I1(Mmax,Nmax), I2(Mmax,Nmax)
      Real*8 Ix1(Mmax), Ix2(Mmax), Iy1(Nmax)
      Real*8 Ix1tmp(Nxp+1), Ix2tmp(Nxp+1), Iy1tmp(Nyp+1)
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      k02 = (1.0 - (bg/g)**2.0)*k**2.0
      cx = (atil/sigx)**2.0/2.0
      cy = (btil/sigy)**2.0/2.0

      erfa = derf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      erfb = derf(a/2.0/sqrt(2.0*pi*sigx**2.0))
      e0 = 2.0*Ql/(pi**3.0*eps*sigx*sigy)/(erfa*erfb)

      dxnp = pi/dfloat(Nxp)
      dynp = pi/dfloat(Nyp)
      bx0 = 0.0
      bxe = pi
      by0 = 0.0
      bye = pi

      Do 11 m = 1, Mmax

         ma(m) = m/atil
         ma2(m) = ma(m)**2.0

         Do 21 j1 = 1, Nxp+1

            xnp = dfloat(j1-1)*dxnp

            rhox = exp(-cx*(xnp-pi/2.0)**2.0)
            sxpm = sin(m*xnp)
            cxpm = cos(m*xnp)

            Ix1tmp(j1) = sxpm*rhox
            Ix2tmp(j1) = cxpm*rhox

 21      Continue

         Ix1(m) = DSIMPS(Ix1tmp,bx0,bxe,Nxp)
         Ix2(m) = DSIMPS(Ix2tmp,bx0,bxe,Nxp)

 11   Continue

      Do 12 n = 1, Nmax

         nb(n) = n/btil
         nb2(n) = nb(n)**2.0

         Do 22 j2 = 1, Nyp+1

            ynp = dfloat(j2-1)*dynp

            rhoy = exp(-cy*(ynp-pi/2.0)**2.0)
            sypn = sin(n*ynp)

            Iy1tmp(j2) = sypn*rhoy

 22      Continue

         Iy1(n) = DSIMPS(Iy1tmp,by0,bye,Nyp)

 12   Continue

      Do 13 m = 1, Mmax

         Do 23 n = 1, Nmax

         kmn2 = ma2(m) + nb2(n)
         kp2 = k02 + kmn2

         I1(m,n) = e0*Ix1(m)*Iy1(n)/kmn2
         I2(m,n) = e0*Ix2(m)*Iy1(n)/kp2

         print*, I1(m,n), I2(m,n)

 23      Continue

 13   Continue

      Return
      End
