      Subroutine Beamgen(NP,xvec,yvec,zvec)

      Integer NP1, id, it, n1, n2
      Integer*8 NP, j1
      Real fspace(200)
      Real Nlx, xvec(NP), xlow, xhigh
      Real Nly, yvec(NP), ylow, yhigh
!     The above variables must be single precisions
      Real*8 zvec(NP)
      Real*8 pi
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      External Nlx, Nly

      pi = 4.0*atan(1.0)

!     Random number generation for the beam distribution in the radial direction

      print*, "Generating the particles"

      xlow = 0.0
      xhigh = pi

      ylow = 0.0
      yhigh = pi

      NP1 = NP

c      Call FUNLXP(Nlx, fspace, xlow, xhigh)
c      Call FUNLUX(fspace, xvec, NP1)

c      Call FUNLXP(Nly, fspace, ylow, yhigh)
c      Call FUNLUX(fspace, yvec, NP1)

      Call RNORML(xvec,NP1)
      Call RNORML(yvec,NP1)

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

      End Do

      Return
      End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!     Function for the initial beam distribution                    !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Real Function Nlx(xhat)

      Real xhat
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real*8 atilp, normx, nx0, nx
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      atilp = atil/sigx/sqrt(2.0)
      normx = sqrt(2.0*pi*sigx**2.0)
      nx0 = derf(a/2.0/normx)
      nx = derf((xhat-a/2.0)/normx)

!      Nlx = (nx+nx0)/nx0/2.0
      Nlx = exp(-atilp**2.0*(xhat-pi/2.0)**2.0)/normx/nx0

      Return
      End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Real Function Nly(yhat)

      Real yhat
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real*8 btilp, normy, ny0, ny
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      btilp = btil/sigy/sqrt(2.0)
      normy = sqrt(2.0*pi*sigy**2.0)
      ny0 = derf(b/2.0/normy)
      ny = derf((yhat-b/2.0)/normy)

!      Nly = (ny+ny0)/ny0/2.0
      Nly = exp(-btilp**2.0*(yhat-pi/2.0)**2.0)/normy/ny0

      Return
      End



