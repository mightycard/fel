      Subroutine Fields(tau,qn,Esc,Bsc)

      Integer*8 Mmax, Nmax
      Integer*8 Nx, Ny, Nz, Nxp, Nyp
      Integer*8 m, n
      Parameter (Nx = 500, Ny = 500, Nz = 100)
      Parameter (Nxp = 1D3, Nyp = 1D3)
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Real*8 ma(100), ma2(100), nb(100), nb2(100)
      Real*8 tau, zn, ztmp, xn, yn 
      Real*8 sxm, cxm, syn, cyn
      Real*8 sxsy, cxsy, sxcy, cxcy
      Real*8 e0x, e1x, e0y, e1y, e1z, b1y, b1z
      Real*8 qn(6), Esc(3), Bsc(3)
      Real*8 I1(100,100), I2(100,100)
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /eigen/ Mmax, Nmax
      Common /intg/ I1, I2

      c2 = c**2.0
      wc2 = w**2.0/c2

      Do m = 1, Mmax

         ma(m) = m/atil
         ma2(m) = ma(m)**2.0

      End Do

      Do n = 1, Nmax

         nb(n) = n/btil
         nb2(n) = nb(n)**2.0

      End Do

      xn = qn(1)
      yn = qn(2)
      zn = qn(3)

      ztmp = zn + tau
      x0c = x0*cos(ztmp)
      x0s = x0*sin(ztmp)

      e0x = 0.0
      e1x = 0.0
      e0y = 0.0
      e1y = 0.0
      e1z = 0.0
      b1y = 0.0
      b1z = 0.0
           
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

            e0x = e0x - cxsy*I1(m,n)*ma(m)
            e1x = e1x - cxsy*I2(m,n)*(ma2(m)-wc2)

            e0y = e0y - sxcy*I1(m,n)*nb(n)
            e1y = e1y - sxcy*I2(m,n)*ma(m)*nb(n)

            e1z = e1z + sxsy*I2(m,n)*ma(m)

            b1y = b1y - cxsy*I2(m,n)

            b1z = b1z + cxcy*I2(m,n)*nb(n)

         End Do

      End Do

      Esc(1) = e0x + e1x*x0c
      Esc(2) = e0y + e1y*x0c
      Esc(3) =       e1z*x0s*k

      Bsc(1) = 0.0
      Bsc(2) =       b1y*x0c*w*k/c2
      Bsc(3) =       b1z*x0s*w/c2


      Return
      End
