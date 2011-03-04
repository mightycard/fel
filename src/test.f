      Program Beamgen
cc      Subroutine Beamgen(NP,xvec,yvec,zvec)

      Integer NP, id, it, n1, n2, j1
      Parameter (NP = 1000)
      Real xvec(NP), yvec(NP), zvec(NP), pi
c      Real g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
c      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      pi = 4.0*atan(1.0)

      bex = 32.0                                   !
      bey = 32.0                                   !
      emx = 1.4D-6                                 !
      emy = 1.4D-6   

      sigx = sqrt(bex*emx)
      sigy = sqrt(bey*emy)

      a = 25.0*sigx
      b = 25.0*sigy
      atil = a/pi
      btil = b/pi

!     Random number generation for the beam distribution in the radial direction

      print*, "Generating the particles"

      Call RNORML(xvec,NP)
      Call RNORML(yvec,NP)

!     Longitudinal distribution of test particles

      n1 = 0
      n2 = 0
      NP2 = 1000

c      Call DATIME(id,it)
c      Call RM48IN(id-it,n1,n2)
c      Call RM48(zvec,NP2)

      Do j1 = 1, NP

         xvec(j1) = xvec(j1)*sigx/atil + pi/2.0
         yvec(j1) = yvec(j1)*sigy/btil + pi/2.0
c         zvec(j1) = zvec(j1)*2.0*pi

         print*, j1, xvec(j1)
c         write(3,*) j1, xvec(j1), yvec(j1)

      End Do

      Return
      End
