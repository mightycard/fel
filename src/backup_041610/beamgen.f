      Subroutine Beamgen(NP,xvec,yvec,zvec)

      Integer NP, id, it, n1, n2, j1
      Real xvec(NP), yvec(NP), zvec(NP), pi
      Real g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      pi = 4.0*atan(1.0)

!     Random number generation for the beam distribution in the radial direction

      print*, "Generating the particles"

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

      End Do

      Return
      End
