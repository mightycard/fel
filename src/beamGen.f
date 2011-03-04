      Subroutine BeamGen (NP,xvec,yvec,zvec)

      Integer NP, IJKLIN, NTOTIN, NTO2IN
      Real xvec(NP), yvec(NP), zvec(NP), pi
      Real g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      pi = 4.0*atan(1.0)

!     Random number generation for the beam distribution in the radial direction

      print*, "Generating the particles"

      Call RNORML(xvec,NP)
      Call RNORML(yvec,NP)

!     Random number generation for the beam distribution in the longitudinal direction

      IJKLIN = 239047
      NTOTIN = 1291023
      NTO2IN = 41234798

c      CALL RMARIN(IJKLIN,NTOTIN,NTO2IN)

      Call RANMAR (zvec,NP)

      Return
      End
