      Subroutine Beamdist(NP,xvec,yvec)

      Integer NP, Npx, Npy
      Parameter (Npx = 100, Npy = 100)
      Integer*8 k1, k2, k3, k4, k5, k6, k7
      Integer*8 countp(Npx,Npy)
      Real xvec(NP), yvec(NP)
      Real*8 pi, xp, xptmp, dxp, yp, yptmp, dyp

      Open(unit=9,file='data/countp.txt')

      print*, ""
      print*, "Counting Particles..."

      pi = 4.0*atan(1.0)

      dxp = pi/dfloat(Npx)
      dyp = pi/dfloat(Npy)
      xptmp = 0.0
      yptmp = 0.0

      Do k1 = 1, Npx

         Do k2 = 1, Npy

            countp(k1,k2) = 0

         End Do

      End Do

      Do k3 = 1, NP

         Do k4 = 1, Npx

            xp = dfloat(k4)*dxp

            Do k5 = 1, Npy

               yp = dfloat(k5)*dyp

               If ((xvec(k3).ge.xptmp).and.(xvec(k3).lt.xp).and.
     &             (yvec(k3).ge.yptmp).and.(yvec(k3).lt.yp)) then

                  countp(k4,k5) = countp(k4,k5) + 1

               End If

               yptmp = yp

            End Do

         xptmp = xp

         End Do

      End Do

      Do k6 = 1, Npx

         xp = dfloat(k6)*dxp

         Do k7 = 1, Npy

            yp = dfloat(k7)*dyp

            Write(9,*) xp, yp, countp(k6,k7)

         End Do

      End Do

      Close(9)

      print*, "done."
      print*, ""
      Return
      End


