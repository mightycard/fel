      Subroutine Beamdist(NP,xvec,yvec)

      Integer NP
      Integer Npx, Npy
      Parameter (Npx = 50, Npy = 50)
      Integer i1, k1, k2, k3, k4, k5, k6, k7
      Integer countp(Npx,Npy)
      Real xvec(NP), yvec(NP)
      Real pi, xp, xptmp, dxp, yp, yptmp, dyp

      Open(unit=9,file='data/countp.txt')

      print*, ""
      print*, "Counting Particles..."

      pi = 4.0*atan(1.0)

      dxp = pi/float(Npx)
      dyp = pi/float(Npy)
      xptmp = 0.0
      yptmp = 0.0

      Do i1 = 1, NP

c          print*, i1, xvec(i1), yvec(i1)
c          write(3,*) xvec(i1), yvec(i1)

      End Do

      Do k1 = 1, Npx

         Do k2 = 1, Npy

            countp(k1,k2) = 0

         End Do

      End Do

      Do k3 = 1, NP

         Do k4 = 1, Npx

            xp = float(k4)*dxp

            Do k5 = 1, Npy

               yp = float(k5)*dyp

               If ((xvec(k3).ge.xptmp).and.(xvec(k3).lt.xp).and.
     &             (yvec(k3).ge.yptmp).and.(yvec(k3).lt.yp)) then

                  countp(k4,k5) = countp(k4,k5) + 1

               End If

               yptmp = yp

c               print*, k3, k4, k5

            End Do

         xptmp = xp

         End Do

      End Do

      Do k6 = 1, Npx

         xp = float(k6)*dxp

         Do k7 = 1, Npy

            yp = float(k7)*dyp

            Write(9,*) xp, yp, countp(k6,k7)

         End Do

      End Do

      print*, "Done."
      print*, "  "
      Return
      End


