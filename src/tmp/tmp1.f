      Program test 

      Integer j1, j2, mmax, nmax
      Real I1(10,10)

      mmax = 10
      nmax = 10

      call test1(mmax,nmax,I1)

      do j1 = 1, mmax

         do j2 = 1, nmax

            print*, j1, j2, I1(j1,j2)

         end do

      end do

      stop
      end


      subroutine test1(mmax, nmax, I1)

      Integer k1, k2, mmax, nmax
      Real I1(10,10)

      do k1 = 1, mmax

         do k2 = 1, nmax

            I1(k1,k2) = k1*k2

         end do

      end do

      return
      end
