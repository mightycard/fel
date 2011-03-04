      Subroutine dqdt(tau, qn, dqn)

      Real*8 tau, qn(6), dqn(6)
      Real*8 Esc(3), Bsc(3), Bext(3)
      Real*8 qmw, awp, bwp, bc
      Real*8 pi, c, e, mk, me, eps, mu
      Real*8 Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real*8 g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      qmw = e/mk/w
      awp = a*w/pi
      bwp = b*w/pi
      bc = bg/g*c

!     Calculate Space Charge Fields

      Call Fields(tau,qn,Esc,Bsc)

!     Caclulate External Magnetic Fields

      Bext(1) = 0.0
      Bext(2) = 0.0
      Bext(3) = Bmax*cos(qn(3)+tau)


!     Differential Equations from Lorentz Force Law

!     dx/dt
      dqn(1) = qn(4)

!     dy/dt
      dqn(2) = qn(5)

!     dz/dt
      dqn(3) = qn(6)

!     dvx/dt
      dqn(4) = (Esc(1) + bwp*qn(5)*Bsc(3) - bc*qn(6)*(Bsc(2) + Bext(2)))
     &         *qmw/awp

!     dvy/dt
      dqn(5) = (Esc(2) - awp*qn(4)*Bsc(3))*qmw/awp

!     dvz/dt
      dqn(6) = (Esc(3) + awp*qn(4)*(Bsc(2)+Bext(2)))*qmw/bc

      Return
      End

