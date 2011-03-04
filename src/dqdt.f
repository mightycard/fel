      Subroutine dqdt(tau, qnp, dqn)

      Real tau, qnp(6), dqn(6)
      Real Esc(3), Bsc(3), Bext(3)
      Real qmw, awp, bwp, bc
      Real pi, c, e, mk, me, eps, mu
      Real Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Real g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil
      Common /phys/ pi, c, e, mk, me, eps, mu
      Common /undulator/ Eu, bex, bey, emx, emy, lb, lu, Bmax, Q
      Common /beam/ g, bg, sigx, sigy, Ql, k, w, x0, a, b, atil, btil

      qmw = e/mk/w
      awp = a*w/pi
      bwp = b*w/pi
      bc = bg/g*c

!     Calculate Space Charge Fields

      Call Fields(tau,qnp,Esc,Bsc)

!     Caclulate External Magnetic Fields

      Bext(1) = 0.0
      Bext(2) = 0.0
      Bext(3) = Bmax*cos(qnp(3)+tau)


!     Differential Equations from Lorentz Force Law

!     dx/dt
      dqn(1) = qnp(4)

!     dy/dt
      dqn(2) = qnp(5)

!     dz/dt
      dqn(3) = qnp(6)

!     dvx/dt
      dqn(4) = (Esc(1) + bwp*qnp(5)*Bsc(3) 
     &                 - bc*qnp(6)*(Bsc(2) + Bext(2)))
     &         *qmw/awp

!     dvy/dt
      dqn(5) = (Esc(2) - awp*qnp(4)*Bsc(3))*qmw/awp

!     dvz/dt
      dqn(6) = (Esc(3) + awp*qnp(4)*(Bsc(2)+Bext(2)))*qmw/bc

      Return
      End

