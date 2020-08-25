module rk45_module
  contains 
  
  Subroutine rkdumb(vstart, nvar, x1, x2, nstep, derivs)
!***********************************************************************
!定步长四阶Runge-Kutta method.《Visual fortran 常用数值算法集》-何光渝 P605
  use common_module
  Implicit Real *8(A-H, O-Z)
  Integer nstep, nvar
  Real *8 x1, x2, vstart(nvar), xx(nstpmx), y(nmax, nstpmx)
  External derivs
  Common /path/xx, y
!U    USES rk4
  Integer i, k
  Real *8 h, x, dv(nmax), v(nmax)
  Do i = 1, nvar
    v(i) = vstart(i)
    y(i, 1) = v(i)
  End Do
  xx(1) = x1
  x = x1
  h = (x2-x1)/nstep
  Do k = 1, nstep
    Call derivs(x, v, dv)
    Call rk4(v, dv, nvar, x, h, v, derivs)
    If (x+h==x) Pause 'stepsize not significant in rkdumb'
    x = x + h
    xx(k+1) = x
    Do i = 1, nvar
      y(i, k+1) = v(i)
    End Do
  End Do
  Return
End Subroutine rkdumb

!***********************************************************************
Subroutine rk4(y, dydx, n, x, h, yout, derivs)
!***********************************************************************
  use common_module
 
Implicit Real *8(A-H, O-Z)
  Integer n
  Real *8 h, x, dydx(n), y(n), yout(n)
  External derivs
  Integer i
  Real *8 h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)
  hh = h*0.5
  h6 = h/6.
  xh = x + hh

  Do i = 1, n
    yt(i) = y(i) + hh*dydx(i)
  End Do

  Call derivs(xh, yt, dyt)
  Do i = 1, n
    yt(i) = y(i) + hh*dyt(i)
  End Do

  Call derivs(xh, yt, dym)
  Do i = 1, n
    yt(i) = y(i) + h*dym(i)
    dym(i) = dyt(i) + dym(i)
  End Do

  Call derivs(x+h, yt, dyt)
  Do i = 1, n
    yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.*dym(i))
  End Do
  Return
End Subroutine rk4

  end module rk45_module