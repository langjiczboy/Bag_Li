module spline_module
  contains 
  Subroutine spline(x, y, n, yp1, ypn, y2)
!***********************************************************************
  Integer n, nmax
  Real *8 :: yp1, ypn, x(n), y(n), y2(n)
  Parameter (nmax=3000)
  Integer :: i, k
  Real *8 :: p, qn, sig, un, u(nmax)
  If (yp1>.99E30) Then
    y2(1) = 0.
    u(1) = 0.
  Else
    y2(1) = -0.5
    u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  End If
  Do i = 2, n - 1
    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p = sig*y2(i-1) + 2.
    y2(i) = (sig-1.)/p
    u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  End Do
  If (ypn>.99E30) Then
    qn = 0.
    un = 0.
  Else
    qn = 0.5
    un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  End If
  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
  Do k = n - 1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  End Do
  Return
  End Subroutine spline
  
Subroutine splint(xa, ya, y2a, n, x, y)
!***********************************************************************
!调用格式，xa ya y2a 原始数据（低精度数据）n,原始数据个数
!x 想要得到的点的x
!y 想要得到的点的y
  Integer n
  Real *8 :: x, y, xa(n), y2a(n), ya(n)
  Integer :: k, khi, klo
  Real *8 :: a, b, h
  klo = 1
  khi = n
  1 If (khi-klo>1) Then
    k = (khi+klo)/2
    If (xa(k)>x) Then
      khi = k
    Else
      klo = k
    End If
    Goto 1
  End If
  h = xa(khi) - xa(klo)
  If (h==0.) Pause 'bad xa input in splint'
  a = (xa(khi)-x)/h
  b = (x-xa(klo))/h
  y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  Return
End Subroutine splint

  end module 