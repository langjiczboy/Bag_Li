  module star_solve_module
  contains

  subroutine star_solve(Pressure_c, M_star, R_star,Lambda,Loven)
  use common_module
  use const_module
  use spline_module
  use rk45_module
  implicit none
  Real *8, Intent (In) :: Pressure_c
  Real *8, Intent (Out) :: M_star, R_star,Lambda,Loven

  integer::n_step,i
  real*8::r_begin,r_end
  real*8::r_max=50d0,r_min=0d0
  real*8::yy
  real*8::yr,rs,M_O_R
  real*8::epsilon_begin,M_begin,epsilon_end
  real*8::v_start(n_var)
  real*8::xx(nstpmx), y(nmax, nstpmx)
  Common /path/xx, y
  !  External derivs

  n_step =10000
  r_begin=(r_max-r_min)/n_step/100000d0
  r_end=r_max
  call splint(EOS_Pressure,EOS_epsilon, D2_epsilon_pressure, n_EOS, Pressure_c, epsilon_begin)  !
      call splint(EOS_Pressure,EOS_epsilon, D2_epsilon_pressure, n_EOS, 0d0, epsilon_end)  !

  M_begin=r_begin**3*epsilon_begin/3.
  v_start(1)=M_begin
  v_start(2)=Pressure_c
  v_start(3)=2d0
  Call rkdumb(v_start, n_var, r_begin, r_end, n_step, derivs)

  Do i = 1, n_step
    If (y(2,i)<=0)  Then
      M_star = y(1, i-1)*cmsun
      R_star = xx(i-1)

      yy = y(3, i-1)
      !yr = yy -  3.0 !quark star
      yr = yy - 4*pi*(R_star)**3*epsilon_end/(M_star )*8.962223634751910d-07

      rs = 2.*cg*M_star*1.98855D0/cc**2
      M_O_R = rs/R_star

      loven = M_O_R**5/20.*(1.-M_O_R)**2*(2.-yr+(yr-1.)*M_O_R)*(M_O_R*(6.-3.*yr+3.*M_O_R*(5.*yr-8.)/2.+M_O_R**2/4.*(26.-22.*yr+M_O_R*(3.*yr-2.)+M_O_R**2*(1.+yr)))+3.*(1.-M_O_R)**2*(2.-yr+(yr-1.)*M_O_R)*dlog(1.-M_O_R))**(-1)

      !lamda = 2./3.*loven*R_star**5/cg/1.d3 ! unit: 10^29m^2*kg*s^2

      lambda = 2./3.*loven*(R_star/(M_star*1.98855D0*cg/cc**2))**5
      Exit
    End If

  End Do
 ! write(*,*) y(2,1:210)
  end   subroutine star_solve

  subroutine derivs(x, y, dydx)
  use spline_module
  use common_module
  use const_module
  Implicit None
  real*8::mass,pressure,epsilon,yy
  Real *8 x, y(*), dydx(*)
  real*8::D_pre_epsilon
  real*8::y_F,y_Q
  mass=y(1)
  pressure=y(2)
  yy=y(3)
  call splint(EOS_pressure, EOS_epsilon, D2_epsilon_pressure, n_EOS, pressure, epsilon)
  call splint(EOS_pressure, EOS_D_pre_epsilon, D2_D_pre_epsilon, n_EOS, pressure, D_pre_epsilon)

  dydx(1)=x**2*epsilon!dM/dr
  dydx(2) = -(epsilon+pressure)*(mass+x**3*pressure)/(crcmg*x**2-2.*x*mass)!dP/dr
  y_F = (crcmg*x-x**3*(epsilon-pressure))/(crcmg*x-2.*mass)
  y_Q = (x*(5.*epsilon+9.*pressure+(epsilon+pressure)/(D_pre_epsilon))-6.*crcmg/x)/(crcmg*x-2.*mass) - 4.*((mass+x**3*pressure)/(crcmg*x*x-2.*x*mass))**2
  dydx(3) = -(yy**2+yy*y_F+x*x*y_Q)/x
  end subroutine derivs
  end module star_solve_module
