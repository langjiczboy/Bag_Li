  module EOS_initialize_module
  contains
  subroutine EOS_initialize_solve
  use common_module
  use spline_module
  implicit none
  integer :: i
  real*8 ::temp


  Open (123, File=file_read) !'fraction.txt'
  Read (123, *)
  do i = 1, n_EOS
   !   Read (123, *) EOS_rho_B(i),temp, temp,temp,temp,temp,temp,temp,EOS_energy(i), EOS_pressure(i)
   ! Read (123, *) EOS_rho_B(i),EOS_energy(i), EOS_pressure(i)
        Read (123, *) EOS_rho_B(i),EOS_energy(i), EOS_pressure(i),EOS_epsilon(i),EOS_D_pre_rho(i),EOS_D_pre_epsilon(i)

  end do
  Close (123)

  call spline(EOS_pressure, EOS_epsilon, n_EOS,3.D30, 3.D30, D2_epsilon_pressure)
  Call spline(EOS_pressure, EOS_rho_B, n_EOS, 3.D30, 3.D30, D2_rho_B_pressure)
    Call spline(EOS_pressure, EOS_D_pre_epsilon, n_EOS, 3.D30, 3.D30, D2_D_pre_epsilon)

  !   Call spline(EOS_rho_B, EOS_epsilon, n_EOS, 3.D30, 3.D30, D2_epsilon)
  !  Call spline(EOS_rho_B, EOS_preassure, n_EOS, 3.D30, 3.D30, D2_pressure)

  end   subroutine EOS_initialize_solve

  end module EOS_initialize_module
