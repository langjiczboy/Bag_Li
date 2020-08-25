module common_module
  implicit none
  integer::n_EOS
  character(len=200)file_read,file_write

  real*8,allocatable::EOS_rho_B(:),EOS_epsilon(:),EOS_pressure(:),EOS_D_pre_rho(:)
  real*8,allocatable::EOS_energy(:),EOS_D_pre_epsilon(:)
  real*8,allocatable::D2_epsilon_pressure(:),D2_rho_B_pressure(:),D2_D_pre_epsilon(:)  
  integer::n_var=3
  integer,Parameter::nmax=50, nstpmx=50000

  
  end module common_module
