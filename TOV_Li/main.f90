  subroutine  main(file_read_in,file_write_in,M_max,R_max,rho_max,P_max)
  !通过已知的状态方程(EOS)，通过解TOV方程，给出M-r  P-r关系
  use const_module
  use EOS_initialize_module
  use star_solve_module
  use spline_module
  use common_module
  implicit none
  real*8::M_max,R_max,rho_max,P_max
  real*8::M_star,R_star,rho_B
  real*8::Lambda,Loven
  integer::ii,jj
  real*8::Pressure_c,Pressure_middle,Pressure_end
  integer::nRow
   Integer,External::GetFileN
       character(len=200)file_read_in,file_write_in
!*************公用区域**********************************
   file_read=file_read_in
   file_write=file_write_in
  Open( 12 , File = file_read )
  nRow = GetFileN( 12 )
  Close( 12 )
  
  n_EOS=nRow-1
  allocate(EOS_rho_B(n_EOS),EOS_epsilon(n_EOS),EOS_pressure(n_EOS),D2_epsilon_pressure(n_EOS),D2_rho_B_pressure(n_EOS))
allocate(EOS_D_pre_rho(n_EOS),EOS_energy(n_EOS),EOS_D_pre_epsilon(n_EOS),D2_D_pre_epsilon(n_EOS))
  !**********导入已知的数据，并给出插值用的二阶导 ********************
  call const_solve()
  call EOS_initialize_solve()
  M_max = 0.D0
  R_max = 0.D0

  !  No_star = 0
  ii = 1
  jj = 1
  !  rho_c_in =   0.56D0 !cpc larger then rho(p=0)  DI1700 0.5  DI255 0.4

  Pressure_c = 1D0
  open(50,file=file_write)
  write(50,'(10(A15))')'Radius','Mass','Pressure','rho','Lambda','loven'
  write(*,'(10(A15))')'Radius','Mass','Pressure','rho','Lambda','loven'
  Pressure_middle=100d0

  Pressure_end=1000d0
  Do

    If (Pressure_c<Pressure_middle) Then
      Pressure_c = ii*1D0
      ii = ii + 1
    Else If (Pressure_c>=Pressure_middle.and. Pressure_c<Pressure_end) Then
      Pressure_c = jj*1D0 +Pressure_middle
      jj = jj + 1
    Else
      Exit
    End If

    Call star_solve(Pressure_c, m_star, r_star,Lambda,loven)
    call splint(EOS_Pressure,EOS_rho_B, D2_rho_B_pressure, n_EOS, Pressure_c, rho_B)  !

    write(50,'(10(G15.8))')R_star,M_star,Pressure_c,rho_B,Lambda,loven
    write(*,'(10(G15.8))')R_star,M_star,Pressure_c,rho_B,Lambda,loven

    If (M_star>=M_max) Then
      M_max = M_star
      R_max = R_star
      rho_max=rho_B
      P_max=Pressure_c
    End If

  End Do

  write(50,'(2A15)') 'R_max','M_max'
  write(50,'(2(F15.8))')R_max,M_max

  write(*,'(2A15)')  'R_max','M_max'
  write(*,'(2(F15.8))')R_max,M_max
  deallocate(EOS_rho_B,EOS_epsilon,EOS_pressure,D2_epsilon_pressure,D2_rho_B_pressure)
  deallocate(EOS_D_pre_rho,EOS_energy,EOS_D_pre_epsilon,D2_D_pre_epsilon)
close(50)
  End subroutine main
  
  
        Integer Function GetFileN( iFileUnit )
  !// 此函数应在打开文件后立即调用。调用后读取位置返回文件起始位置
    Implicit None
    Integer , Intent( IN ) :: iFileUnit
    character( Len = 1 ) :: cDummy
    integer :: ierr
    GetFileN = 0
    Rewind( iFileUnit )
    Do
      Read( iFileUnit , * , ioStat = ierr ) cDummy
      If( ierr /= 0 ) Exit
      GetFileN = GetFileN + 1
    End Do
    Rewind( iFileUnit )
    End Function GetFileN