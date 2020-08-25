  module const_module
  implicit none

  !****************  Some constants  *************************************
  real*8::pi,cc,cg,cmevg,crcmg,crcg,crcs,cmsun
  contains 
  subroutine const_solve
  implicit none
  pi=asin(1.)*2
  cc =     2.99792458D0 !(10^(10) cm.s-1

  cg =6.67259D0  ! 6.67408d0            !(10^(-8) cm3.g-1.s-2

  cmevg =  1.782662712D0 !1MeV/c2=1.782662712 *10^(-27) g

  crcmg = cc*cc*1.D6/4.D0/pi/cmevg/cg
  crcg = 0.5*1.0D43*cc*cc/cg
  crcs = 1.0D-10/cc**2
  cmsun = 4.0D0*pi*cmevg/1.98855D6
  end   subroutine const_solve

  !***********************************************************************

  end module const_module
