#include "fintrf.h" !必须有的头文件,里面有mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix等函数的申明

  SubRoutine mexFunction(out_num,p_out,in_num,p_in)!函数接口名称必须为mexFunction,

  !out_num:输出参数个数

  !p_out:输出参数数组指针

  !in_num:输入参数个数

  !p_in:输入参数数组指针

  !参数顺序不能随意更改



  Integer in_num,out_num
  mwPointer p_in(*),p_out(*)                           !mwPointer专门用于表示指针变量,这个不能随意用Integer代替

  mwPointer mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix, M_ptr ,R_ptr,rho_ptr,P_ptr,mxGetDoubles !这个对返回指针函数的再次申明,

  integer::status,mxGetstring,mxCreateFull
  character*200::file_read,file_wirte
   Real*8::M_max,R_max,rho_max,P_max
   integer::n_EOS
  !****************输入******************

  !*************Char*/**************
  status=mxGetstring(p_in(1),file_read,200)
  status=mxGetstring(p_in(2),file_wirte,200)
    if(status/=0) then
    call mexErrMsgTxt('字符串赋值出错')
  end if
   !******integer
  Call main(file_read,file_wirte,M_max,R_max,rho_max,P_max)!调用内部函数,max_M,max_R

  !!*******************输出******************
  p_out(1) = mxCreateDoubleMatrix(1,1,0)
  p_out(2) = mxCreateDoubleMatrix(1,1,0)
  p_out(3) = mxCreateDoubleMatrix(1,1,0)
  p_out(4) = mxCreateDoubleMatrix(1,1,0)

  M_ptr = mxGetPr(p_out(1))
  R_ptr = mxGetPr(p_out(2))
  rho_ptr=mxGetPr(p_out(3))
  P_ptr=mxGetPr(p_out(4))

  call mxCopyReal8ToPtr(M_max,M_ptr,1)
  call mxCopyReal8ToPtr(R_max,R_ptr,1)
  call mxCopyReal8ToPtr(rho_max,rho_ptr,1)
  call mxCopyReal8ToPtr(P_max,P_ptr,1)
  Return

  End SubRoutine
