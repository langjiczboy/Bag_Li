#include "fintrf.h" !�����е�ͷ�ļ�,������mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix�Ⱥ���������

  SubRoutine mexFunction(out_num,p_out,in_num,p_in)!�����ӿ����Ʊ���ΪmexFunction,

  !out_num:�����������

  !p_out:�����������ָ��

  !in_num:�����������

  !p_in:�����������ָ��

  !����˳�����������



  Integer in_num,out_num
  mwPointer p_in(*),p_out(*)                           !mwPointerר�����ڱ�ʾָ�����,�������������Integer����

  mwPointer mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix, M_ptr ,R_ptr,rho_ptr,P_ptr,mxGetDoubles !����Է���ָ�뺯�����ٴ�����,

  integer::status,mxGetstring,mxCreateFull
  character*200::file_read,file_wirte
   Real*8::M_max,R_max,rho_max,P_max
   integer::n_EOS
  !****************����******************

  !*************Char*/**************
  status=mxGetstring(p_in(1),file_read,200)
  status=mxGetstring(p_in(2),file_wirte,200)
    if(status/=0) then
    call mexErrMsgTxt('�ַ�����ֵ����')
  end if
   !******integer
  Call main(file_read,file_wirte,M_max,R_max,rho_max,P_max)!�����ڲ�����,max_M,max_R

  !!*******************���******************
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
