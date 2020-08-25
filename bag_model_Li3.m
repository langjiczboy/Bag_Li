%%%º∆À„ Zhou En-Ping 2018 PRD
%%% P.Haensel,J.L.Zdunik 1986 Strange quark stars


clear
close all
clc
%% path
TOV_main_path='D:\sj\code\Fortran\Called\TOV_MAT';
bag_main_path='D:\sj\code\MATLAB\bag';
%% const
load('const_num.mat')

%% ini
TOV_head='QS-MR-';
col_TOV=6;
M_care=1.4;
%% par
% a_4=0.754;
% B_eff=131^4;
% m_s=100;
a_4=1;
B_eff=0;
m_s=0;
N_c=3;
%% calculate
%%% vec mu 
syms n_b
n_b_vec=flip([0:0.01:10]');
n_b_vec=3./exp(n_b_vec);
P_all=NaN;
epsilon_all=NaN;

epsilon_sym=(9*pi^(2/3)/4)*n_b^(4/3)/a_4^(1/3)+B_eff;
for num_b=1:length(n_b_vec)
    n_b_ans=n_b_vec(num_b)*fmm32MeV3;
    epsilon=subs(epsilon_sym,n_b,n_b_ans);
    epsilon_all(num_b,1)=epsilon/fmm32MeV3;

    P=subs((epsilon-4*B_eff)/3,n_b,n_b_ans);
    P_all(num_b,1)=P/fmm32MeV3;

end
E_all=epsilon_all./n_b_vec;
[min_E,min_E_posi]=min(E_all);
min_n=n_b_vec(min_E_posi);
plot(n_b_vec,E_all,'.')
%% com
cd(bag_main_path)
bag_filename=['Bag Li ','a_4 ',num2str(a_4),' B ',num2str(B_eff^(1/4)),' m_s ', num2str(m_s),'.txt'];
diff_P_epcilon=diff_vec(epsilon_all,P_all)';
diff_P_rho=diff_vec(n_b_vec,P_all)';

com_bag_name={'n_B','E','Pressure','epsilon','diff_P_rho','diff_P_epcilon'};
fid=fopen(['com_',bag_filename],'w');
row_bag=size(com_bag_name,2);
for i=1:row_bag
    fprintf(fid,'%15s',char(com_bag_name{i}));
end
fprintf(fid,'\n');
bag_data=[n_b_vec,E_all,P_all,epsilon_all,diff_P_rho,diff_P_epcilon];
fprintf(fid,[repmat('%15.8f',[1,row_bag]),'\n'],bag_data');
fclose(fid);
%% TOV
cd(TOV_main_path)

[M_max_bag,R_max_bag,rho_max_bag,P_max_bag]=TOV_MAT([bag_main_path,'\com_',bag_filename],[TOV_head,bag_filename]);

%% Lambda_tilde
fid=fopen([TOV_head,bag_filename]);
TOV_txt=textscan(fid,repmat('%f',[1,col_TOV]),'headerlines',1);
fclose(fid);
Radius=TOV_txt{:,1};
Mass=TOV_txt{:,2};
Pressure=TOV_txt{:,3};
rho=TOV_txt{:,4};
Lambda=TOV_txt{:,5};
loven=TOV_txt{:,6};
[~,position_1d4]=min(abs(Mass-M_care));
Lambda_1d4=Lambda(position_1d4);
[lambda_tild17_low,lambda_tild17_high]=GW17(Mass,Lambda);
if isempty(lambda_tild17_low)
    lambda_tild17_low=NaN;
end
%              lambda_tild17_low=NaN;
lambda_tild17_low_min=min(lambda_tild17_low);
lambda_tild17_low_max=max(lambda_tild17_low);
