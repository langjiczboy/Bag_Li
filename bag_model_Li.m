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
syms mu mu_e
mu_s=mu;
mu_d=mu;
mu_u=mu-mu_e;
mu_b=mu_u+mu_d+mu_s;
TOV_head='QS-MR-';
col_TOV=6;
M_care=1.4;
%% par
a_4=0.61;
B_eff=138^4;
m_s=0;
N_c=3;
%% Omega
nu_s=sqrt(mu_s^2-m_s^2);

Omega_u=-N_c*mu_u^4/(12*pi^2);
Omega_d=-N_c*mu_d^4/(12*pi^2);
Omega_e=-mu_e^4/(12*pi^2);

Omega_s0=mu_s*nu_s*(mu_s^2-5/2*m_s^2)+3/2*m_s^4*acosh(mu_s/m_s);
Omega_s=(Omega_s0)/(-4*pi^2);
% n_s=(mu_s^2-m_s^2)/pi^2*((mu_s^2-m_s^2)^(1/2));
Omega=(Omega_u+Omega_d+Omega_s+Omega_e)+3/(4*pi^2)*(1-a_4)*(mu_b/3)^4+B_eff;
%% n
%Free
% n_u=N_c*mu_u^3/(3*pi^2);
% n_d=N_c*mu_d^3/(3*pi^2);
% n_s=-diff(Omega_s,mu_s);
% Account
n_u=N_c/(3*pi^2)*mu_u^3-3/pi^2*(1-a_4)*(mu_b/3)^3;
n_d=N_c/(3*pi^2)*mu_d^3-3/pi^2*(1-a_4)*(mu_b/3)^3;
n_s=-diff(Omega_s,mu_s)-3/pi^2*(1-a_4)*(mu_b/3)^3;


n_e=-diff(Omega_e,mu_e);
% clear
% syms mu_u mu_d mu_s pi alpha_c pi m_s rho mu_s p_s ln_s mu_e nu_s
%% calculate
% %%% vec n_B 
% n_B_vec=(0.3:0.01:3)';
% P_all=NaN;
% epsilon_all=NaN;
% n=NaN;
% for num_B=1:length(n_B_vec)
%     n_B=n_B_vec(num_B)*fmm32MeV3;
%     [mu_ans,mu_e_ans]=vpasolve([n_B==(n_u+n_d+n_s)/3,2*n_u==n_d+n_s+3*n_e],[mu,mu_e],[-1e4 1e4; -1e4 1e4]);
%     P=subs(-Omega,[mu,mu_e],[mu_ans,mu_e_ans]);
%     epsilon=subs(Omega+mu_u*n_u+mu_d*n_d+mu_s*n_s+mu_e*n_e,[mu,mu_e],[mu_ans,mu_e_ans]);
%     P_all(num_B,1)=P/fmm32MeV3;
%     epsilon_all(num_B,1)=epsilon/fmm32MeV3;
%     n(num_B,1)=subs(n_u,[mu,mu_e],[mu_ans,mu_e_ans])/fmm32MeV3;
%     n(num_B,2)=subs(n_d,[mu,mu_e],[mu_ans,mu_e_ans])/fmm32MeV3;
%     n(num_B,3)=subs(n_s,[mu,mu_e],[mu_ans,mu_e_ans])/fmm32MeV3;
%     n(num_B,4)=subs(n_e,[mu,mu_e],[mu_ans,mu_e_ans])/fmm32MeV3;
%     
% end

%%% vec mu 
mu_vec=(100:1:630)';
P_all=NaN;
epsilon_all=NaN;
n=NaN;
n_u_sym=n_u;
n_d_sym=n_d;
n_s_sym=n_s;
n_e_sym=n_e;
for num_mu=1:length(mu_vec)
    mu_ans=mu_vec(num_mu);
    n_u=subs(n_u_sym,mu,mu_ans);
        n_d=subs(n_d_sym,mu,mu_ans);
    n_s=subs(n_s_sym,mu,mu_ans);
    n_e=subs(n_e_sym,mu,mu_ans);

    [mu_e_ans]=vpasolve([2*n_u==n_d+n_s+3*n_e],mu_e,[-1e4 1e4]);
    P=subs(-Omega,[mu,mu_e],[mu_ans,mu_e_ans]);
    epsilon=subs(Omega+mu_u*n_u+mu_d*n_d+mu_s*n_s+mu_e*n_e,[mu,mu_e],[mu_ans,mu_e_ans]);
    P_all(num_mu,1)=P/fmm32MeV3;
    epsilon_all(num_mu,1)=epsilon/fmm32MeV3;
    n(num_mu,1)=subs(n_u,mu_e,mu_e_ans)/fmm32MeV3;
    n(num_mu,2)=subs(n_d,mu_e,mu_e_ans)/fmm32MeV3;
    n(num_mu,3)=subs(n_s,mu_e,mu_e_ans)/fmm32MeV3;
    n(num_mu,4)=subs(n_e,mu_e,mu_e_ans)/fmm32MeV3;
    
end
n_B_vec=sum(n(:,1:3),2)/3;
E_all=epsilon_all./n_B_vec;
min(E_all)
plot(n_B_vec,E_all,'.')
%% com
cd(bag_main_path)
bag_filename=['Bag Li ','a_4 ',num2str(a_4),' B ',num2str(B_eff^(1/4)),' m_s ', num2str(m_s),'.txt'];
diff_P_epcilon=diff_vec(epsilon_all,P_all)';
diff_P_rho=diff_vec(n_B_vec,P_all)';

com_bag_name={'n_B','E','Pressure','epsilon','diff_P_rho','diff_P_epcilon'};
fid=fopen(['com_',bag_filename],'w');
row_bag=size(com_bag_name,2);
for i=1:row_bag
    fprintf(fid,'%15s',char(com_bag_name{i}));
end
fprintf(fid,'\n');
bag_data=[n_B_vec,E_all,P_all,epsilon_all,diff_P_rho,diff_P_epcilon];
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
