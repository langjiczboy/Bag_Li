%%%计算 Zhou En-Ping 2018 PRD
%%% P.Haensel,J.L.Zdunik 1986 Strange quark stars
%%%以mu为自变量

clear
close all
clc
tic
%% path
TOV_main_path='D:\sj\code\Fortran\Called\TOV_MAT';
bag_main_path='D:\sj\code\MATLAB\bag';
%% 编译TOV
cd(TOV_main_path)
delete *.mexw64
mex -largeArrayDims TOV_MAT.f90 common_module.f90 const_module.f90 spline_module.f90 EOS_initialize_module.f90 rk45_module.f90  star_solve_module.f90  main.f90

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
m_s=100;
N_c=3;

%% calculate

%%% vec
a_4_vec=[0.5:0.1:1]';
B_vec=[130:5:145]';
mu_vec=(200:2:700)';
bag_Li_mat=NaN;
num_mat=1;
for num_B=1:length(B_vec)
    B_eff=B_vec(num_B)^4;
    disp('--------------------')
    disp(['B = ',num2str(B_vec(num_B))])
    for num_a_4=1:length(a_4_vec)
        a_4=a_4_vec(num_a_4);
        
        disp(['a_4 = ',num2str(a_4)])
        P_all_tran=NaN;
        epsilon_all_tran=NaN;
        n_all_tran=NaN;
        %% Omega
        nu_s=sqrt(mu_s^2-m_s^2);
        Omega_u=-N_c*mu_u^4/(12*pi^2);
        Omega_d=-N_c*mu_d^4/(12*pi^2);
        Omega_e=-mu_e^4/(12*pi^2);
        
        Omega_s0=mu_s*nu_s*(mu_s^2-5/2*m_s^2)+3/2*m_s^4*acosh(mu_s/m_s);
        Omega_s=(Omega_s0)/(-4*pi^2);
        Omega=(Omega_u+Omega_d+Omega_s+Omega_e)+3/(4*pi^2)*(1-a_4)*(mu_b/3)^4+B_eff;
        %% n
        n_u=N_c/(3*pi^2)*mu_u^3-1/pi^2*(1-a_4)*(mu_b/3)^3;
        n_d=N_c/(3*pi^2)*mu_d^3-1/pi^2*(1-a_4)*(mu_b/3)^3;
        n_s=-diff(Omega_s,mu_s)-1/pi^2*(1-a_4)*(mu_b/3)^3;
        n_e=-diff(Omega_e,mu_e);
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
            [mu_e_ans]=vpasolve([2*n_u==n_d+n_s+3*n_e],mu_e,[0 1e3]);
            
            P=subs(-Omega,[mu,mu_e],[mu_ans,mu_e_ans]);
            epsilon=subs(Omega+mu_u*n_u+mu_d*n_d+mu_s*n_s+mu_e*n_e,[mu,mu_e],[mu_ans,mu_e_ans]);
            P_all_tran(num_mu,1)=P/fmm32MeV3;
            epsilon_all_tran(num_mu,1)=epsilon/fmm32MeV3;
            n_all_tran(num_mu,1)=subs(n_u,mu_e,mu_e_ans)/fmm32MeV3;
            n_all_tran(num_mu,2)=subs(n_d,mu_e,mu_e_ans)/fmm32MeV3;
            n_all_tran(num_mu,3)=subs(n_s,mu_e,mu_e_ans)/fmm32MeV3;
            n_all_tran(num_mu,4)=subs(n_e,mu_e,mu_e_ans)/fmm32MeV3;
        end
        
        n_B_vec_tran=sum(n_all_tran(:,1:3),2)/3;
        E_all_tran=epsilon_all_tran./n_B_vec_tran;
        [min_E,min_E_posi]=min(E_all_tran);
        
        min_n=n_B_vec_tran(min_E_posi);
        E_all=E_all_tran(min_E_posi-1:end);
        P_all=P_all_tran(min_E_posi-1:end);
        n_B_vec=n_B_vec_tran(min_E_posi-1:end);
        epsilon_all=epsilon_all_tran(min_E_posi-1:end);
        if isempty(find(P_all<0, 1))
            error('没有零压点')
        end
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
        fprintf(fid,[repmat('%15.8g',[1,row_bag]),'\n'],bag_data');
        fclose(fid);
        %% TOV
        cd(TOV_main_path)
        
        [M_max_bag,R_max_bag,rho_max_bag,P_max_bag]=TOV_MAT([bag_main_path,'\com_',bag_filename],[TOV_head,bag_filename]);
        if P_max_bag>max(P_all)
            warning('最大化学势过小')
        end
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
        R_1d4=Radius(position_1d4);
        [lambda_tild17_low,lambda_tild17_high]=GW17(Mass,Lambda);
        if isempty(lambda_tild17_low)
            lambda_tild17_low=NaN;
        end
        %              lambda_tild17_low=NaN;
        lambda_tild17_low_min=min(lambda_tild17_low);
        lambda_tild17_low_max=max(lambda_tild17_low);
        bag_Li_mat(num_mat,1)=m_s;
        bag_Li_mat(num_mat,2)=B_vec(num_B);
       bag_Li_mat(num_mat,3)=a_4_vec(num_a_4);

        bag_Li_mat(num_mat,4)=min(E_all);
        bag_Li_mat(num_mat,5)=R_max_bag;%M_max_CLF,R_max_CLF
        bag_Li_mat(num_mat,6)=M_max_bag;
        bag_Li_mat(num_mat,7)=Lambda_1d4;
        bag_Li_mat(num_mat,8)=lambda_tild17_low_min;
        bag_Li_mat(num_mat,9)=lambda_tild17_low_max;
        num_mat=num_mat+1;
    end
end
%% ans
bag_Li_mat_set=mat2dataset(bag_Li_mat,'VarNames',{'m_s','B','a_4','E_min','R_max','M_max','Lambda_1d4','lambda_tild17_low_min','lambda_tild17_low_max'});
scatter(bag_Li_mat_set.a_4,bag_Li_mat_set.B,log(bag_Li_mat_set.M_max)*15,(bag_Li_mat_set.lambda_tild17_low_min),'filled')
colormap(jet)
caxis([0 720])
hold on

h=colorbar;
ylabel('B^{1/4}_{eff}(MeV)')
xlabel('a_4')
set(get(h,'title'),'string','$\tilde{\Lambda}(q=0.73)$','Interpreter','latex')
saveas(gcf,[bag_main_path,'\Bag-Li2-B-a4-Lambda.fig'])
save([bag_main_path,'\Bag_Li2.mat'],'bag_Li_mat_set')
figure
scatter(bag_Li_mat_set.a_4,bag_Li_mat_set.B,[],(bag_Li_mat_set.M_max),'filled')
colormap(jet)
hold on

h=colorbar;
ylabel('B^{1/4}_{eff}(MeV)')
xlabel('a_4')
set(get(h,'title'),'string','$M_{TOV}$','Interpreter','latex')
saveas(gcf,[bag_main_path,'\Bag-Li2-B-a4-Mass.fig'])

% toc