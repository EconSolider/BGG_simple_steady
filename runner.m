cd C:\Users\lzkzz\Desktop\matlearn\Dynare_train\BGG_ver2
clc,clear;
%% parameters
sigm=0.48; 
mu=0.12; 
bet=0.99; 
zeta=3.33; 
alph=0.35;
chi=1-0.64/(1-alph); 
delt=0.025; 
thet=0.8;
rhoa=0.9;
phipi=1.5; 
rho=0.9;
epsilon=11;
Omega=0.8;
xi=1.4;    %设定方法:使得H尽量接近1/3

%% steady state
% 调整gamma_e，校准L的稳态，并求解omegabar的稳态
L_ss=3/2;
x0=0.4;
para=[sigm mu L_ss];
omegabar_ss=fsolve(@(x)comp_omegabar(x,para),x0);
F_ss=normcdf((log(omegabar_ss))/sigm+sigm/2);
Fp_ss=1/(sigm*omegabar_ss)*normpdf((log(omegabar_ss))/sigm+sigm/2);
G_ss=normcdf((log(omegabar_ss))/sigm-sigm/2);
Gp_ss=1/sigm*normpdf(log(omegabar_ss)/sigm+sigm/2);
Gamma_ss=G_ss+omegabar_ss*(1-F_ss);
Gammap_ss=1-F_ss;
S_ss=(L_ss-1)/L_ss/(Gamma_ss-mu*G_ss);

A_ss=1;
pi_ss=1;
pistar_ss=1;
pm_ss=(epsilon-1)/epsilon;
D_ss=1;
R_ss=1/bet;
in_ss=R_ss;
Rk_ss=S_ss*R_ss;
I_K=delt;
Q_ss=1;
Ye_K=1/(alph*pm_ss)*(Rk_ss-1+delt);
Y_K=Ye_K;
V_K=(1-Gamma_ss)*S_ss*R_ss*Q_ss;
He_ss=1;
we_K=(1-alph)*(1-Omega)*pm_ss*Y_K;
%注意此处调整gamma_e
gamma_e=(1/L_ss-we_K)/V_K;

N_K=1/L_ss;
Ce_K=(1-gamma_e)*V_K;
C_K=Y_K-Ce_K-I_K-mu*G_ss*Rk_ss;
H_ss=(xi*C_K*((1-alph)*Omega*pm_ss*Ye_K)^(-1)+1)^-1;
K_ss=H_ss^Omega * (Ye_K)^(-1/(1-alph));

I_ss=I_K*K_ss;
Ye_ss=Ye_K*K_ss;
Y_ss=Ye_ss;
V_ss=V_K*K_ss;
we_ss=we_K*K_ss;
N_ss=N_K*K_ss;
Ce_ss=Ce_K*K_ss;
C_ss=C_K*K_ss;
w_ss=xi/(1-H_ss)*C_ss;
lambd_ss=1/C_ss;
x1_ss=lambd_ss*pm_ss*Y_ss/(1-thet*bet);
x2_ss=lambd_ss*Y_ss/(1-thet*bet);
M_ss=zeta/(lambd_ss-bet*lambd_ss);
save _param;

%% IRF
dynare BGG_ver2;



