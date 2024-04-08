var
C ${C}$ (long_name='Consumption of workers')
lambd ${\lambda}$ (long_name='Largelange multiplier')
M ${M}$ (long_name='Money demand')
pi ${\pi}$ (long_name='Inflation rate')
R ${R}$ (long_name='Risk-free interest rate')
H ${H}$ (long_name='Labor of workers')
w ${w}$ (long_name='Real wage of workers')
F ${F}$ (long_name='Cdf of risk distribution')
G ${G}$ (long_name='Expected risk under cut-off value')
Gamma ${\Gamma}$ (long_name='Function of Ft and Gt')
omegabar ${\bar{\omega}}$ (long_name='Cut-off value')
S ${S}$ (long_name='Risk premium')
Rk ${Rk}$ (long_name='Real interest rate under risk')
L ${L}$ (long_name='Leverage of enterprenuer')
Q ${Q}$ (long_name='Price of capital goods')
K ${K}$ (long_name='Capital stock') 
N ${N}$ (long_name='Net wealth of enterprenuer')
Gammap ${{\Gamma}_p}$ (long_name='Gamma prime')
Gp ${G_p}$ (long_name='G prime')
V ${V}$ (long_name='Value of firms')
Ce ${C^e}$ (long_name='Consumption of enterprenuer')
He ${H^e}$ (long_name='Labor of enterprenuer')
we ${w^e}$ (long_name='Real wage of enterprenuer')
Ye ${Y^e}$ (long_name='output of enterprenuer firm')
A ${A}$ (long_name='TFP')
I ${I}$ (long_name='Investment of capital producer')
pistar ${\pi^*}$ (long_name='ratio of optimal and current price')
x1 x2 
D ${D}$ (long_name='Price dispersion')
Y ${Y}$ (long_name='Final goods')
in ${i}$ (long_name='Nominal interest rate controlled by CB')
pm ${p_m}$ (long_name='Real price of output of enterprenuer firm')
;

varexo
epsr ${\epsilon^r}$ (long_name='Monetary policy shock')
epsa ${\epsilon^a}$ (long_name='TFP shock')
;

parameters
sigm ${\sigma}$ 
mu ${\mu}$
bet ${\beta}$
zeta ${\zeta}$
alph ${\alpha}$
chi ${\chi}$
delt ${\delta}$
thet ${\theta}$
rhoa ${\rho^a}$
phipi ${\phi_{\pi}}$
rho ${\rho}$
epsilon ${\epsilon}$ 
Omega ${\Omega}$
xi ${\xi}$
gamma_e ${\gamma_e}$
%steady_state_values
C_ss lambd_ss M_ss pi_ss R_ss H_ss w_ss F_ss G_ss Gamma_ss 
omegabar_ss S_ss Rk_ss L_ss Q_ss K_ss N_ss Gammap_ss Gp_ss V_ss
Ce_ss He_ss we_ss Ye_ss A_ss I_ss pistar_ss x1_ss x2_ss D_ss
Y_ss in_ss pm_ss
;

load '_param.mat';
set_param_value('sigm',sigm);
set_param_value('mu',mu);
set_param_value('bet',bet);
set_param_value('zeta',zeta);
set_param_value('alph',alph);
set_param_value('chi',chi);
set_param_value('delt',delt);
set_param_value('thet',thet);
set_param_value('rhoa',rhoa);
set_param_value('phipi',phipi);
set_param_value('rho',rho);
set_param_value('epsilon',epsilon);
set_param_value('Omega',Omega);
set_param_value('xi',xi);
set_param_value('gamma_e',gamma_e);
set_param_value('C_ss',C_ss);
set_param_value('lambd_ss',lambd_ss);
set_param_value('M_ss',M_ss);
set_param_value('pi_ss',pi_ss);
set_param_value('R_ss',R_ss);
set_param_value('H_ss',H_ss);
set_param_value('w_ss',w_ss);
set_param_value('F_ss',F_ss);
set_param_value('G_ss',G_ss);
set_param_value('Gamma_ss',Gamma_ss);
set_param_value('omegabar_ss',omegabar_ss);
set_param_value('S_ss',S_ss);
set_param_value('Rk_ss',Rk_ss);
set_param_value('L_ss',L_ss);
set_param_value('Q_ss',Q_ss);
set_param_value('K_ss',K_ss);
set_param_value('N_ss',N_ss);
set_param_value('Gammap_ss',Gammap_ss);
set_param_value('Gp_ss',Gp_ss);
set_param_value('V_ss',V_ss);
set_param_value('Ce_ss',Ce_ss);
set_param_value('He_ss',He_ss);
set_param_value('we_ss',we_ss);
set_param_value('Ye_ss',Ye_ss);
set_param_value('A_ss',A_ss);
set_param_value('I_ss',I_ss);
set_param_value('pistar_ss',pistar_ss);
set_param_value('x1_ss',x1_ss);
set_param_value('x2_ss',x2_ss);
set_param_value('D_ss',D_ss);
set_param_value('Y_ss',Y_ss);
set_param_value('in_ss',in_ss);
set_param_value('pm_ss',pm_ss);

model;
1/C-lambd=0;
zeta/M-lambd+bet*lambd(+1)/pi(+1)=0;
-lambd+bet*lambd(+1)*R=0;
-xi/(1-H)+lambd*w=0;
F=normcdf(log(omegabar)/sigm+sigm/2);
G=normcdf(log(omegabar)/sigm-sigm/2);
Gamma=G+omegabar*(1-F);
 %由于模型中的变量不能只有(+1)时点的，因此S的定义需要换一下
S=Rk(+1)/R;
L=Q*K/N;                        %资本存量是状态变量，净资产是状态变量
(Gammap-mu*Gp)/Gammap*(1-Gamma)+(Gamma-mu*G)=1/S;
(Gamma-mu*G)*S=1-1/L;
V=(1-Gamma)*S*R*Q*K;           %资本存量是状态变量
Ce=(1-gamma_e)*V;
N=gamma_e*V+He*we; %企业家净资产是状态变量
Ye=A*K(-1)^alph*(H^Omega*He^(1-Omega))^(1-alph); %资本存量是状态变量
w*H=(1-alph)*Omega*pm*Ye;
we*He=(1-alph)*(1-Omega)*pm*Ye;
Rk=alph*pm*Ye/Q(-1)/K(-1)+Q*(1-delt)/Q(-1); %资本存量是状态变量
pistar=epsilon/(epsilon-1)*x1/x2;
x1=lambd*pm*Y+thet*bet*x1(+1)*(pi(+1))^epsilon;
x2=lambd*Y+thet*bet*x2(+1)*(pi(+1))^(epsilon-1);
1=thet*pi^(epsilon-1)+(1-thet)*pistar^(1-epsilon);
K=(1-delt)*K(-1)+I-chi/2*(I/K(-1)-delt)^2*K(-1); %资本存量是状态变量
Q=1/(1-chi*(I/K(-1)-delt)); %资本存量是状态变量
Ye=D*Y;
D=thet*pi^epsilon*D(-1)+(1-thet)*pistar^(-epsilon);
Y=C+Ce+I+chi/2*(I/K(-1)-delt)^2*K(-1)+mu*G*Rk*Q(-1)*K(-1); %资本存量是状态变量
R=in/pi(+1);
log(in/steady_state(in))=rho*log(in/steady_state(in))+phipi*log(pi/steady_state(pi))+epsr;
log(A)=rhoa*log(A(-1))+epsa;
He=1;
Gammap=1-F;
Gp=normpdf(log(omegabar)/sigm+sigm/2)/sigm;
end;

steady_state_model;
C=C_ss;
lambd=lambd_ss;
M=M_ss;
pi=pi_ss;
R=R_ss;
H=H_ss;
w=w_ss;
F=F_ss;
G=G_ss;
Gamma=Gamma_ss;
omegabar=omegabar_ss;
S=S_ss;
Rk=Rk_ss;
L=L_ss;
Q=Q_ss;
K=K_ss;
N=N_ss;
Gammap=Gammap_ss;
Gp=Gp_ss;
V=V_ss;
Ce=Ce_ss;
He=He_ss;
we=we_ss;
Ye=Ye_ss;
A=A_ss;
I=I_ss;
pistar=pistar_ss;
x1=x1_ss;
x2=x2_ss;
D=D_ss;
Y=Y_ss;
in=in_ss;
pm=pm_ss;
end;
steady;
check;

shocks;
var epsa=0.04^2;
%var epsr=0.01^2;
end;
stoch_simul(order=1,irf=80,periods=0);










