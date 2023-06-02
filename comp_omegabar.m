function [err,F,G,Gamma,F_p,G_p,Gamma_p]=comp_omegabar(x,param)
sigma=param(1);
mu=param(2);
L=param(3);
F=normcdf((log(x))/sigma+sigma/2);
F_p=1/(sigma*x)*normpdf((log(x))/sigma+sigma/2);
G=normcdf((log(x))/sigma-sigma/2);
G_p=1/sigma*normpdf(log(x)/sigma+sigma/2);
Gamma=G+x*(1-F);
Gamma_p=1-F;

err=(Gamma_p-mu*G_p)/Gamma_p*(1-Gamma)-(Gamma-mu*G)/(L-1);
end