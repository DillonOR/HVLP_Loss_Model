% [N,Lm_er,Ipk_e, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = erickson(12,150000,1.7,0.9,8,0.9,6E-6,2000,1.5,'N87','EFD30')
function [N,Lm_er,Ipk_e, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = erickson(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge,CoreM, CoreS)


C           = Constants();
R_eq         = (V_in^2)/P_out;

[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM) 
B_max       = B_peak*0.65;

N           = ((V_out+Vd)*(1-D_max))/(V_in*D_max);   
Lm_e        = (((V_in^2)*(D_max^2))/(P_out))*((1/f)/2); 
Lm_er       = ((R_eq*((1-D_max)^2))/(2*(N^2)*f));
Ipk_e       = (V_in*D_max)/(Lm_e*f);


[Core,Ac, lc, Wa, Ap, Vc, MLT, Al, G] = Geometry(CoreS, 'Core_Size')

Lg              = (C.mu_0 * Lm_e * (Ipk_e^2) / (B_max^2 * Ac)) * 1e4;
N_pm            = sqrt((Lm_e*Lg)/(C.mu_0*(Ac*1E4)));

end