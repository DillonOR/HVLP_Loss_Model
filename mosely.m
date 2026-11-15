function [AgLg, N_pm, N_sm, Lm_e, B_opt ] = mosely(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge, Lg, CoreM, CoreS)

%% The Energy Storage Approach %%
%   CITATION:
%   https://drive.google.com/file/d/1fAjinR0-encuzV9i8NZz0HiL_FuuEbAW/view
%   [AgLg, N_pm, N_sm, Lm_e, B_opt ] = mosely(12,150000,1.7,0.9,8,0.9,6E-6,2000,0, 0.001,'N87','EFD30' )
% With this approach - you either need to estimate Lm or Ipk to start the
% design process. 
C = Constants();
[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM) 
B_max       = B_peak*0.65;

N               = ((V_out+Vd)*(1-D_max))/(V_in*D_max);                                                              % Call Core Size Parameters 
Ipk_e           =   3;                                      % An estimated peak current based on P_out requirement 
E_C             =   0.5*C_out*(V_out^2);
E_pulse         =   E_C/f;                                 % Estimated number of energy pulses based on a BCM frequency and 6 joule thruster
Lm_e            =   (2*E_pulse)/(Ipk_e^2);                  % Primary side inductance calculation 
E_num_pulses    =   E_pulse*f;                              % The Number of pulses required to reach 6 joules of energy 
AgLg            =   C.mu_0*Lm_e*((Ipk_e/B_max)^2)*1E6;      % Defines required gap volume Calculation of the Gap volume multiply by 10^9 to get in mm3 or 10^6 for cm3
E_AgLg          =   ((1/(2*C.mu_0))*AgLg*(B_max^2))*1E-6;
B_opt           =   (sqrt(E_AgLg/((1/(2*C.mu_0))*(AgLg*1E-6))));
%Lg_range        =   0.001:0.01:0.05;
[Core, Ac, lc, Wa, Ap, Vc, MLT, Al, G] = Geometry(CoreS, 'Core_Size')  
%Convert Ac to mmsquared
Ac              = Ac*100;

N_pm            = sqrt((Lm_e*Lg)/(C.mu_0*Ac));
N_sm            = N_pm*N;
%[candidates, iCore, iGap] = sweep_core_gap(AgLg, E_pulse, B_max, Lg_range, Lm_e, N)
%[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM)      % Call Core Material Parameters
%  








end 