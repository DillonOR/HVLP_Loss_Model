% This is the O'Reilly approach
% [N,Lm_cc,Ipk_cc, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = OReilly(12,150000,1.7,0.9,8,0.9,6E-6,2000,1,'N87','EFD30' )

function [N,Lm_cc,Ipk_cc, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = OReilly(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge, CoreM, CoreS)
%% The O'Reilly Approach %% 

C           = Constants();
N           = ((V_out+Vd)*(1-D_max))/(V_in*D_max);  
[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM) 
B_max       = B_peak*0.65;

Ipk_cc = ((2*N*V_in+V_out)*C_out*V_out)/(n*V_in*(t_charge));     % Peak Current in LT3751 Datasheet
Lm_cc = 1/(f*Ipk_cc*((1/V_in)+(N/V_out)));                 % This is the O'Reilly Design Method for BCM inductor value

E_cc = 0.5*Lm_cc*(Ipk_cc^2);

T_on_cc  = (Lm_cc * Ipk_cc) / V_in; T_off_cc = (N * Lm_cc * Ipk_cc) / V_out;   Ts_cc  = T_on_cc + T_off_cc;                                            

f_cc   = 1 / Ts_cc;                                                   
D_on_cc  = T_on_cc  / Ts_cc; D_off_cc = T_off_cc / Ts_cc;  Dsum_cc  = D_on_cc + D_off_cc;                                            % [-]  (should be 1 if Ts_cc = Ton+Toff)
f_check_cc = 1 / ( Lm_cc * Ipk_cc * ( (1/V_in) + (N/V_out) ) );     

[Core,Ac, lc, Wa, Ap, Vc, MLT, Al, G] = Geometry(CoreS, 'Core_Size')
Pcu = 0.4; % In Ericksons Method for calculating core size - he talks about 'allowed' loss#
Ku = 0.3;
Kg = ((C.rho_w*(Lm_cc^2)*(Ipk_cc^4))/(3*(B_max^2)*Pcu*Ku))*1E8;

Lg_cc       = ((C.mu_0*Lm_cc*(Ipk_cc^2))/((B_max^2)*Ac))*1E5;
Fr          = (1+(Lg_cc/sqrt((Ac*10000)))*log((2*G/Lg_cc)));                                                  % This is the fringing flux factor 
Np          = sqrt((Lm_cc*(Lg_cc+(lc/perm)))/((Ac*10000)*Fr*(C.mu_0)));
%Np_cc       = ((Lm_cc*Ipk_cc)/(Ac*B_max));
Ns_cc       = Np*N;


% Calculation of the secondary winding parasitic capacitance - I'm only
% using the Bank winding capacitance estimate method. 
% CITATION:
% https://ieeexplore.ieee.org/abstract/document/7299303
Ns_l = 1;           % Number of secondary layers
qs = 1;             % This is the number of sections used. 
ins_sin = 1.10;     % single-build enamel (NEMA Type 1)
ins_dou = 1.15;     % heavy-build enamel (NEMA Type 2)
ins_tri = 1.70;     % triple-insulated wire (TIW) typical

bare_wireD = sqrt(4*A_ws/pi)*1E-3;  
out_wireD = (bare_wireD*ins_tri);
T_pl  = ceil(Ns/Ns_l);                   % turns per layer
pitch = tape_t+out_wireD;                   % layer to layer pitch

layer_w =(T_pl - 1)*pitch + out_wireD;      % [m] axial layer width

eff_di = pitch - 1.15*bare_wireD + 0.26*out_wireD;

% Use bank winding - we know it works from paper 
L_a = (Ns_l-1)*out_wireD+out_wireD;
Cs =  mu_0fm*kapton*(out_wireD*(MLT*1E-2)/eff_di)*(L_a/layer_w);
Pcs = 0.5*(Cs/qs)*(V_out^2)*f; % This is the loss associated to the parasitic capacitance 


%% Integrated General Steinmetz Equations (For Non-Sinusoidal Excitations)
%       CITATION:
%       https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1196712
        Vref = V_out/N;
        S1 = abs(((V_in/(Np*Ac))^alpha)*D1); % I need to confirm values for Ac and that Ae = Ac - Ac here is actually Ae
        S2 = abs(((Vref/(Np*Ac))^alpha)*D2); 
        ki = Kc / ( 2^(beta+1)*pi^(alpha-1)*(0.2761 + (1.7061/(alpha + 1.354))));
        P_loss = (ki * (delta_B^(beta - alpha)) * (S1 + S2))*Vc;
end