% This Code strictly computes transformer paramters using only Hurley's
% design principles, 
% Use this to try: [N,L_mh,Ipk_h, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = HurleyMethod(12,150000,1.7,0.9,8,0.9,6E-6,2000,1.5,'N87','EFD30')
    
function [N,L_mh,Ipk_h, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = HurleyMethod(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge, CoreM, CoreS)

C           = Constants();


[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM) 
B_max       = B_peak*0.65;

N           = ((V_out+Vd)*(1-D_max))/(V_in*D_max);   
L_mh        = (((V_in^2)*(D_max^2))/(P_out))*((1/f)/2);   %For DCM must be less than and for CCM must be greater than
Ipk_h       = (V_in*D_max)/(L_mh*f);


T_on        = (L_mh*Ipk_h)/V_in;        T_off = (N*L_mh*Ipk_h)/(V_out+Vd);                              % We similarly find the off time. 
Ts          = T_on + T_off;             f     = 1/Ts;                                                   % The above gives us the total switching period, we can now estimate frequency
D1          = T_on/Ts;                  D2    = T_off/Ts;                   Dc  = D1+D2;                % This is the conduction fraction of the cycle. In BCM it equals 1, in DCM less than 1.

Ilm_rms     = Ipk_h*sqrt(Dc/3); Ilm_avg = 0.5*Ipk_h*Dc;                                                 % In Hurley's method, he calls Ipk as Îp but this is for CCM, so in my case
                                                                                                        % Ipk maps to Îp and since change in current is simply Ipk, then yp = 1
ki          = 2/(sqrt(3*Dc));   Kip     = sqrt(D1/3);                                                   % This is the equivalent current waveform factor derived by Hurley but for DCM/BCM mode
Ip_rms      = Ipk_h*Kip;                                                                                % Since Îp = Ipk, then this is the equivalent Ip_rms primary in DCM/BCM

Is          = P_out/(D2*V_out); Is_avg = P_out/V_out;                                                   % The ripple of the secondary current is then a fraction of the turns
Ispk        = (2*P_out)/(D2*V_out);                                                                     % ratio referred to the primary side, in my case this is 1/N, so: 
Kis         = sqrt(D2/3);
Is_rms      = Ispk*Kis;
K_up        = C.Ku*(1/(1+(Is_rms/(N*Ip_rms))));                                                         %  need to identify what this value should be. a is equal to Np/Ns
E           = 0.5*L_mh*(Ipk_h^2);


Ap_calc     = ((((sqrt(1+C.y))*Kip*E)/(B_max*C.Kt*(K_up/sqrt(C.Ku))*sqrt(C.delta_T)))^(8/7))*1E8;       % Area Product Calculation Hurley 
[Core,Ac, lc, Wa, Ap, Vc, MLT, Al, G] = Geometry(CoreS, 'Core_Size')

R_theta     = 0.06/sqrt(Vc);                                                                            % R_theta = 0.06/sqrt(Vc); Vc should be in m^3 and R_theta is given in Celcius/Watt
P_md        = C.delta_T/R_theta;                                                                        % P_md is the max dissipoation 
Pcup        = (K_up/C.Ku)*(P_md);                                                                       % Not sure how y works here 
mu_opt      = (B_max*(lc*1E-2)*Kip)/( 4*pi*1e-7*sqrt((Pcup*K_up*(Wa*1E-4))/(C.rho_w*(MLT*1E-2))));          % Optimal Core Permeability 
Lg          = 0.038;%(((lc*10^(-2))/mu_opt)*1E3)/10;              % cm    - use 0.071 to check                                         % Use 0.071 to match with datasheet for check this is in cm
Ag          = Lg*Ac;                                      %                                              % I need to verify this formula %  area of the air gap
Fr          = (1+(Lg/sqrt((Ac*10000)))*log((2*G/Lg)));                                                  % This is the fringing flux factor 
Np          = sqrt((L_mh*(Lg+(lc/perm)))/((Ac*10000)*Fr*(C.mu_0)));
%Np          = sqrt(L_mh/(160E-9));  %This is from the datasheet, so it has the actual gap lenght specified and the Al value at 0.071 is 160E-9 
Np_b        = (L_mh*Ipk_h)/(Ac*B_max);                                                                  % I need to verify if Ae and Ac are the same thing
Ns          = Np*N;
Jo          = (C.Kt*(sqrt(C.delta_T))/((sqrt(C.Ku*(1+C.y)))*((Ap*1E-8)^(1/8))))*1E-4;                   % This comes out in A/cm^2 - Look at Hurley pg 82 for re-assurance - I am pretty confident in this now 
A_wp        = (Ip_rms/Jo)*100;
A_ws        = (Is_rms/Jo)*100;

[wire_rp, wire_rs] = getwire_resistance(A_wp, A_ws)

Rdc         = ((Np*MLT*(wire_rp*1E-6))*(1+C.alpha_20*(C.T_max-20)))*1E3;
P_cu        = (Rdc*1E-3)*(Ip_rms^2);
Rdc_s       = ((Ns*MLT*(wire_rs*1E-6))*(1+C.alpha_20*(C.T_max-20)))*1E3;
P_cu_s      = (Rdc_s*1E-3)*(Is_rms^2);


delta_B = (V_in*D1)/(Np*Ac*f);
P_fe = Vc*Kc*(f^alpha)*(delta_B^beta);  %Simple Steinmetz GSE power loss equation 
P_loss = P_cu +P_cu_s+P_fe;







end