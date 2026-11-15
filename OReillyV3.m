function [N,L_mh,Ipk_h, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss, LossModel] = OReillyV3(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge, CoreM, CoreS)
%% The OReilly Approach %%
% CITATION 
% Transformers and Inductors for Power Electronics Book by W.G. Hurley
% [N,L_mh,Ipk_h, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss, LossModel] = OReillyV3(12,150000,1.7,0.9,8,0.9,6E-6,2000,1.5,'N87','EFD30')
C           = Constants();

[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM) 
B_max       = B_peak*0.65;
N           = ((V_out+Vd)*(1-D_max))/(V_in*D_max);                                                      % This is the turns ratio referred to the secondary side

%V_ref       = (V_out+Vd)/N;
%V_spike     = Vclamp - (V_in+V_ref);
%N_max       = (V_out+Vd)/(Vds-V_in-V_spike);

L_mh        = (((V_in^2)*(D_max^2))/(P_out))*((1/f)/2);                                                 % L_m must be less than or equal to this formula or it will not operate in DCM/BCM
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
Rdc_s       = ((Ns*MLT*(wire_rs*1E-6))*(1+C.alpha_20*(C.T_max-20)))*1E3;
% AC resistance correction factors.
% These are explicit first-pass assumptions until measured Rac/Rdc values
% are available. Replace them with measured values when possible.
FR_p = 1.20;
FR_s = 1.50;

% If you later measure Rac at the operating frequency, use:
% FR_p = Rp_ac_measured/Rp_dc;
% FR_s = Rs_ac_measured/Rs_dc;


P_cu        = (Rdc*1E-3)*(Ip_rms^2)*FR_p;
P_cu_s      = (Rdc_s*1E-3)*(Is_rms^2)*FR_s;






Ns_l = 1;           % Number of secondary layers
qs = 4;             % This is the number of sections used. 
ins_sin = 1.10;     % single-build enamel (NEMA Type 1)
ins_dou = 1.15;     % heavy-build enamel (NEMA Type 2)
ins_tri = 1.70;     % triple-insulated wire (TIW) typical

bare_wireD = sqrt(4*A_ws/pi)*1E-3;  
out_wireD = (bare_wireD*ins_tri);
T_pl  = ceil(Ns/Ns_l);                   % turns per layer
pitch = C.tape_t+out_wireD;                   % layer to layer pitch

layer_w =(T_pl - 1)*pitch + out_wireD;      % [m] axial layer width

eff_di = pitch - 1.15*bare_wireD + 0.26*out_wireD;

% Use bank winding - we know it works from paper 
L_a = (Ns_l-1)*out_wireD+out_wireD;
Cs =  C.mu_0fm*C.kapton*(out_wireD*(MLT*1E-2)/eff_di)*(L_a/layer_w);
Pcs = 0.5*(Cs)*(V_out^2)*f; % This is the loss associated to the parasitic capacitance 



 %% General Steinmetz Equations (Assumes Sinusoidal Excitation) %% 
 % The flux density ripple can be calculated using Faradays law with Vi as
 % the voltage across the inductor using the T_on time frame
 
% delta_B = (V_in*D1)/(Np*Ac*f);
% P_fe = Vc*Kc*(f^alpha)*(delta_B^beta);  %Simple Steinmetz GSE power loss equation 
% 
% 
%         Vref = V_out/N;
%         S1 = abs(((V_in/(Np*Ac))^alpha)*D1); % I need to confirm values for Ac and that Ae = Ac - Ac here is actually Ae
%         S2 = abs(((Vref/(Np*Ac))^alpha)*D2); 
%         ki_gse = Kc / ( 2^(beta+1)*pi^(alpha-1)*(0.2761 + (1.7061/(alpha + 1.354))));
%         P_feigse = (ki_gse * (delta_B^(beta - alpha)) * (S1 + S2))*Vc;
%         P_loss = P_cu +P_cu_s+P_feigse+Pcs;

%% ========================================================================
%  PREDICTED LOSS-MODEL CURVES
%
%  Models generated:
%   1. GSE baseline:
%          GSE core loss + primary copper + secondary copper
%
%   2. iGSE model:
%          iGSE core loss + primary copper + secondary copper
%
%   3. HV transformer model:
%          iGSE core loss + primary copper + secondary copper
%          + estimated secondary capacitance loss
%
%   4. Proposed full PPU loss model:
%          HV transformer model + MOSFET conduction/switching/Coss/gate
%          + diode conduction/junction-capacitance
%          + leakage/clamp + current-sense + output-capacitor ESR
%          + incremental control loss
%
%  These are predicted loss-power curves versus output voltage.
%  They can later be interpolated onto the measured Vout(t) trajectory
%  and integrated using the measured time vector.
% ========================================================================


%% ---------------- Endpoint loss calculation ----------------------------

% Flux-density excursion at the final operating point.
delta_B = (V_in*D1)/(Np*Ac*f);

% General Steinmetz equation.
P_fe = Vc*Kc*(f^alpha)*(delta_B^beta);

% Your existing iGSE formulation.
Vref = V_out/N;

S1 = abs((V_in/(Np*Ac))^alpha)*D1;
S2 = abs((Vref/(Np*Ac))^alpha)*D2;

ki_iGSE = Kc / ...
    (2^(beta+1)*pi^(alpha-1)* ...
    (0.2761 + 1.7061/(alpha + 1.354)));

P_feigse = ...
    ki_iGSE * ...
    delta_B^(beta-alpha) * ...
    (S1 + S2) * Vc;


%% ---------------- Effective capacitance assumption ---------------------

% Your Cs equation estimates capacitance from the complete winding
% geometry. Therefore it is presently treated as the effective
% capacitance seen by the converter.
Cs_eff = Cs;

% If Cs is later redefined as the capacitance of one secondary section,
% change the previous line to:
%
% Cs_eff = Cs/qs;

Pcs = 0.5*Cs_eff*V_out^2*f;


%% ---------------- Endpoint model hierarchy -----------------------------

P_loss_GSE_end = ...
    P_cu + ...
    P_cu_s + ...
    P_fe;

P_loss_iGSE_end = ...
    P_cu + ...
    P_cu_s + ...
    P_feigse;

P_loss_HV_end = ...
    P_cu + ...
    P_cu_s + ...
    P_feigse + ...
    Pcs;

% Preserve a scalar endpoint output for compatibility.
% It is overwritten below after the complete PPU model is generated.
P_loss = P_loss_HV_end;


%% ========================================================================
%  GENERATE LOSS CURVES THROUGHOUT THE CHARGE
% ========================================================================

% The model is evaluated over output voltage rather than time.
% Later, it can be interpolated onto the measured Vout(t) waveform.

Nv = 2000;

% Avoid exactly zero volts because Toff contains Vout in its denominator.
V_start = 1;

Vvec = linspace(V_start,V_out,Nv);


%% ---------------- BCM switching trajectory -----------------------------

% Assumptions:
%   - constant peak current
%   - constant input voltage
%   - complete demagnetisation every switching cycle
%   - Ts = Ton + Toff

T_on_vec = ...
    (L_mh*Ipk_h/V_in) .* ones(size(Vvec));

T_off_vec = ...
    (N*L_mh*Ipk_h) ./ (Vvec + Vd);

Ts_vec = ...
    T_on_vec + T_off_vec;

f_vec = ...
    1 ./ Ts_vec;

D1_vec = ...
    T_on_vec ./ Ts_vec;

D2_vec = ...
    T_off_vec ./ Ts_vec;


%% ---------------- Current trajectory -----------------------------------

Ip_rms_vec = ...
    Ipk_h .* sqrt(D1_vec/3);

% Secondary peak current referred through the turns ratio.
Is_pk_vec = ...
    (Ipk_h/N) .* ones(size(Vvec));

Is_rms_vec = ...
    Is_pk_vec .* sqrt(D2_vec/3);

Is_avg_vec = ...
    0.5 .* Is_pk_vec .* D2_vec;


%% ---------------- Copper-loss trajectory -------------------------------

% Rdc and Rdc_s are stored numerically in milliohms in your code.
Rp_dc = Rdc*1e-3;
Rs_dc = Rdc_s*1e-3;

P_cu_p_vec = ...
    Rp_dc .* Ip_rms_vec.^2 .* FR_p;

P_cu_s_vec = ...
    Rs_dc .* Is_rms_vec.^2 .* FR_s;

P_cu_total_vec = ...
    P_cu_p_vec + P_cu_s_vec;


%% ---------------- GSE core-loss trajectory -----------------------------

delta_B_vec = ...
    (V_in .* D1_vec) ./ ...
    (Np .* Ac .* f_vec);

P_fe_GSE_vec = ...
    Vc .* Kc .* ...
    f_vec.^alpha .* ...
    delta_B_vec.^beta;


%% ---------------- iGSE core-loss trajectory ----------------------------

% Use your existing iGSE formulation at every output-voltage point.

Vref_vec = ...
    Vvec ./ N;

% To include the diode drop instead, use:
%
% Vref_vec = (Vvec + Vd)./N;

S1_vec = ...
    abs((V_in/(Np*Ac))^alpha) .* D1_vec;

S2_vec = ...
    abs((Vref_vec./(Np*Ac)).^alpha) .* D2_vec;

P_fe_iGSE_vec = ...
    Vc .* ki_iGSE .* ...
    delta_B_vec.^(beta-alpha) .* ...
    (S1_vec + S2_vec);


%% ---------------- Secondary-capacitance trajectory ---------------------

P_Cs_vec = ...
    0.5 .* Cs_eff .* ...
    Vvec.^2 .* ...
    f_vec;


%% ---------------- Converter-loss assumptions ---------------------------
%
% User-specified values:
Rds_on_25 = 0.9;     % MOSFET Rds(on) at 25 C [ohm]
L_lk      = 0.3e-6;     % transformer leakage inductance [H]
%
% Initial engineering assumptions. Replace these with datasheet or measured
% values as they become available.
Rds_temp_factor   = 1.60;      % hot/cold Rds(on) multiplier
Vds_overshoot     = 30;        % residual drain overshoot above reflected voltage [V]
t_fall            = 20e-9;     % MOSFET turn-off current fall time [s]
t_rise            = 10e-9;     % MOSFET turn-on voltage/current overlap time [s]
I_turnon_fraction = 0.05;      % BCM residual turn-on current as fraction of Ipk
Coss_eq           = 80e-12;    % effective MOSFET output capacitance [F]
Qg                = 10e-9;     % MOSFET total gate charge [C]
Vgate             = 6;         % gate-drive voltage [V]
Rsense            = 25e-3;     % primary current-sense resistance [ohm]
R_diode           = 0.50;      % diode dynamic resistance [ohm]
Cj_diode          = 12e-12;     % effective HV diode junction capacitance [F]
ESR_out           = 50e-3;     % output capacitor ESR [ohm]
P_control_inc     = 0.20;      % incremental controller/driver power [W]
%
% Note:
% P_control_inc should not include the static quiescent input power removed
% by your pre-charge baseline-current correction.


%% ---------------- MOSFET losses ----------------------------------------

Rds_on_hot = Rds_on_25 * Rds_temp_factor;

Vds_vec = ...
    V_in + ...
    (Vvec + Vd)./N + ...
    Vds_overshoot;

P_mos_cond_vec = ...
    Rds_on_hot .* ...
    Ip_rms_vec.^2;

I_turnon_vec = ...
    I_turnon_fraction .* Ipk_h .* ...
    ones(size(Vvec));

P_mos_overlap_vec = ...
    0.5 .* Vds_vec .* ...
    (Ipk_h .* t_fall + I_turnon_vec .* t_rise) .* ...
    f_vec;

P_mos_Coss_vec = ...
    0.5 .* Coss_eq .* ...
    Vds_vec.^2 .* ...
    f_vec;

P_gate_vec = ...
    Qg .* Vgate .* ...
    f_vec;

P_switch_total_vec = ...
    P_mos_cond_vec + ...
    P_mos_overlap_vec + ...
    P_mos_Coss_vec + ...
    P_gate_vec;


%% ---------------- Leakage/clamp loss -----------------------------------

% Assumes the leakage energy is dissipated in the clamp/snubber each cycle.
P_lk_vec = ...
    0.5 .* L_lk .* ...
    Ipk_h.^2 .* ...
    f_vec;


%% ---------------- Diode losses -----------------------------------------

P_diode_cond_vec = ...
    Vd .* Is_avg_vec + ...
    R_diode .* Is_rms_vec.^2;

P_diode_cap_vec = ...
    0.5 .* Cj_diode .* ...
    Vvec.^2 .* ...
    f_vec;

P_diode_total_vec = ...
    P_diode_cond_vec + ...
    P_diode_cap_vec;


%% ---------------- Sense, ESR and incremental control --------------------

P_sense_vec = ...
    Rsense .* ...
    Ip_rms_vec.^2;

P_cap_ESR_vec = ...
    ESR_out .* ...
    Is_rms_vec.^2;

P_control_vec = ...
    P_control_inc .* ...
    ones(size(Vvec));


%% ---------------- Create the model curves ------------------------------

% Model 1: classical GSE and AC-corrected copper.
P_loss_GSE_vec = ...
    P_cu_total_vec + ...
    P_fe_GSE_vec;

% Model 2: iGSE and AC-corrected copper.
P_loss_iGSE_vec = ...
    P_cu_total_vec + ...
    P_fe_iGSE_vec;

% Model 3: HV transformer only.
P_loss_HV_vec = ...
    P_cu_total_vec + ...
    P_fe_iGSE_vec + ...
    P_Cs_vec;

% Model 4: proposed complete PPU loss model.
P_loss_PPU_vec = ...
    P_loss_HV_vec + ...
    P_switch_total_vec + ...
    P_lk_vec + ...
    P_diode_total_vec + ...
    P_sense_vec + ...
    P_cap_ESR_vec + ...
    P_control_vec;

% Use the trajectory endpoint for all reported scalar values so the endpoint
% and vector models are calculated from the same current definitions.
P_loss_GSE_end  = P_loss_GSE_vec(end);
P_loss_iGSE_end = P_loss_iGSE_vec(end);
P_loss_HV_end   = P_loss_HV_vec(end);
P_loss_PPU_end  = P_loss_PPU_vec(end);

% Return the full PPU endpoint as the main scalar loss output.
P_loss = P_loss_PPU_end;


%% ---------------- Store all results in one structure -------------------

LossModel = struct();

% Independent variable.
LossModel.Vout = Vvec;

% Predicted operating trajectory.
LossModel.frequency = f_vec;
LossModel.Ton = T_on_vec;
LossModel.Toff = T_off_vec;
LossModel.D1 = D1_vec;
LossModel.D2 = D2_vec;

% Current quantities.
LossModel.Ip_rms = Ip_rms_vec;
LossModel.Is_peak = Is_pk_vec;
LossModel.Is_rms = Is_rms_vec;
LossModel.Is_avg = Is_avg_vec;

% Magnetic quantities.
LossModel.delta_B = delta_B_vec;

% Individual loss terms.
LossModel.P_cu_primary = P_cu_p_vec;
LossModel.P_cu_secondary = P_cu_s_vec;
LossModel.P_cu_total = P_cu_total_vec;

LossModel.P_fe_GSE = P_fe_GSE_vec;
LossModel.P_fe_iGSE = P_fe_iGSE_vec;

LossModel.Cs_estimated = Cs;
LossModel.Cs_effective = Cs_eff;
LossModel.P_Cs = P_Cs_vec;

LossModel.P_mos_cond = P_mos_cond_vec;
LossModel.P_mos_overlap = P_mos_overlap_vec;
LossModel.P_mos_Coss = P_mos_Coss_vec;
LossModel.P_gate = P_gate_vec;
LossModel.P_switch_total = P_switch_total_vec;

LossModel.P_leakage_clamp = P_lk_vec;

LossModel.P_diode_cond = P_diode_cond_vec;
LossModel.P_diode_cap = P_diode_cap_vec;
LossModel.P_diode_total = P_diode_total_vec;

LossModel.P_sense = P_sense_vec;
LossModel.P_cap_ESR = P_cap_ESR_vec;
LossModel.P_control_incremental = P_control_vec;

% Complete model levels.
LossModel.P_loss_GSE = P_loss_GSE_vec;
LossModel.P_loss_iGSE = P_loss_iGSE_vec;
LossModel.P_loss_HV = P_loss_HV_vec;
LossModel.P_loss_PPU = P_loss_PPU_vec;

% Endpoint values.
LossModel.Endpoint.Vout = V_out;
LossModel.Endpoint.frequency = f;

LossModel.Endpoint.P_loss_GSE = P_loss_GSE_end;
LossModel.Endpoint.P_loss_iGSE = P_loss_iGSE_end;
LossModel.Endpoint.P_loss_HV = P_loss_HV_end;
LossModel.Endpoint.P_loss_PPU = P_loss_PPU_end;

LossModel.Endpoint.P_fe_GSE = P_fe_GSE_vec(end);
LossModel.Endpoint.P_fe_iGSE = P_fe_iGSE_vec(end);
LossModel.Endpoint.P_cu_primary = P_cu_p_vec(end);
LossModel.Endpoint.P_cu_secondary = P_cu_s_vec(end);
LossModel.Endpoint.P_Cs = P_Cs_vec(end);

LossModel.Endpoint.P_switch_total = P_switch_total_vec(end);
LossModel.Endpoint.P_leakage_clamp = P_lk_vec(end);
LossModel.Endpoint.P_diode_total = P_diode_total_vec(end);
LossModel.Endpoint.P_sense = P_sense_vec(end);
LossModel.Endpoint.P_cap_ESR = P_cap_ESR_vec(end);
LossModel.Endpoint.P_control_incremental = P_control_vec(end);

% Important model parameters.
LossModel.Parameters.Vin = V_in;
LossModel.Parameters.Vout_final = V_out;
LossModel.Parameters.Cout = C_out;
LossModel.Parameters.turns_ratio = N;
LossModel.Parameters.Lm = L_mh;
LossModel.Parameters.Ipk = Ipk_h;
LossModel.Parameters.Np = Np;
LossModel.Parameters.Ns = Ns;
LossModel.Parameters.CoreMaterial = CoreM;
LossModel.Parameters.CoreSize = CoreS;
LossModel.Parameters.qs = qs;
LossModel.Parameters.FR_p = FR_p;
LossModel.Parameters.FR_s = FR_s;

LossModel.Parameters.Rds_on_25 = Rds_on_25;
LossModel.Parameters.Rds_temp_factor = Rds_temp_factor;
LossModel.Parameters.L_lk = L_lk;
LossModel.Parameters.Vds_overshoot = Vds_overshoot;
LossModel.Parameters.t_fall = t_fall;
LossModel.Parameters.t_rise = t_rise;
LossModel.Parameters.I_turnon_fraction = I_turnon_fraction;
LossModel.Parameters.Coss_eq = Coss_eq;
LossModel.Parameters.Qg = Qg;
LossModel.Parameters.Vgate = Vgate;
LossModel.Parameters.Rsense = Rsense;
LossModel.Parameters.R_diode = R_diode;
LossModel.Parameters.Cj_diode = Cj_diode;
LossModel.Parameters.ESR_out = ESR_out;
LossModel.Parameters.P_control_inc = P_control_inc;


%% ---------------- Automatically save the model data --------------------

CoreM_text = char(string(CoreM));
CoreS_text = char(string(CoreS));

save_name = sprintf( ...
    'OReillyV2_FullPPU_%dS_%gVin_%gVout_%s_%s.mat', ...
    qs,V_in,V_out,CoreM_text,CoreS_text);

% Remove any characters unsuitable for a filename.
save_name = regexprep(save_name,'[^a-zA-Z0-9_.-]','_');

LossModel.FileName = save_name;

save(save_name,'LossModel');

fprintf('\nPredicted loss-model data saved to:\n%s\n',save_name);

fprintf('\nEndpoint predicted losses at %.1f V:\n',V_out);
fprintf('GSE model:         %.6f W\n',P_loss_GSE_end);
fprintf('iGSE model:        %.6f W\n',P_loss_iGSE_end);
fprintf('HV transformer:    %.6f W\n',P_loss_HV_end);
fprintf('Full PPU model:    %.6f W\n',P_loss_PPU_end);





end