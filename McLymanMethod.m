%  [N,Lm_er,Ipk_e, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = McLymanMethod(12,150000,1.7,0.9,8,0.9,6E-6,2000,1.5,'N87','EFD30')
function [N,L_mm,Ipk_m,mu_eff,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = McLymanMethod(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge,CoreM,CoreS)

%% McLyman-Based Flyback / Coupled-Inductor Sizing
% n is assumed to be efficiency eta.
% N is the secondary-to-primary turns ratio: N = Ns/Np.
%
% Core geometry assumptions from Geometry():
% Ac  = effective core area [cm^2]
% lc  = magnetic path length [cm]
% Wa  = window area [cm^2]
% Ap  = area product [cm^4]
% Vc  = core volume [cm^3] or matched to getMatProps() core-loss units
% MLT = mean length per turn [cm]
% G   = window/gap dimension used for fringing [cm]
%
% Your Constants() uses:
% C.mu_0 = 4*pi*1e-9, which is effectively H/cm.
% Therefore:
% mu0_SI = C.mu_0*1e2 gives H/m.

C = Constants();

[perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM);
B_max = B_peak*0.65;

eta = n;

%% ------------------------------------------------------------------------
%  Core geometry
% -------------------------------------------------------------------------
[Core,Ac,lc,Wa,Ap,Vc,MLT,Al,G] = Geometry(CoreS,'Core_Size');

% Local SI conversions
mu0_SI = C.mu_0*1e2;       % H/m, because C.mu_0 is effectively H/cm
Ac_m2  = Ac*1e-4;          % cm^2 -> m^2
lc_m   = lc*1e-2;          % cm -> m
Wa_m2  = Wa*1e-4;          % cm^2 -> m^2
MLT_m  = MLT*1e-2;         % cm -> m
G_m    = G*1e-2;           % cm -> m

%% ------------------------------------------------------------------------
%  Electrical flyback sizing using McLyman triangular-current energy method
% -------------------------------------------------------------------------

% Secondary-to-primary turns ratio
N = ((V_out + Vd)*(1 - D_max))/(V_in*D_max);

% Switching period and on-time
Ts_input = 1/f;
T_on     = D_max*Ts_input;

% McLyman / energy-transfer peak current
% Ipk = 2PoT / (eta Vin Ton)
% Since Ton = D*T, this becomes:
Ipk_m = (2*P_out)/(eta*V_in*D_max);

% Magnetising inductance from primary current ramp
L_mm = (V_in*T_on)/Ipk_m;

% Equivalent closed form:
% L_mm = eta*V_in^2*D_max^2/(2*P_out*f);

%% ------------------------------------------------------------------------
%  Timing check using the calculated Lm and Ipk
% -------------------------------------------------------------------------

T_on  = (L_mm*Ipk_m)/V_in;
T_off = (N*L_mm*Ipk_m)/(V_out + Vd);

Ts = T_on + T_off;
f  = 1/Ts;

D1 = T_on/Ts;
D2 = T_off/Ts;
Dc = D1 + D2;

%% ------------------------------------------------------------------------
%  RMS currents for triangular DCM/BCM flyback waveforms
% -------------------------------------------------------------------------

Ip_rms  = Ipk_m*sqrt(D1/3);
Ilm_avg = 0.5*Ipk_m*Dc;

% Secondary peak current from output power during secondary conduction
Is_avg = P_out/V_out;
Ispk   = (2*P_out)/(D2*V_out);
Is_rms = Ispk*sqrt(D2/3);

%% ------------------------------------------------------------------------
%  McLyman magnetic sizing
% -------------------------------------------------------------------------

% Stored magnetic energy
E = 0.5*L_mm*Ipk_m^2;

% Primary turns from flux-density constraint
Np_flux = (L_mm*Ipk_m)/(B_max*Ac_m2);

% Use continuous turns for comparison.
% For a buildable design, you can use ceil(Np_flux) instead.
Np = Np_flux;

% Secondary turns
Ns = Np*N;

% Air gap from energy / flux relation, no fringing initially
Lg_m_noFr = mu0_SI*L_mm*Ipk_m^2/(B_max^2*Ac_m2);

% Fringing correction.
% McLyman-style fringing expression:
% Fr = 1 + (lg/sqrt(Ac))*ln(2G/lg)
% Use metres consistently here.
if Lg_m_noFr > 0 && G_m > Lg_m_noFr
    Fr = 1 + (Lg_m_noFr/sqrt(Ac_m2))*log((2*G_m)/Lg_m_noFr);
else
    Fr = 1;
end

% Include finite core reluctance and fringing.
% Magnetic reluctance form:
% L = mu0*Ac*Np^2 / (lg/Fr + lc/mur)
%
% Rearranged:
% lg = Fr*(mu0*Ac*Np^2/L - lc/mur)
Lg_m = Fr*((mu0_SI*Ac_m2*Np^2/L_mm) - (lc_m/perm));

% Avoid negative gap if the finite-core correction dominates
if Lg_m <= 0
    Lg_m = Lg_m_noFr;
    Fr   = 1;
end

% Return gap in cm to remain consistent with your Hurley function output
Lg = Lg_m;

% Effective permeability of the gapped structure
% Returned as mu_eff for compatibility with the mu_opt output slot.
mu_eff = lc_m / ((Lg_m/Fr) + (lc_m/perm));

%% ------------------------------------------------------------------------
%  McLyman core geometry coefficient
% -------------------------------------------------------------------------

% Core geometry coefficient, SI form [m^5]
Kg_core = (Ac_m2^2 * Wa_m2)/MLT_m;

% Thermal dissipation estimate.
% McLyman commonly uses surface-area/thermal-resistance style estimates.
% This keeps your existing thermal structure.
R_theta = 0.06/sqrt(Vc);
P_diss_allowed = C.delta_T/R_theta;

% If C.y = P_fe/P_cu, then:
% P_total = P_cu + P_fe = P_cu*(1 + y)
P_cu_allowed = P_diss_allowed/(1 + C.y);

% Required McLyman Kg for the primary magnetic design.
% Derived from:
% Pcu = Irms^2 * rho * MLT * Np^2 / (Ku*Wa)
% and
% Np = Lm*Ipk/(Bmax*Ac)
%
% SI form:
Kg_required = (C.rho_w * L_mm^2 * Ipk_m^2 * Ip_rms^2) / ...
              (B_max^2 * C.Ku * P_cu_allowed);

%% ------------------------------------------------------------------------
%  Wire areas from copper-loss allocation
% -------------------------------------------------------------------------

% Use the allowed copper loss to determine current density.
% With A_wire = I_rms/J:
%
% Pcu_total = rho*MLT*J*(Np*Ip_rms + Ns*Is_rms)
%
% Therefore:
J = P_cu_allowed/(C.rho_w*MLT_m*(Np*Ip_rms + Ns*Is_rms));   % A/m^2

% Wire copper areas [m^2]
A_wp_m2 = Ip_rms/J;
A_ws_m2 = Is_rms/J;

% Return in mm^2, matching your getwire_resistance() convention
A_wp = A_wp_m2*1e6;
A_ws = A_ws_m2*1e6;

% Window utilization check
Ku_used = (Np*A_wp_m2 + Ns*A_ws_m2)/Wa_m2;

if Ku_used > C.Ku
    warning('McLymanMethod:WindowFillExceeded', ...
        'Window fill exceeded: Ku_used = %.3f, Ku_allowed = %.3f', Ku_used, C.Ku);
end

if Kg_required > Kg_core
    warning('McLymanMethod:KgTooSmall', ...
        'Core Kg may be too small: Kg_required = %.3e m^5, Kg_core = %.3e m^5', Kg_required, Kg_core);
end

%% ------------------------------------------------------------------------
%  Copper loss
% -------------------------------------------------------------------------

[wire_rp, wire_rs] = getwire_resistance(A_wp,A_ws);

% Your getwire_resistance() appears to return resistance in micro-ohm/cm
% or a similar cm-based unit. This preserves your previous convention.
Rdc_p = ((Np*MLT*(wire_rp*1e-6))*(1 + C.alpha_20*(C.T_max - 20)))*1e3;
Rdc_s = ((Ns*MLT*(wire_rs*1e-6))*(1 + C.alpha_20*(C.T_max - 20)))*1e3;

P_cu_p = (Rdc_p*1e-3)*Ip_rms^2;
P_cu_s = (Rdc_s*1e-3)*Is_rms^2;

%% ------------------------------------------------------------------------
%  Core loss
% -------------------------------------------------------------------------

% Flux-density swing using Faraday's law:
% deltaB = Vin*Ton/(Np*Ac)
delta_B = (V_in*T_on)/(Np*Ac_m2);

% Core loss.
% If your Kc is fitted for Vc in cm^3, use Vc directly.
% If your Kc is fitted for SI volume, replace Vc with Vc*1e-6.
P_fe = Vc*Kc*(f^alpha)*(delta_B^beta);

%% ------------------------------------------------------------------------
%  Total loss
% -------------------------------------------------------------------------

P_loss = P_cu_p + P_cu_s + P_fe;

end