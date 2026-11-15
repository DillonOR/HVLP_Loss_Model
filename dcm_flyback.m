%% Flyback Capacitor Charger by Dillon O'Reilly %%
% This function will take a list of defined converter specifications to 
% make a flyback transformer that will meet the provided specifications for
% capacitor charging

%   INPUTS:     Input Voltage - V_in, Max Duty Cycle - D_max, Frequency -f, Desired output power
%               P_out, efficiency - n, Output Capacitance - C_out and Desired Output Voltage - V_out 

%   OUTPUTS:    primary side inductance - Lm, peak current - Ipk, secondary side inductance - Ls, the number of
%               turns of wire on the primary side - Np and secondary side Ns and the
%               respective wire sizes - A_wp and A_ws and estimated power loss - P_loss 

% 
% This code is created during a PhD project to design a DC-DC power supply
% for a capacitor charger of a pulsed plasma thruster. 
% It is a heuristic model of multiple design methods meshed into one optimised solution. 

% The code is meant to be a useful tool in deriving the discrete values for the construction of
% a high voltage flyback transformer. Required in the case of a pulsed
% plasma thruster. The choice of magnetic material is at the designers /
% users discretion. In this case, some magnetic materials which the author
% believes are useful for this application are provided. 

% Author: Dillon O'Reilly Contact: dillon.oreilly@community.isunet.edu

% When calling the function use
% [N,L_m,Ipk, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = dcm_flyback(12,150000,1.7,0.9,8,0.9,6E-6,2000,1.5,'N87','EFD30')
% Make sure the function is saved to the working folder path. 


function [N,L_mh,Ipk_h, mu_opt,Ip_rms,Lg,Np,Ns,A_wp,A_ws,P_loss] = dcm_flyback(V_in,f,Vd,D_max,P_out,n,C_out,V_out,t_charge, CoreM, CoreS)
%% Select the Magnetic Material and the Core Sizing Properties %%
% These are two seperate functions I created to call a core material and size from an excel spreadsheet database

 [perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(CoreM)      % Call Core Material Parameters
 [Ac, lc, Wa, Ap, Vc, MLT, Al] = Geometry(CoreS, 'Core_Size')               % Call Core Size Parameters

%% Constants Throughout the Script %% 
C = Constants();

%% Turns Ratio Estimation 
N               =   ((V_out+Vd)*(1-D_max))/(V_in*D_max);  % This is the turns ratio referred to the secondary side
N_lim           =   V_out/V_in;                           % This is a limit on how big N should be.

 
end





% Since in the case of a capacitor charger operating in BCM mode, the
% frequency is variable and therefore the primary side current is
% voltage-dependent, then we must account for the change in peak current
% over the whole charging waveform. We can achieve this by tying the
% instantanous output power of the flyback to the power into the capacitor
% and then integrating from V_init to V_out. Re-arranging then gives us a
% value for the peak current below: 


% To determine the ON time of the switch we need to find the slope of the
% change of current during the ON time, which is given by:


%Having the above gives us the duty ratio for the portion of the time which
%the switch remains on, D1 = T_on/Ts

% The basic concept here is that we are transferring energy from one
% electrical component to another. One is simply emptying and refilling,
% while the other is storing all of the energy. The capacitor. A capacitors
% energy is simply a half times CV squared. If we are to charge the
% capacitor to provide 5-6J of energy into a thruster head and we have a
% desired frequency of operation based on orbitl maneouvering requirements
% of 1Hz. Then we need 5-6 watts on the output. Assuming typical efficiency
% losses, one can say that we should be putting 8-10 watts into the
% converter. If the final energy in the capacitor will be 6 joules after 1
% second, then the transformer must provide that energy in packets over the
% 1 second. 

% Since we know that in the beginning of the cycle, the capacitor
% will charge quickly and rapidly with high currents and towards the end,
% it will charge slowly, then we can say, that most of the energy is
% transferred in the firt few milliseconds of operation. Over the total
% charging cycle, if 6 joules is simply divided by a frequency of
% operation, say 100kHz, then we would need 60uJ per cycle to be delivered
% to the capacitor. If we increase the frequency we reduce the energy per
% cycle transfer.
 

% The total energy delivered in one second would therefore be the amount of
% energy per pulse * the frequency. If we assumed the final frequency: 

% In reality, we will have losses, so that number likely
% needs to be higher. Knowing that E is also equal to 0.5*Lm*Ipk^2. We can
% say that 2*Epulse = Lm*Ipk^2 which would equate to 120uJ. Adding some
% margin for efficiency etc we could round it off to 150uJ per cycle. With
% that said, we can now determine how changes in Lm affect changes in Ipk
% and for a given Lm what value of Ipk would be appropriate. Now if I
% wanted to estimate Lm, I could choose Ipk in accordance with what power
% is available on my input supply. For example, say 2-4 amps. Therefore: 

% This will clearly show that as you decrease the current on the primary
% side, you increase the primary side inductor, likewise, if you increase
% the current, you decrease the inductance. For specific power delivery,
% you can only use small ampage for power supplies. Now using the stored
% energy equation from the energy approach we can approximate an airgap and
% area of airgap. 


% Hurley defines the secondary current as power divided by the off time
% times the output voltage which makes sense. In this case, off time is
% T_off or D2 so we substitute D2 instead of 1-D



% With this information, we can now calculate Kup


% Now we sneed to specify the material properties of the core that we are
% going to use. - I'd like to create a simple function that returns values
% of Kc (material parameter), alpha (a), Beta (b) (material constants) and B_sat based of the material I choose. So I
% could say material N87 in the function and it returns Kc, a, b and B_sat
% values associeated to N87. 
 

 
 % Now a wire size needs to be selected and a wire resistivity needs to be
 % associated to it in order to calculate Rdc - the thernmal loss due to
 % wire heat
 

 % On the ground hc can be assumed a 10W/m but in space there is no
 % convective heating, therefore this number will be higher. I assume 30W/m
 % Assuming an efficiency of 90%, I can estimate thermal change due to
 % power loss 
 
 %P_est_loss = P_out*n;
 
 %delta_T = R_theta*P_est_loss;
 % Better thermal calilcations can be estimated by using hc =
 % 1.42*((delta_T/H)^0.25); where H is the height of the trafo
 % Furthermore, thermal resistance of the core is sometimes given by the
 % manufacturer and can be used to calculate the max dissioation

% Formulas for primary side inductance of a coupled inductor operating in
% DCM/BCM vary slightly from textbook to textbook and from datasheet to datasheet. 
% The general physics is the same and usually each formula presented is
% simply a different way of describing the same fundamental relationships.
% Said another way, each formula presented is simply an algebraically re-arranged 
% version of the same fundamental mathematical relationships. The reason for re-arranging is usually to view the circuit analysis 
% from a perspective that the designer is most interested in. Whether thats
% power, current, frequency, resistance or potentially capacitance. Hurley
% does not include efficiency in his calculation of inductance. 
