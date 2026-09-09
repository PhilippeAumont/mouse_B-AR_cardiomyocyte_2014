function dxdt = modelLIH_RA2019(X,p)

%Unpack sate vector
Aeff    = X(1);
H      = X(2)/1e6;
pH     = X(3);
K      = X(4)/1e6;
Na     = X(5)/1e6;
Cl     = X(6)/1e6;
Ca_T   = X(7)/1e6;
Ca_F   = X(8)/1e6;
%{
Ca_C   = X(9)/1e6;
Na_C   = X(10)/1e6;
K_C    = X(11)/1e6;
%}
Ca_C = 1e-7;
Na_C = 0.01;
K_C = 0.145;

B = 0.081;
CAX_Ca = 1;
CAX_H = 3;
CLC_Cl = 2;
CLC_H = 1;
CLC_type = 'WT';
F = 96485;  %C/mol
NA = 6.02e23; %u/mol
NCX_Ca = 3;
NCX_Na = 1;
NHE_H = 1;
NHE_Na = 1;
P_Ca = 1.49e-10;  %cm/ms
P_Cl = 1.2e-8;  %cm/ms
P_H = 6e-8; %cm/ms
P_K = 7.1e-10; %cm/ms
P_Na = 9.6e-10; %cm/ms
R = 0.34;
RTF = 25.690;
beta_pH = 0.04;
cap = 1.4527e-17;
cap_0 = 1e-6;       %F/cm2

N_Lys = 1;
init_V = 1.6464e-16*N_Lys;    %L
S = 1.4527e-8*N_Lys;    %cm2
%
N_CAX = 10*N_Lys;
N_CLC = 300*N_Lys;
N_NCX = 10*N_Lys;
N_NHE = 0*N_Lys;
N_VATP = 300*N_Lys;

pH_C = 7.2;
psi_in = 0;   %mV
psi_out = -50;    %mV
q = 2.2;
r = 0.1;
tau_act = 1000;   %ms
tau_deact = 250;    %ms
Cl_C = 0.01;    %M

#######################################3

%Membrane Potential
psi     = (F/(cap_0*S)*1000)*init_V*(H + K + Na - Cl + 2*Ca_T - B);

%Modified Cytoplasmic Surface Concentrations
pH_C0   = (pH_C+psi_out/(RTF*2.3));
K_C0    = K_C*exp(-psi_out/RTF);
Na_C0   = Na_C*exp(-psi_out/RTF);
Cl_C0   = Cl_C*exp(psi_out/RTF);
Ca_F_C0 = Ca_C*exp(-2*psi_out/RTF);

%Modified Luminal Surface Concentrations
pH_L0   = (pH+psi_in/(RTF*2.3));
K_L0    = K*exp(-psi_in/RTF);
Na_L0   = Na*exp(-psi_in/RTF);
Cl_L0   = Cl*exp(psi_in/RTF);
Ca_F_L0 = Ca_F*exp(-2*psi_in/RTF);

delta_pH    = pH_C0-pH_L0;

%Treatment of singular terms for passive ion flux
if abs(psi) > 0.01
    gg      =  psi / (1 - exp (- psi / RTF)) / RTF;
    gg_Ca   = 2*psi/(1-exp(-2*psi/RTF))/RTF;

else
    gg      =  1 / 1 - (psi / RTF)/2 + (psi / RTF)^2/6 - (psi / RTF)^3 / 24 + (psi / RTF) ^ 4 / 120;
    gg_Ca   = 1/(1 - (psi/RTF) + (2/3)*(psi/RTF)^2 - (1/3)*(psi/RTF)^3 +(2/15)*(psi/RTF)^4);

end

%V-ATPase performance
V_ATPASE    = lys.find_VATP_rate(p.v_flux,psi,pH);
J_VATPASE   = N_VATP*real(V_ATPASE);

%ClC-7 Antiporter {H out, Cl in}
CLC_mu      = (CLC_H + CLC_Cl)*psi + RTF*(CLC_H*2.3*delta_pH + CLC_Cl*log(Cl_C0/Cl_L0));
%Switching function
x           = 0.5 + 0.5*tanh((CLC_mu + 250)/75);
%Activity
A           = 0.3*x + 1.5E-5*(1-x)*CLC_mu^2;

if strcmp(CLC_type,'fast')
    J_CLC    = N_CLC*A*CLC_mu;
    dAeffdt     = 0;
elseif strcmp(CLC_type,'WT')
    if A < Aeff
        tau = tau_deact;
    else
        tau = tau_act;
    end

    dAeffdt   = (1/tau)*(A - Aeff);
    J_CLC    = N_CLC*Aeff*CLC_mu;

end

%CAX Antiporter {H out, Ca in}
CAX_mu   = (CAX_H - 2*CAX_Ca)*psi + RTF*(CAX_H*2.3*delta_pH + CAX_Ca/2*log(Ca_F_L0/Ca_F_C0));
J_CAX    = N_CAX*CAX_mu/1000;

%Passive flux
J_H    = P_H*S*(10^(-pH_C0)*exp(-psi/RTF)-10^(-pH_L0))*gg*NA/1000;
J_K    = P_K*S*(K_C0*exp(-psi/RTF)-K_L0)*gg*NA/1000;
J_Na   = P_Na*S*(Na_C0*exp(-psi/RTF)-Na_L0)*gg*NA/1000;
J_Cl_unc   = P_Cl*S*(Cl_C0-Cl_L0*exp(-psi/RTF))*gg*NA/1000;
J_Ca   = P_Ca*S*(Ca_F_C0*exp(-2*psi/RTF)-Ca_F_L0)*gg_Ca*NA/1000;

%TRPML1 channel
y = 0.5 - 0.5*tanh(psi + 40);
P_trpml1 = 3.88e-9*(y*abs(psi) + (1-y)*(abs(psi + 40)^3)/(pH_L0^q));
J_Ca_trpml1 = P_trpml1/1000*S*(Ca_F_C0*exp(-2*psi/RTF)-Ca_F_L0)*gg_Ca*NA/1000;

###############################################


%Time Dependent Quantities

dNHdt   = J_H + (J_VATPASE) - (CLC_H*J_CLC) - (CAX_H*J_CAX);
dHdt  = dNHdt/(p.NA*p.V_lys_L); %uM

dpHdt   = (-dNHdt/init_V/NA)/beta_pH;

dNKdt   = J_K;
dKdt = dNKdt/(p.NA*p.V_lys_L);%uM

dNNadt  = J_Na;
dNadt = dNNadt/(p.NA*p.V_lys_L);%uM

dNCldt  = J_Cl_unc + (CLC_Cl*J_CLC);
dCldt = dNCldt/(p.NA*p.V_lys_L);%uM

dNCaTdt = J_Ca + (CAX_Ca*J_CAX) + J_Ca_trpml1;
dCaTdt = dNCaTdt/(p.NA*p.V_lys_L);%uM

dNCaFdt = dNCaTdt*r;
dCaFdt = dNCaFdt/(p.NA*p.V_lys_L);%uM

J_N = [J_K; J_Na; J_Cl_unc; J_CLC; J_Ca; J_CAX; J_Ca_trpml1];
J_uM = J_N/(p.NA*(p.V_cyt/1e6)); %Fluxes for cytoplasm


%OUTPUT
dxdt = [dAeffdt; dHdt; dpHdt; dKdt; dNadt; dCldt; dCaTdt ; dCaFdt; J_uM];

