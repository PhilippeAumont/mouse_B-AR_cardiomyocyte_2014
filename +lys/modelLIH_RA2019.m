function dxdt = modelLIH_RA2019(X,p)

%Unpack sate vector
Aeff    = X(1);
H      = X(2)/1e6;
pH     = X(3);
K      = X(4)/1e6;
Na     = X(5)/1e6;
Cl     = X(6)/1e6;
Ca_T   = X(7)/1e6;
%Ca_F   = X(8)/1e6;
Ca_F = p.r*Ca_T;

Ca_C   = X(8)/1e6;
Na_C   = X(9)/1e6;
K_C    = X(10)/1e6;
Ca_md = X(11)/1e6;
NAADP = X(12);


#######################################3

%Membrane Potential
psi     = (p.F_Lys/(p.cap_0*p.S)*1000)*p.init_V*(H + K + Na - Cl + 2*Ca_T - p.B);

%Modified Cytoplasmic Surface Concentrations
pH_C0   = (p.pH_C+p.psi_out/(p.RTF*2.3));
K_C0    = K_C*exp(-p.psi_out/p.RTF);
Na_C0   = Na_C*exp(-p.psi_out/p.RTF);
Cl_C0   = p.Cl_C*exp(p.psi_out/p.RTF);
Ca_F_C0 = Ca_C*exp(-2*p.psi_out/p.RTF);
Ca_md_C0 = Ca_md*exp(-2*p.psi_out/p.RTF);

%Modified Luminal Surface Concentrations
pH_L0   = (pH+p.psi_in/(p.RTF*2.3));
K_L0    = K*exp(-p.psi_in/p.RTF);
Na_L0   = Na*exp(-p.psi_in/p.RTF);
Cl_L0   = Cl*exp(p.psi_in/p.RTF);
Ca_F_L0 = Ca_F*exp(-2*p.psi_in/p.RTF);

delta_pH    = pH_C0-pH_L0;

%Treatment of singular terms for passive ion flux
if abs(psi) > 0.01
    gg      =  psi / (1 - exp (- psi / p.RTF)) / p.RTF;
    gg_Ca   = 2*psi/(1-exp(-2*psi/p.RTF))/p.RTF;

else
    gg      =  1 / (1 - (psi / p.RTF)/2 + (psi / p.RTF)^2/6 - (psi / p.RTF)^3 / 24 + (psi / p.RTF) ^ 4 / 120);
    gg_Ca   = 1/(1 - (psi/p.RTF) + (2/3)*(psi/p.RTF)^2 - (1/3)*(psi/p.RTF)^3 +(2/15)*(psi/p.RTF)^4);

end

%V-ATPase performance
V_ATPASE    = lys.find_VATP_rate(p.v_flux,psi,pH);
J_VATPASE   = p.N_VATP*real(V_ATPASE);

%ClC-7 Antiporter {H out, Cl in}
CLC_mu      = (p.CLC_H + p.CLC_Cl)*psi + p.RTF*(p.CLC_H*2.3*delta_pH + p.CLC_Cl*log(Cl_C0/Cl_L0));
%Switching function
x           = 0.5 + 0.5*tanh((CLC_mu + 250)/75);
%Activity
A           = 0.3*x + 1.5E-5*(1-x)*CLC_mu^2;

if strcmp(p.CLC_type,'fast')
    J_CLC    = p.N_CLC*A*CLC_mu;
    dAeffdt     = 0;
elseif strcmp(p.CLC_type,'WT')
    if A < Aeff
        tau = p.tau_deact;
    else
        tau = p.tau_act;
    end

    dAeffdt   = (1/tau)*(A - Aeff);
    J_CLC    = p.N_CLC*Aeff*CLC_mu/1000;

end

%CAX Antiporter {H out, Ca in}
CAX_mu   = (p.CAX_H - 2*p.CAX_Ca)*psi + p.RTF*(p.CAX_H*2.3*delta_pH + p.CAX_Ca/2*log(Ca_F_L0/Ca_F_C0));
J_CAX    = p.N_CAX*CAX_mu/1000;

%Passive flux
J_H    = p.P_H*p.S*(10^(-pH_C0)*exp(-psi/p.RTF)-10^(-pH_L0))*gg*p.NA/1000;
J_K    = p.P_K*p.S*(K_C0*exp(-psi/p.RTF)-K_L0)*gg*p.NA/1000;
J_Na   = p.P_Na*p.S*(Na_C0*exp(-psi/p.RTF)-Na_L0)*gg*p.NA/1000;
J_Cl_unc   = p.P_Cl*p.S*(Cl_C0-Cl_L0*exp(-psi/p.RTF))*gg*p.NA/1000;
J_Ca   = p.P_Ca*p.S*(Ca_F_C0*exp(-2*psi/p.RTF)-Ca_F_L0)*gg_Ca*p.NA/1000;

%TRPML1 channel - REVIEW TO MAKE RELEVANT
%{
y = 0.5 - 0.5*tanh((psi + 40)/15);
dv = abs(psi+40)/((1 + (abs(psi+40)/200)^8)^(1/8));
P_trpml1 = 3.88e-9*(y*abs(psi) + (1-y)*dv^3/(pH_L0^p.q));
#Adjust below to run
J_Ca_trpml1 = 0.001* P_trpml1/1000*p.S*(Ca_md_C0*exp(-2*psi/p.RTF)-Ca_F_L0)*gg_Ca*p.NA/1000;
%}
f_NAADP = (NAADP^p.n_TPC/(p.Ka_TPC^p.n_TPC + NAADP^p.n_TPC));
J_Ca_trpml1 = p.v_TPC*f_NAADP*(Ca_md_C0*exp(-2*psi/p.RTF)-Ca_F_L0)*gg_Ca*p.NA/1000;
###############################################


%Time Dependent Quantities

dNHdt   = J_H + (J_VATPASE) - (p.CLC_H*J_CLC) - (p.CAX_H*J_CAX);
dHdt  = dNHdt/(p.NA*p.init_V)*1e6; %uM

dpHdt   = (-dNHdt/p.init_V/p.NA)/p.beta_pH;

dNKdt   = J_K;
dKdt = dNKdt/(p.NA*p.init_V)*1e6;%uM

dNNadt  = J_Na;
dNadt = dNNadt/(p.NA*p.init_V)*1e6;%uM

dNCldt  = J_Cl_unc + (p.CLC_Cl*J_CLC);
dCldt = dNCldt/(p.NA*p.init_V)*1e6;%uM

dNCaTdt = J_Ca + (p.CAX_Ca*J_CAX) + J_Ca_trpml1;
dCaTdt = dNCaTdt/(p.NA*p.init_V)*1e6;%uM

dNCaFdt = dNCaTdt*p.r;
dCaFdt = dNCaFdt/(p.NA*p.init_V)*1e6;%uM

J_N = [J_K; J_Na; J_Cl_unc; J_CLC; J_Ca; J_CAX; J_Ca_trpml1];
J_uM = [J_N(1:6)/(p.NA*(p.V_cyt/1e6))*1e6; J_N(7)/(p.NA*(p.V_md/1e6))*1e6];   %Fluxes for cytoplasm
J_uM(7)

%OUTPUT
dxdt = [dAeffdt; dHdt; dpHdt; dKdt; dNadt; dCldt; dCaTdt ; dCaFdt; J_uM];

