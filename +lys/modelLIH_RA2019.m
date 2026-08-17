function dxdt = modelLIH_RA2019(X,p)

%global i Y Z

%= i+1;


%dxdt    = zeros(size(X));
%{
Aeff    = X(1);
NH      = X(2);
pH      = X(3);
NK      = X(4);
NNa     = X(5);
NCl     = X(6);
NCa_T   = X(7);
NCa_F   = X(8);


%Luminal Concentrations
H       = NH/init_V/NA;
K       = NK/init_V/NA;
Na      = NNa/init_V/NA;
Cl      = NCl/init_V/NA;
Ca_F    = NCa_F/init_V/NA;
Ca_T    = NCa_T/init_V/NA;
r       = Ca_F/Ca_T;
%}
Aeff    = X(1);
H      = X(2);
pH     = X(3);
K      = X(4);
Na     = X(5);
Cl     = X(6);
Ca_T   = X(7);
Ca_F   = X(8);
Ca_C   = X(9);
Na_C   = X(10);
K_C    = X(11);

%Membrane Potential
psi     = (p.F/p.C_m_lys/1e6)*p.V_lys_L*(H + K + Na - Cl + 2*Ca_T - p.B);

%Modified Cytoplasmic Surface Concentrations
pH_C0   = (p.pH_C+p.psi_out/(p.RTF*2.3));
K_C0    = K_C*exp(-p.psi_out/p.RTF);
Na_C0   = Na_C*exp(-p.psi_out/p.RTF);
Cl_C0   = p.Cl_C*exp(p.psi_out/p.RTF);
Ca_F_C0 = Ca_C*exp(-2*p.psi_out/p.RTF);

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
    gg      =  1 / 1 - (psi / p.RTF)/2 + (psi / p.RTF)^2/6 - (psi / p.RTF)^3 / 24 + (psi / p.RTF) ^ 4 / 120;
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

if A < Aeff
    tau = p.tau_deact;
else
    tau = p.tau_act;
end

dAeffdt = (1/tau)*(A - Aeff);
J_CLC    = p.N_CLC*Aeff*CLC_mu;

%CAX Antiporter {H out, Ca in}
CAX_mu      = (p.CAX_H - 2*p.CAX_Ca)*psi + p.RTF*(p.CAX_H*2.3*delta_pH + p.CAX_Ca/2*log(Ca_F_L0/Ca_F_C0));
J_CAX    = p.N_CAX*CAX_mu;

%Passive flux
J_H    = p.P_H*p.S*(10^(-pH_C0)*exp(-psi/p.RTF)-10^(-pH_L0))*gg*p.NA*1e3;
J_K    = p.P_K*p.S*(K_C0*exp(-psi/p.RTF)-K_L0)*gg*p.NA*1e3;
J_Na   = p.P_Na*p.S*(Na_C0*exp(-psi/p.RTF)-Na_L0)*gg*p.NA*1e3;
J_Cl_unc   = p.P_Cl*p.S*(Cl_C0-Cl_L0*exp(-psi/p.RTF))*gg*p.NA*1e3;
J_Ca   = p.P_Ca*p.S*(Ca_F_C0*exp(-2*psi/p.RTF)-Ca_F_L0)*gg_Ca*p.NA*1e3;

%TRPML1 channel
y = 0.5 - 0.5*tanh(psi + 40);
P_trpml1 = p.p_trpml1*(y*abs(psi) + (1-y)*(abs(psi + 40)^3)/(pH_L0^p.q));
J_Ca_trpml1 = P_trpml1*p.S*(Ca_F_C0*exp(-2*psi/p.RTF)-Ca_F_L0)*gg_Ca*p.NA*1e3;
J_Ca_trpml1 = 0;
%Time Dependent Quantities

dNHdt   = J_H + (J_VATPASE) - (p.CLC_H*J_CLC) - (p.CAX_H*J_CAX);

dpHdt   = (-dNHdt/p.V_lys_uL/p.NA)/p.beta_pH;

dNKdt   = J_K;

dNNadt  = J_Na;

dNCldt  = J_Cl_unc + (p.CLC_Cl*J_CLC);

dNCaTdt = J_Ca + (p.CAX_Ca*J_CAX) + J_Ca_trpml1;

dNCaFdt = dNCaTdt*p.r;


%Y = [time, H, K, Na, Cl, Ca_F, Ca_T, psi, ...
%    V_ATPASE, J_VATPASE, CLC_mu, A, J_CLC, CAX_mu, J_CAX, ...
%    J_H, J_K, J_Na, J_Cl_unc, J_Ca, J_Ca_trpml1, P_trpml1];

%save_values (Y);

%OUTPUT
%disp([J_H; J_VATPASE; (p.CLC_H*J_CLC); (p.CAX_H*J_CAX)]);
dxdt = [dAeffdt; dNHdt; dpHdt; dNKdt; dNNadt; dNCldt; dNCaTdt ; dNCaFdt; J_K; J_Na; J_Cl_unc; J_CLC; J_Ca; J_CAX; J_Ca_trpml1];

