function dxdt = odes(t, X, p)

  %Unpack state vector X
  V = X(1);
  Ca_i = X(2);
  Ca_ss = X(3);
  Ca_JSR = X(4);
  Ca_NSR = X(5);
  LTRPNCa = X(6);
  HTRPNCa = X(7);
  P_RyR = X(8);
  Na_i = X(9);
  K_i = X(10);
  a_to_f = X(11);
  i_to_f = X(12);
  n_Ks = X(13);
  a_to_s = X(14);
  i_to_s = X(15);
  a_ur = X(16);
  i_ur = X(17);
  a_Kss = X(18);
  i_Kss = X(19);
  S_LCC = X(20:27);
  S_RyR = X(28:31);
  S_Na = X(32:40);
  S_IKr = X(41:45);

  %Normalize Markov states
  S_LCC = max(S_LCC, 0); S_LCC = S_LCC/sum(S_LCC);
  S_RyR = max(S_RyR, 0); S_RyR = S_RyR/sum(S_RyR);
  S_Na = max(S_Na, 0); S_Na = S_Na/sum(S_Na);
  S_IKr = max(S_IKr, 0); S_IKr = S_IKr/sum(S_IKr);

  %Build Q matrix for MSM
  Q_RyR = rates.Q_ryr(Ca_ss);
  Q_LCC = rates.Q_lcc(V, Ca_ss);
  Q_Na = rates.Q_na(V);
  Q_IKr = rates.Q_ikr(V);

  %Get state change for MSM
  dS_RyR = Q_RyR * S_RyR;
  dS_LCC = Q_LCC * S_LCC;
  dS_Na = Q_Na * S_Na;
  dS_IKr = Q_IKr * S_IKr;

  %Get open probability for MSM
  P_open_RyR = S_RyR(1) + S_RyR(2);
  P_open_LCC = S_LCC(1);
  P_open_Na = S_Na(1);
  P_open_IKr = S_IKr(1);

  %Calculating Factors
  B_i = rates.buffering(p.CMDN_tot, Ca_i, p.K_CMDN);
  B_ss = rates.buffering(p.CMDN_tot, Ca_ss, p.K_CMDN);
  B_JSR = rates.buffering(p.CSQN_tot, Ca_JSR, p.K_CSQN);

  %Get fluxes
  J_rel = p.v1 * P_open_RyR * (Ca_JSR-Ca_ss) * P_RyR;
  J_tr = (Ca_NSR - Ca_JSR)/p.T_tr;
  J_xfer = (Ca_ss - Ca_i)/p.T_xfer;
  J_leak = p.v2*(Ca_NSR - Ca_i);
  J_up = p.v3*Ca_i^2/(p.K_m_up^2 + Ca_i^2);
  J_trpn = p.k_htrpn_on*Ca_i*(p.HTRPN_tot-HTRPNCa) - p.k_htrpn_off*HTRPNCa...
          + p.k_ltrpn_on*Ca_i*(p.LTRPN_tot-LTRPNCa) - p.k_ltrpn_off*LTRPNCa;



  %Get currents

      %Calcium
  I_CaL = p.G_CaL*P_open_LCC*(V-p.E_Ca_L);
  I_pCa = p.I_pCa_max*Ca_i^2 / (p.K_m_pCa^2 + Ca_i^2);
  I_NaCa = rates.NCX(Na_i, Ca_i, V, p);
  I_Cab = rates.Cab(V, Ca_i, p);

      %Sodium
  E_Na = p.R*p.T/p.F * log((0.9*p.Na_o+0.1*p.K_o)/(0.9*Na_i+0.1*K_i));
  I_Na = p.G_Na*P_open_Na*(V-E_Na);
  I_Nab = p.G_Nab*(V-E_Na);
  I_NaK = rates.NaK(V, Na_i, p);

      %Potassium
  E_K = p.R*p.T/p.F * log(p.K_o/K_i);
  I_Kto_f = p.G_Kto_f*a_to_f^3*i_to_f*(V-E_K);
  I_Kto_s = p.G_Kto_s*a_to_s*i_to_s*(V-E_K);
  I_K1 = 0.2938*(p.K_o/(p.K_o + 210.0))*((V-E_K)/(1+e^(0.0896*(V-E_K))));
  I_Ks = p.G_Ks*n_Ks^2*(V-E_K);
  I_Kur = p.G_Kur*a_ur*i_ur*(V-E_K);
  I_Kss = p.G_Kss*a_Kss*i_Kss*(V-E_K);
  I_Kr = P_open_IKr*p.G_Kr*(V - p.R*p.T*inv(p.F)...
          *log((0.98*p.K_o + 0.02*p.Na_o)/(0.98*K_i + 0.02*Na_i)));

      %Chloride
  I_Cl_Ca = rates.Cl_Ca(V, Ca_i);

      %Stimulation Current
  I_stim = odes.I_stim(t, p);

  %Voltage ODE
  dV = -(I_CaL + I_pCa + I_NaCa + I_Cab + I_Na + I_Nab + I_NaK + I_Kto_f...
        +I_Kto_s + I_K1 + I_Ks + I_Kur + I_Kss + I_Kr + I_Cl_Ca + I_stim)/p.C_m;

  %Concentration ODEs
  dCa_i = B_i*(J_leak + J_xfer - J_up - J_trpn...
           - (I_Cab - 2*I_NaCa + I_pCa)*(p.A_cap*p.C_m/2*p.V_myo*p.F));
  dCa_ss = B_ss * (J_rel*p.V_JSR/p.V_ss - J_xfer*p.V_myo/p.V_ss...
             - I_CaL*p.A_cap*p.C_m/(2*p.V_ss*p.F));
  dCa_JSR = B_JSR * (J_tr - J_rel);
  dCa_NSR = (J_up - J_leak)*p.V_myo/p.V_NSR - J_tr*p.V_JSR/p.V_NSR;
  dLTRPNCa = p.k_ltrpn_on*Ca_i*(p.LTRPN_tot-LTRPNCa) - p.k_ltrpn_off*LTRPNCa;
  dHTRPNCa = p.k_htrpn_on*Ca_i*(p.HTRPN_tot-HTRPNCa) - p.k_htrpn_off*HTRPNCa;

  dNa_i = -(I_Na + I_Nab + 3*I_NaCa + 3*I_NaK)*p.A_cap*p.C_m/(p.V_myo*p.F);

  dK_i = -(I_Kto_f + I_Kto_s + I_K1 + I_Ks + I_Kss + I_Kur + I_Kr...
          - 2*I_NaK)*p.A_cap*p.C_m/(p.V_myo*p.F);

  %Other ODEs
  dP_RyR = -0.04*P_RyR - 0.1*(I_CaL/p.I_CaL_max)*e^(-(V-5.0)^2 / 648.0);

      %Fast Transient Outward
  a_a = 0.18064*e^(0.03577*(V+30.0));
  b_a = 0.3956*e^(-0.06237*(V+30.0));
  a_i = 0.000152*e^(-(V+13.5)/7.0)/(0.067083*e^(-(V+33.5)/7.0) + 1);
  b_i = 0.00095*e^((V+33.5)/7.0)/(0.051335*e^((V+33.5)/7.0) + 1);
  da_to_f = a_a*(1 - a_to_f) - b_a*a_to_f;
  di_to_f = a_i*(1 - i_to_f) - b_i*i_to_f;
      %Slow Transient Outward
  a_ss = inv(1+e^(-(V+22.5)/7.7));
  i_ss = inv(1+e^((V+45.2)/5.7));
  T_ta_s = 0.493*e^(-0.0629*V) + 2.058;
  T_ti_s = 270.0 + 1050.0/(1 + e^((V+45.2)/5.7));
  da_to_s = (a_ss - a_to_s)/T_ta_s;
  di_to_s = (i_ss - i_to_s)/T_ti_s;
      %Slowed Delayed Rectifier
  a_n = 0.00000481333*(V+26.5)/(1 - e^(-0.128*(V+26.5)));
  b_n = 0.0000953333*e^(-0.038*(V+26.5));
  dn_Ks = a_n*(1-n_Ks) - b_n*n_Ks;
      %Ultrarapid activating delayed rectifier
  T_aur = 0.493*e^(-0.0629*V) + 2.058;
  T_iur = 1200.0 - 170.0/(1 + e^((V+45.2)/5.7));
  da_ur = (a_ss - a_ur)/T_aur;
  di_ur = (i_ss - i_ur)/T_iur;
      %Noninactivating steady-state
  T_Kss = 39.3*e^(-0.0862*V) + 13.17;
  da_Kss = (a_ss - a_Kss)/T_Kss;
  di_Kss = 0;

  %Pack output
  dxdt = [dV; dCa_i; dCa_ss; dCa_JSR; dCa_NSR; dLTRPNCa; dHTRPNCa;
          dP_RyR; dNa_i; dK_i; da_to_f; di_to_f; dn_Ks; da_to_s;
          di_to_s; da_ur; di_ur; da_Kss; di_Kss; dS_LCC; dS_RyR;
          dS_Na; dS_IKr];
endfunction
