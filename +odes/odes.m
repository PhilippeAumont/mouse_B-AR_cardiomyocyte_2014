function dxdt = odes(t, X, p)

  %Unpack state vector X
  V = X(1);
  Ca_i = X(2);
  Ca_ss = X(3);
  Ca_JSR = X(4);
  Ca_NSR = X(5);
  LTRPNCa = X(6);
  HTRPNCa = X(7);
  Na_i = X(8);
  K_i = X(9);
  a_to_f = X(10);
  i_to_f = X(11);
  a_ur = X(12);
  i_ur = X(13);
  a_Kss = X(14);
  f_cav_PLM_p = X(15);
  f_ecav_IKur = X(16);
  a_urp = X(17);
  i_urp = X(18);
  f_ecav_IKto_f = X(19);
  a_to_fp = X(20);
  i_to_fp = X(21);
  f_cyt_PLB_p = X(22);
  f_cyt_Tnl_p = X(23);
  R_cav_PKA = X(24);
  R_cav_GRK2 = X(25);
  Gs_cav_aGTP = X(26);
  Gs_cav_By = X(27);
  Gs_cav_aGDP = X(28);
  R_ecav_PKA = X(29);
  R_ecav_GRK2 = X(30);
  Gs_ecav_aGTP = X(31);
  Gs_ecav_By = X(32);
  Gs_ecav_aGDP = X(33);
  R_cyt_PKA = X(34);
  R_cyt_GRK2 = X(35);
  Gs_cyt_aGTP = X(36);
  Gs_cyt_By = X(37);
  Gs_cyt_aGDP = X(38);
  cAMP_cav_AC56 = X(39);
  cAMP_ecav_AC47 = X(40);
  cAMP_cyt_AC56 = X(41);
  cAMP_cyt_AC47 = X(42);
  PDE3_cav_p = X(43);
  PDE4_cav_p = X(44);
  cAMP_cav_PDE2 = X(45);
  cAMP_cav_PDE3 = X(46);
  cAMP_cav_PDE4 = X(47);
  PDE4_ecav_p = X(48);
  cAMP_ecav_PDE2 = X(49);
  cAMP_ecav_PDE4 = X(50);
  PDE3_cyt_p = X(51);
  PDE4_cyt_p = X(52);
  cAMP_cyt_PDE2 = X(53);
  cAMP_cyt_PDE3 = X(54);
  cAMP_cyt_PDE4 = X(55);
  cAMP_cav_PKA = X(56);
  ARC_cav = X(57);
  A2RC_cav = X(58);
  A2R_cav = X(59);
  C_cav = X(60);
  PKIC_cav = X(61);
  cAMP_ecav_PKA = X(62);
  ARC_ecav = X(63);
  A2RC_ecav = X(64);
  A2R_ecav = X(65);
  C_ecav = X(66);
  PKIC_ecav = X(67);
  cAMP_cyt_PKA = X(68);
  ARC_cyt = X(69);
  A2RC_cyt = X(70);
  A2R_cyt = X(71);
  C_cyt = X(72);
  PKIC_cyt = X(73);
  Inhib1_cyt_p_tot = X(74);
  cAMP_cav = X(75);
  cAMP_ecav = X(76);
  cAMP_cyt = X(77);
  P_RyR  = X(78);
  S_LCC_cav = X(79:96);
  S_LCC_ecav = X(97:114);
  S_RyR = X(115:122);
  S_Na = X(123:140);
  S_IKr = X(141:145);


  %Normalize Markov states
  S_LCC_cav = max(S_LCC_cav, 0); S_LCC_cav = S_LCC_cav/sum(S_LCC_cav);
  S_LCC_ecav = max(S_LCC_ecav, 0); S_LCC_ecav = S_LCC_ecav/sum(S_LCC_ecav);
  S_RyR = max(S_RyR, 0); S_RyR = S_RyR/sum(S_RyR);
  S_Na = max(S_Na, 0); S_Na = S_Na/sum(S_Na);
  S_IKr = max(S_IKr, 0); S_IKr = S_IKr/sum(S_IKr);

%============================== Signalling =====================================

  C = [C_cav; C_ecav; C_cyt];
  cAMP = [cAMP_cav; cAMP_ecav; cAMP_cyt];

  %Beta1-Adrenoreceptor module
  S_B_AR = X(24:38);
  dB_AR = rates.B_AR(S_B_AR, C, p);

  %Adenylyl Cyclase Module
  S_AC = [Gs_cav_aGTP;Gs_cav_By;Gs_ecav_aGTP;Gs_ecav_By;Gs_cyt_aGTP;Gs_cyt_By];
  dAC = rates.AC(S_AC, p);

  %Phosphodiesterase Module
  S_PDE = [PDE3_cav_p;PDE4_cav_p;PDE4_ecav_p;PDE3_cyt_p;PDE4_cyt_p];
  dPDE = rates.PDE(S_PDE, C, cAMP, p);

  %cAMP-Protein Kinase A Module
  S_PKA = [X(57:61);X(63:67);X(69:73)];
  dPKA = rates.PKA(S_PKA, C, cAMP, p);

  %Protein Phosphatase & Inhibitor-1 Module
  S_PP = [Inhib1_cyt_p_tot];
  dPP = rates.PP(S_PP, C(3));
  dInhib1_cyt_p_tot = dPP(1);
  PP1_cyt_f = dPP(2);

  %cAMP flux calculation
  S_cAMP = [dAC; dPDE(3:5); dPDE(7:8); dPDE(11:13); dPKA(1); dPKA(7); dPKA(13)];
  dcAMP = rates.cAMP(S_cAMP, cAMP, p);

%========================== Electrochemical Section ============================
  %MSM
  %LCC
  dS_LCC_cav = rates.LCC_cav(S_LCC_cav, C(1), Ca_i, V, p);
  I_cav_CaL = 0.2*(0.3772*S_LCC_cav(1) +  0.7875*S_LCC_cav(10))*(V-52.0);

  dS_LCC_ecav = rates.LCC_ecav(S_LCC_ecav, C(2), Ca_ss, V, p);
  I_ecav_CaL = 0.8*(0.3772*S_LCC_ecav(1) + 0.7875*S_LCC_ecav(10))*(V-52.0);

  I_CaL = I_cav_CaL + I_ecav_CaL;

  %Fast Na
  dS_Na = rates.Fast_Na(S_Na, C(1), V);
  E_Na = (p.R*p.T/p.F) * log((0.9*p.Na_o + 0.1*p.K_o)/(0.9*Na_i + 0.1*K_i));
  I_Na = (14.4*S_Na(1) + 18.0*S_Na(10))*(V-E_Na);

  %RyR
  dS_RyR = rates.RyR(S_RyR, Ca_ss, C(2), p);

  %IKr
  E_Kr = (p.R*p.T/p.F) * log((0.98*p.K_o + 0.02*p.Na_o)/(0.98*K_i + 0.02*Na_i));
  dS_IKr = rates.IKr(S_IKr, V, p);
  I_Kr = p.G_Kr*S_IKr(1)*(V-E_Kr);

  %Others
  %PLB
  S_PLB = [f_cyt_PLB_p; PP1_cyt_f];
  dPLB = rates.PP(S_PLB, C(3));
  Km_up = dPLB(1);
  df_cyt_PLB_p = dPLB(2);

  %Troponin
  S_Tn1 = [LTRPNCa;HTRPNCa;f_cyt_Tnl_p];
  dTn1 = rates.Tn1(S_Tn1, Ca_i, C(3));
  k_off_ltrpn = dTn1(1);
  df_cyt_Tn1_p = dTn1(2);
  dLTRPNCa = dTn1(3);
  dHTRPNCa = dTn1(4);

  %Potassium channels (Kur, Kto,f, K1, Kss)
  S_K = [f_ecav_IKur;a_ur;i_ur;a_urp;i_urp;f_ecav_IKto_f;a_to_f;i_to_f;a_to_fp;i_to_fp;a_Kss];
  dS_K = rates.K_channels(S_K, C(2), V, K_i, p);
  I_Kur = dS_K(1);
  I_Kto_f = dS_K(2);
  I_K1 = dS_K(3);
  I_Kss = dS_K(4);

  %Other channels (NaK, NCX, CaCl, Ca_bckgrd, Na_bckgrd, Ca pump)
  dOthers = rates.Misc_channels(C(1), V, Na_i, Ca_i, E_Na, f_cav_PLM_p, p);
  df_cav_PLM_p = dOthers(1);
  I_NaK = dOthers(2);
  I_pCa = dOthers(3);
  I_NaCa = dOthers(4);
  I_Cab = dOthers(5);
  I_Nab = dOthers(6);
  I_ClCa = dOthers(7);

%================================ Compiling ====================================
  %Fluxes
  J_rel = p.v1*(S_RyR(1)+S_RyR(2)+S_RyR(5)+S_RyR(6))*(Ca_JSR-Ca_ss*P_RyR);
  J_tr = (Ca_NSR - Ca_JSR)/p.t_tr;
  J_xfer = (Ca_ss - Ca_i)/p.t_xfer;
  J_leak = p.v2*(Ca_NSR - Ca_i);
  J_up = p.v3*Ca_i^2/(Km_up^2 + Ca_i^2);
  J_trpn = 2.37*Ca_i*(140.0-HTRPNCa) - 0.032*HTRPNCa + 32.7*Ca_i*(70.0 - LTRPNCa)...
          - k_off_ltrpn*LTRPNCa;  %Params from Tn1 module

  dP_RyR = -0.04*P_RyR - 0.01*(I_ecav_CaL/7.0)*e^-((V+5.0)^2/648.0);

  %Concentration changes and buffers
  Ca = [Ca_i; Ca_ss; Ca_JSR; Ca_NSR];
  J = [J_rel; J_tr; J_xfer; J_leak; J_up; J_trpn];
  I = [I_cav_CaL; I_ecav_CaL; I_Na; I_Kr; I_Kur; I_Kto_f; I_K1; I_Kss; I_NaK; I_pCa; I_NaCa; I_Cab; I_Nab; I_ClCa];
  dC = rates.Concentrations(Ca, J, I, p);

  %membrane potential
  I_stim = 0;
  dV = -(I_CaL + I_pCa + I_NaCa + I_Cab + I_Na + I_Nab + I_NaK + I_Kto_f...
        + I_K1 + I_Kur + I_Kss + I_Kr + I_ClCa - I_stim)/p.C_m;

  %Pack output
  dxdt = [
  dV;
  dC(1:4);
  dLTRPNCa; dHTRPNCa;
  dC(5:6);
  dS_K(5:9);
  f_cav_PLM_p;
  dS_K(10:15);
  df_cyt_PLB_p; df_cyt_Tn1_p;
  dB_AR;
  dAC;
  dPDE;
  dPKA;
  dInhib1_cyt_p_tot;
  dcAMP;
  dP_RyR;
  dS_LCC_cav;
  dS_LCC_ecav;
  dS_RyR;
  dS_Na;
  dS_IKr
  ];

endfunction
