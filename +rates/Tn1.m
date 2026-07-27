function dS = Tn1(S, Ca_i, C_cyt)
  %Unpack state vector
  LTRPNCa = S(1);
  HTRPNCa = S(2);
  f_cyt_Tn1_p = S(3);

  %Parameteres
  LTRPN_tot = 70.0;
  HTRPN_tot = 140.0;
  k_on_htrpn = 2.37;
  k_off_htrpn = 0.032;
  k_on_ltrpn = 32.7;
  k_off_ltrpn_np = 19.6;
  k_off_ltrpn_p = 29.4;
  k_Tn1_PKA = 0.0247254;
  K_Tn1_PKA = 2.71430;
  k_Tn1_PP2A = 0.0865898;
  K_Tn1_PP2A = 0.801420;

  PP2A_cyt = 0.0607843;

  %Calculations
  k_off_ltrpn = k_off_ltrpn_np*(1-f_cyt_Tn1_p) + k_off_ltrpn_p*f_cyt_Tn1_p;
  df_cyt_Tn1_p = k_Tn1_PKA*C_cyt*(1-f_cyt_Tn1_p)/(K_Tn1_PKA+(1-f_cyt_Tn1_p))...
                - k_Tn1_PP2A*PP2A_cyt*f_cyt_Tn1_p/(K_Tn1_PP2A+f_cyt_Tn1_p);

  dLTRPNCa = k_on_ltrpn*Ca_i*(LTRPN_tot-LTRPNCa) - k_off_ltrpn*LTRPNCa;
  dHTRPNCa = k_on_htrpn*Ca_i*(HTRPN_tot-HTRPNCa) - k_off_htrpn*HTRPNCa;

  %Package output
  dS = [k_off_ltrpn; df_cyt_Tn1_p; dLTRPNCa; dHTRPNCa];
endfunction
