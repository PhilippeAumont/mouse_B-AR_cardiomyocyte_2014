function dS = PLB-SERCA(S)
  %Unpack state vector
  f_cyt_PLB_p
  %Parameters
  Km_up_np = 0.41;
  Km_up_p = 0.31;
  k_PLB_PKA = 0.108917;
  K_PLB_PKA = 4.90970;
  k_PLB_PP1 = 4.41956e-2;
  K_PLB_PP1 = 1.69376e-2;

  %Calculations
  Km_up = Km_up_np*(1-f_cyt_PLB_p) + Km_up_p*f_cyt_PLB_p;
  df_cyt_PLB_p = k_PLB_PKA*C_cyt*(1-f_cyt_PLB_p)/(K_PLB_PKA+(1-f_cyt_PLB_p))...
                - k_PLB_PP1*PP1_cyt*f_cyt_PLB_p/(K_PLB_PP1+f_cyt_PLB_p);


  %Pack outputs
endfunction
