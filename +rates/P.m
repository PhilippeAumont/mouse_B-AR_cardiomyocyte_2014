function k = P(S, C_cav, I_CaL)
  %This script is used to simplify the LCC calculations
  k_ICaL_PKA = 1.74e-2;       %1/s
  K_ICaL_PKA = 0.5;           %uM

  k = k_ICaL_PKA*C_cav/(K_ICaL_PKA+I_CaL*S);
endfunction
