function k = DP(S, PP, I_CaL)
  %This script is used to simplify the LCC calculations
  k_ICaL_PP = 2.325e-4;       %1/s
  K_ICaL_PP = 0.2;            %uM

  k = k_ICaL_PP*PP/(K_ICaL_PP+I_CaL*S);
endfunction
