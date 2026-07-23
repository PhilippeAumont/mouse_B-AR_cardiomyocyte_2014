function dS = Concentrations(S, J, I, p)
  %Unpack state vector

  %Calcium
  B_i = inv(1 + p.CMDN_tot*p.Km_CMDN/(Km_CMDN+Ca_i)^2);
  B_ss = inv(1 + p.CMDN_tot*p.Km_CMDN/(Km_CMDN+Ca_ss)^2);
  B_JSR = inv(1 + p.CSQN_tot*p.Km_CSQN/(Km_CSQN+Ca_JSR)^2);

  dCa_i = B_i*(J_leak + J_xfer - J_up - J_trpn - (I_Cab - 2*I_NaCa + I_pCa...
                + I_cav_CaL)*(p.A_cap*p.C_m/(2*p.V_cyt*p.F)));
  dCa_ss = B_ss*(J_rel*p.V_JSR/p.V_ss - J_xfer*p.V_cyt/p.V_ss...
          - I_ecav_CaL*p.A_cap*p.C_m/(2*p.V_ss*F))
  dCa_JSR = B_JSR*(J_tr - J_rel);
  dCa_NSR = (J_up - J_leak)*p.V_cyt/p.V_NSR - J_tr*p.V_JSR/p.V_NSR;

  %Sodium
  dNA_i = -(I_Na + I_Nab + 3*I_NaCa + 3*I_NaK)*p.A_cap*p.C_m/(p.V_cyt*p.F);

  %Potassium
  dK_i = -(I_Kto_f + I_Kto_s + I_Kur + I_Kss + I_K1 + I_Kr + I_Ks - 2*I_NaK)...
            *p.A_cap*p.C_m/(p.V_cyt*p.F);

endfunction
