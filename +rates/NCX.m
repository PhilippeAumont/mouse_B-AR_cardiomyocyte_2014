function I_NaCa = NCX(Na_i, Ca_i, V, p)
  k_NaCa = 292.8;     %pA/pF
  K_m_Na = 87500;     %uM
  K_m_Ca = 1380;      %uM
  k_sat = 0.1;
  n = 0.35;           %Voltage dependence of NCX

  I_NaCa = k_NaCa*inv((K_m_Na^3+p.Na_o^3)*(K_m_Ca+p.Ca_o)*(1+k_sat*e^((n-1)*V*p.F/(p.R*p.T))))...
            *(e^((n*V*p.F)/(p.R*p.T))*Na_i^3*p.Ca_o - e^((n-1)*V*p.F/(p.R*p.T))*p.Na_o^3*Ca_i);

endfunction
