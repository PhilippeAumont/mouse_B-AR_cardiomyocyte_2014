function I_NaK = NaK(V, Na_i, p)
  I_NaK_max = 0.88;   %pA/pF
  K_m_Nai = 21000;    %uM
  K_m_Ko = 1500;      %uM

  s = (e^(p.Na_o/67300) - 1)/7;
  f_NaK = inv(1 + 0.1245*e^(-(0.1*V*p.F)/(p.R*p.T)) + 0.0365*s*e^(-V*p.F/(p.R*p.T)));

  I_NaK = I_NaK_max*f_NaK*inv(1 + (K_m_Nai/Na_i)^1.5)*(p.K_o/(p.K_o + K_m_Ko));

endfunction
