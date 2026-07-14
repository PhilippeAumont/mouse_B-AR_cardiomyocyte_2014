function I_Cl_Ca = Cl_Ca(V, Ca_i)
  G_Cl_Ca = 10.0;   %mS/uF
  K_m_Cl = 10.0;    %uM
  E_Cl = -40.0;     %mV

  O_Cl_Ca = 0.2/(1 + e^(-(V-46.7)/7.8));
  I_Cl_Ca = G_Cl_Ca*O_Cl_Ca*(Ca_i/(Ca_i + K_m_Cl))*(V-E_Cl);
endfunction
