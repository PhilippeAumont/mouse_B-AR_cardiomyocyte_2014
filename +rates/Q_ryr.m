function Q = Q_ryr(Ca_ss)
  k_a_on = 0.006075;  %uM^-4/ms
  k_a_off = 0.07125;   %1/ms
  k_b_on = 0.00405;   %uM^-3/ms
  k_b_off = 0.965;     %1/ms
  k_c_on = 0.009;     %1/ms
  k_c_off = 0.0008;    %1/ms


  Q = [
  -(k_a_off+k_b_on*Ca_ss+k_c_on), k_b_off,    k_a_on*Ca_ss,    k_c_off;
  k_b_on*Ca_ss,                  -(k_b_off), 0,               0;
  k_a_off,                        0,         -(k_a_on*Ca_ss), 0;
  k_c_on,                        0,         0,               -(k_c_off)
  ];

endfunction
