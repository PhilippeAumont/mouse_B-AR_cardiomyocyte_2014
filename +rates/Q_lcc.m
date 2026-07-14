function Q = Q_lcc(V, Ca)

  K_pc_max = 0.23324;   %1/ms
  K_pc_half = 20.0;     %uM
  K_pcb = 0.0005;       %1/ms

  a = 0.4*e^((V+12.0)/10.0) * (1 + 0.7*e^(-(V+40)^2/10.0) - 0.75*e^(-(V+20.0)^2/400.0))...
      / (1 + 0.12*e^((V+12.0)/10.0));
  b = 0.05*e^(-(V+12.0)/13.0);
  g = K_pc_max*Ca/(K_pc_half+Ca);
  K_pcf = 13.0*(1 - e^(-(V+14.5)^2/100.0));

  %Build Rows
  O = [-(g + 0.001*K_pcf + 4*b), 0, 0, 0, a, K_pcb, 0.001*a, 0];
  C1 = [0, -(4*a),b ,0 ,0 ,0 ,0 ,0 ];
  C2 = [0, 4*a, -(3*a + b), 2*b, 0, 0, 0, 0];
  C3 = [0, 0, 3*a, -(2*a + 2*b), 3*b, 0, 0, 0];
  C4 = [4*b, 0, 0, 2*a, -(g*K_pcf + 0.01*a*g + a + 0.002*K_pcf + 3*b), 0.04*b*K_pcb, 0.008*b, 4*b*K_pcb];
  I1 = [g, 0, 0, 0, 0.01*a*g, -(0.001*K_pcf + K_pcb + 0.04*b*K_pcb), 0, 0.001*a];
  I2 = [0.001*K_pcf, 0, 0, 0, 0.002*K_pcf, 0, -(g + 0.001*a + 0.008*b), K_pcb];
  I3 = [0, 0, 0, 0, g*K_pcf, 0.001*K_pcf, g, -(4*b*K_pcb + K_pcb + 0.001*a)];

  %Assemble rows in matrix
  Q = [
  O,
  C1,
  C2,
  C3,
  C4,
  I1,
  I2,
  I3
  ];

endfunction
