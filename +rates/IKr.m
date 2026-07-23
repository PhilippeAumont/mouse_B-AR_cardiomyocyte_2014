function dS_IKr = IKr(S, V, p)
  %Calculations
  a_a0 = 0.022348*e^(0.01176*V);
  B_a0 = 0.047002*e^(-0.0631*V);
  a_a1 = 0.013733*e^(0.038198*V);
  B_a1 = 0.0000689*e^(-0.04178*V);
  a_ir = 0.090821*e^(0.023391*(V+5.0));
  B_ir = 0.006497*e^(-0.03268*(V+5.0));


  %Prepare MSM Transition matrix

  Q =[
    -(a_ir+B_a1),0,      0,            a_a1,         B_ir;
    0,           -(a_a0),B_a0,         0,            0;
    0,           a_a0,   -(p.k_f+B_a0),p.k_b,        0;
    B_a1,        0,      p.k_f,        -(a_a1+p.k_b),0;
    a_ir,        0,      0,            0,            -(B_ir)
  ];

  %Calculate new states
  dS_IKr = Q*S;

endfunction
