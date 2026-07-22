function dS_LCC = LCC_ecav(S, C_ecav, Ca, V, p)
  %Unpack state vector
  o = S(1);
  c1 = S(2);
  c2 = S(3);
  c3 = S(4);
  c4 = S(5);
  cp = S(6);
  i1 = S(7);
  i2 = S(8);
  i3 = S(9);
  op = S(10);
  c1p = S(11);
  c2p = S(12);
  c3p = S(13);
  c4p = S(14);
  cpp = S(15);
  i1p = S(16);
  i2p = S(17);
  i3p = S(18);

  C = C_ecav;

  %Parameters
  f_ecav_ICaL = 0.8;
  I_CaL_tot = 0.0273;         %uM
  K_pc_max = 233.24;          %1/s
  K_pc_half = 10.0;           %uM
  K_pcf = 40000;              %1/s
  K_pcb = 2.4;                %1/s
  k_co = 1000;                %1/s
  k_cop = 4000;               %1/s
  k_oc = 1000;                %1/s

  PP = 0.1;                   %uM , CTE from PP and Inhbitor 1 module

  %Calculations
  a = 0.4*e^((V+15.0)/15.0);
  ap = 0.4*e^((V+15.0+20.0)/15.0);
  B = 0.13*e^(-(V+15.0)/18.0);
  I_CaL_ecav_tot = f_ecav_ICaL*I_CaL_tot*p.V_cell/p.V_ecav;
  y = K_pc_max*Ca/(K_pc_half+Ca);

  c1_c1p = rates.P(c1,C,I_CaL_ecav_tot);
  c1p_c1 = rates.DP(c1p,PP,I_CaL_ecav_tot)*ap^3*k_cop/(a^3*k_co);
  c2_c2p = rates.P(c2,C,I_CaL_ecav_tot);
  c2p_c2 = rates.DP(c2p,PP,I_CaL_ecav_tot)*ap^2*k_cop/(a^2*k_co);
  c3_c3p = rates.P(c3,C,I_CaL_ecav_tot);
  c3p_c3 = rates.DP(c3p,PP,I_CaL_ecav_tot)*ap*k_cop/(a*k_co);
  c4_c4p = rates.P(c4,C,I_CaL_ecav_tot);
  c4p_c4 = rates.DP(c4p,PP,I_CaL_ecav_tot)*k_cop/k_co;
  cp_cpp = rates.P(cp,C,I_CaL_ecav_tot);
  cpp_cp = rates.DP(cpp,PP,I_CaL_ecav_tot)*a*k_cop/(ap*k_co);
  o_op = rates.P(o,C,I_CaL_ecav_tot);
  op_o = rates.DP(op,PP,I_CaL_ecav_tot)*a/ap;
  i1_i1p = rates.P(i1,C,I_CaL_ecav_tot);
  i1p_i1 = rates.DP(i1p,PP,I_CaL_ecav_tot)*a/ap;
  i2_i2p = rates.P(i2,C,I_CaL_ecav_tot);
  i2p_i2 = rates.DP(i2p,PP,I_CaL_ecav_tot);
  i3_i3p = rates.P(i3,C,I_CaL_ecav_tot);
  i3p_i3 = rates.DP(i3p,PP,I_CaL_ecav_tot);

  %Prepare MSM Transition matrix rows
  O =   [-(k_oc+y+0.001*K_pcf+o_op),0,0,0,0,k_co,K_pcb,0.001*a,0,op_o,0,0,0,0,0,0,0,0];
  C1 =  [0,-(4*a+c1_c1p),B,0,0,0,0,0,0,0,c1p_c1,0,0,0,0,0,0,0];
  C2 =  [0,4*a,-(B+3*a+c2_c2p),2*B,0,0,0,0,0,0,0,c2p_c2,0,0,0,0,0,0];
  C3 =  [0,0,3*a,-(2*a+c3_c3p+2*B),3*B,0,0,0,0,0,0,0,c3p_c3,0,0,0,0,0];
  C4 =  [0,0,0,2*a,-((y*K_pcf*k_co/k_oc)+0.01*a*y*k_co/k_oc+a+0.002*K_pcf*k_co/k_oc+c4_c4p+3*B),4*B,0.04*B*K_pcb,0.008*B,4*B*K_pcb,0,0,0,0,c4p_c4,0,0,0,0];
  Cp =  [k_oc,0,0,0,a,-(k_co+cp_cpp+4*B),0,0,0,0,0,0,0,0,cpp_cp,0,0,0];
  I1 =  [y,0,0,0,0.01*a*y*k_co/k_oc,0,-(0.001*K_pcf+K_pcb+i1_i1p+0.04*B*K_pcb),0,0.001*a,0,0,0,0,0,0,i1p_i1,0,0];
  I2 =  [0.001*K_pcf,0,0,0,0.002*K_pcf*k_co/k_oc,0,0,-(0.001*a+y+i2_i2p+0.008*B),K_pcb,0,0,0,0,0,0,0,i2p_i2,0];
  I3 =  [0,0,0,0,y*K_pcf*k_co/k_oc,0,0.001*K_pcf,y,-(4*B*K_pcb+i3_i3p+K_pcb+0.001*a),0,0,0,0,0,0,0,0,i3p_i3];
  Op =  [o_op,0,0,0,0,0,0,0,0,-(y+op_o+0.001*K_pcf+k_oc),0,0,0,0,k_cop,K_pcb,0.001*ap,0];
  C1p = [0,c1_c1p,0,0,0,0,0,0,0,0,-(4*ap+c1p_c1),B,0,0,0,0,0,0];
  C2p = [0,0,c2_c2p,0,0,0,0,0,0,0,4*ap,-(c2p_c2+3*ap+B),2*B,0,0,0,0,0];
  C3p = [0,0,0,c3_c3p,0,0,0,0,0,0,0,3*ap,-(c3p_c3+2*ap+2*B),3*B,0,0,0,0];
  C4p = [0,0,0,0,c4_c4p,0,0,0,0,0,0,0,2*ap,-(c4p_c4+y*K_pcf*k_cop/k_oc+0.01*ap*y*k_cop/k_oc+ap+0.002*K_pcf*k_cop/k_oc+3*B),4*B,0.04*B*K_pcb,0.008*B,4*B*K_pcb];
  Cpp = [0,0,0,0,0,cp_cpp,0,0,0,k_oc,0,0,0,ap,-(cpp_cp+k_cop+4*B),0,0,0];
  I1p = [0,0,0,0,0,0,i1_i1p,0,0,y,0,0,0,0.01*ap*y*k_cop/k_oc,0,-(i1p_i1+0.001*K_pcf+K_pcb+0.04*B*K_pcb),0,0.001*ap];
  I2p = [0,0,0,0,0,0,0,i2_i2p,0,0.001*K_pcf,0,0,0,0.002*K_pcf*k_cop/k_oc,0,0,-(i2p_i2+0.001*ap+y+0.008*B),K_pcb];
  I3p = [0,0,0,0,0,0,0,0,i3_i3p,0,0,0,0,y*K_pcf*k_cop/k_oc,0,0.001*K_pcf,y,-(4*B*K_pcb+i3p_i3+K_pcb+0.001*ap)];

  %Assemble Rows
  Q = [O;C1;C2;C3;C4;Cp;I1;I2;I3;
       Op;C1p;C2p;C3p;C4p;Cpp;I1p;I2p;I3p];

  %Calculate new states
  dS_LCC = Q*S;

endfunction
