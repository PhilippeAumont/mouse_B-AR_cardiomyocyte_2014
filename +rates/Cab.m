function I_Cab = Cab(V, Ca_i, p)
  E_CaN = p.R*p.T/(2*p.F) * log(p.Ca_o/Ca_i);
  I_Cab = p.G_Cab*(V - E_CaN);
endfunction
