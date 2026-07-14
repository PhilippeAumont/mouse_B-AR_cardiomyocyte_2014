function B = buffering(buffer_conc, Ca, buffer_K)
 B = 1/(1 + (buffer_conc*buffer_K/(buffer_K+Ca)^2));
 endfunction
