function plot = line_plot(t,X,i)
  %figure('units', 'normalized', 'outerposition', [0 0 1 1]);
  figure('Position', [100, 100, 1000, 1000]);  % [left, bottom, width, height]
  plot(t, X(:,i));
  hold on;
  xline(1.0534);
endfunction
