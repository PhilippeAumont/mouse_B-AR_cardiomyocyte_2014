;1
clear;
%{
This is the main file of a mouse cardiomyocyte model.
It should be in a directory with other directories such as +odes, +plotting and
+rates. This allows it to call scripts/function from those files
Refer to the original paper for details on parameters. I tried to stay as consistent
with the notation as possible, expect when it was too verbose.
The end of this file contains the settings to run simulations and instructions.
Philippe Aumont
Bondarenko VE (2014) A Compartmentalized Mathematical Model of the
B1-Adrenergic Signaling System in Mouse Ventricular Myocytes. PLoS
ONE 9(2): e89113. https://doi.org/10.1371/journal.pone.0089113

Units
- Time: ms (Published model is in s. It was converted to ms)
- Concentration: uM
- Volume: uL
- Voltage: mV
- Currents: pA/pF
- Other rates: mS/uF
%}

%Parameters
% %Cell Parameters
load("parameters.mat");

%  %Load Starting Conditions
X0 = load("X0_OG.mat").ans(:);

%This functions is used to live track progress.
function status = progressBar(t, y, flag, tf)
  persistent lastPercent tstart
  switch flag
    case 'init'
      lastPercent = -1;
      tstart = tic;
      fprintf('Progress:   0%%');
    case ''
      if isempty(t)
        status = 0;
        return;
      end
      currT = t(end);
      percent = floor(100*currT/tf);
      if percent > lastPercent
        lastPercent = percent;
        elapsed = toc(tstart);
        fprintf('\rProgress: %3d%%  (t = %8.2f / %d, elapsed %5.1fs)', ...
                percent, currT, tf, elapsed);
      end
    case 'done'
      fprintf('\rProgress: 100%%  done.                                   \n');
  end
  status = 0;  % return 0 to keep integrating, 1 would stop it
end


%=================================== Set up ====================================
%First run generate_parameters() in the terminal to create the parameters file.
%Then fill out the setting below, and run this main file to run the simulation.

%Stimulation protocols:
p.protocol = "none";
p.NAADP_protocol = "cte";
p.OCaR = 0.00;   %TPC control by OCaR. Scalar controlling J_TPC.
p.L = 0;      %B_AR ligand (ISO) concentration [uM]
p.IBMX = 0;   %PDE inhibitor concentration [uM]

%{
Patch-clamp protocols:
- none: No stimulation
- spike: sigmoidal spike. Needed: start, dur, amp, k
- two_spikes: two single spikes. Needed: start, 2nd_start, dur, amp, k
- train: regular smoothed spikes. Needed: start, period, dur, amp, k

NAADP protocols:
- cte: Constant concentration (can be 0)
- sig: sigmoidal increase
- spike: NAADP spike
Warning: Very short smoothed spikes may not reach full amplitude
Warning: protocol misspell leads to "value on the right hand side of assignment is undefined".
%}
p.stim_start = 1000;%[ms]
p.stim_2nd_start = 10000; %for two_spikes_smooth [ms]
p.stim_period = 100; %For train_smooth [ms]
p.stim_dur = 1;%[ms]
p.stim_amp = -80;%[mV]
p.stim_k = 100000;    %Sharpness of the stimulation current curve

p.NAADP_C = 0;%[uM]

p.NAADP_k1 = 0.01;%Used for sig - if == 1, whole switch happens within 20 ms. If 10, about 1 ms
p.NAADP_t0 = 1000;%[ms], midpoint of the curve
p.NAADP_k2 = 100; %Used for spike
p.NAADP_stim_start = 1000; %[ms]
p.NAADP_stim_dur = 10;%[ms]

p.extra_var = false(); %whether to calculate currents and fluxes, can be time consuming

tspan = [0,5000];
pt_interval = 0.1;%[ms]

%===============================================================================
%Solver
%Might have to adjust tolerances for the lysosome model
abstol_vect = 1e-9*ones(1,155); abstol_vect(79) = 1e-12; abstol_vect(149:155) = 1e-6;
options = odeset('RelTol', 1e-6, 'AbsTol', abstol_vect,"NonNegative", [2:155], 'MaxStep', 1);
options = odeset(options, 'OutputFcn', @(t,y,flag) progressBar(t,y,flag,tspan(2)));

[t,X] = ode15s(@(t,x) odes.odes(t,x,p), tspan, X0, options);

%Time point interpolation
tquery = 0:pt_interval:tspan(2);
X_interp = interp1(t,X,tquery,'spline');  %'spline' assumes continuity, 'pchip' does not.

%Calculate and save the Currents and Fluxes at each timepoint
if (p.extra_var == true())
  I = zeros(size(tquery)(1), 15); J = zeros(size(tquery)(1), 15);
  for k = 1:length(tquery)
    [~, I_k, J_k] = odes.odes(tquery(k), X_interp(k,:).', p);
    I(k,:) = I_k(:).';
    J(k,:) = J_k(:).';
  end
endif

%Plotting - Can be run in the command window
%Line plot:
  %plotting.line_plot(tquery, X_interp, [1]);




