;1

% This is the main file of a mouse cardiomyocite model.
% It should be in a directory with other directories such as +odes, +plotting and
% +rates. This allows it to call scripts/function from those files
% Philippe Aumont
% Model from Bondarenko, V.E. 2014.


%Parameters
% %Cell Parameters

p.A_cap = 1.534e-4;   %cm2
p.V_cell = 38.00e-4;  %uL
p.V_myo = 25.84e-6;   %uL
p.V_JSR = 0.12e-6;    %uL
p.V_NSR = 2.098e-6;   %uL
p.V_ss = 1.485e-9;    %uL
p.V_cav = 7.600e-05;  %uL
p.V_ecav = 1.520e-04; %uL

% %Extracellular Ion concentrations
p.K_o = 5400;    %uM
p.Na_o = 140000; %uM
p.Ca_o = 1800;   %uM

% %SR Parameters
p.v1 = 4.5;       %ms-1
p.v2 = 1.74e-5;   %ms-1
p.v3 = 0.45;      %uM/ms
p.K_m_up = 0.5;   %uM
p.T_tr = 20.0;    %ms
p.T_xfer = 8.0;   %ms

% %LCC parameters
p.G_CaL = 0.1729;     %mS/uF
p.E_Ca_L = 63.0;      %mV
p.I_CaL_max = 7.0;    %pA/pF

% %Buffering
p.LTRPN_tot = 70.0;   %uM
p.HTRPN_tot = 140.0;  %uM
p.k_htrpn_on = 0.00237; %1/uM*ms
p.k_htrpn_off = 3.2e-5;  %ms^-1
p.k_ltrpn_on = 0.0327;  %1/uM*ms
p.k_ltrpn_off = 0.0196;  %ms^-1
p.CMDN_tot = 50.0;    %uM
p.CSQN_tot = 15000.0; %uM
p.K_CMDN = 0.238;     %uM
p.K_CSQN = 800.0;     %uM

% %Membrane current parameters
p.C_m = 1.0;        %uF/cm2
p.F = 96.5;         %C/mmol
p.T = 298;          %K
p.R = 8.314;        %J/mol*K
p.k_NaCa = 292.8;   %pA/pF
p.I_pCa_max = 1.0;  %pA/pF
p.K_m_pCa = 0.5;    %uM
p.G_Cab = 0.000367; %mS/uF
p.G_Na = 13.0;      %mS/uF
p.G_Nab = 0.0026;   %mS/uF
p.G_Kto_f = 0.4067; %Apex = 0.4067 | Septum = 0.0798
p.G_Ks = 0.00575;   %mS/uF
p.G_Kto_s = 0.0;    %Apex = 0.0   | Septum = 0.0629
p.G_Kur = 0.160;    %Apex = 0.160 | Septum = 0.0975
p.G_Kss = 0.050;    %Apex = 0.050 | Septum = 0.0324
p.G_Kr = 0.078;     %mS/uF

% %Temporary Parameters
p.I_stim = 0;

%Initial Markov State conditions
S_LCC_0 = [
0.930308e-18,   %O
0.999876,       %C1
0.124216e-3,    %C2
0.578679e-8,    %C3
0.119816e-12,   %C4
0.497923e-18,   %I1
0.345847e-13,   %I2
0.185106e-13    %I3
];

S_RyR_0 = [
0.149102e-4,    %O1
0.951726e-10,   %O2
0.999817,       %C1
0.167740e-3     %C2
];

S_Na_0 = [
0.713483e-6,    %O
0.279132e-3,    %C1
0.020752,       %C2
0.624646,       %C3
0.153176e-3,    %IF
0.673345e-6,    %I1
0.155787e-8,    %I2
0.0113879,      %IC2
0.342780        %IC3
];

S_IKr_0 = [
0.175298e-3,    %O
0.998159,       %C1
0.992513e-3,    %C2
0.641229e-3,    %C3
0.319129e-4     %I
];


%Initial Conditions
X0 = [
-82.4202,   %membrane potential %mV
0.115001,   %myoplasmic Ca (Ca_i)%uM
0.115001,   %subspace Ca (Ca_ss)%uM
1299.50,    %JSR Ca (Ca_JSR)    %uM
1299.50,    %NSR Ca (Ca_NSR)    %uM
11.2684,    %LTRPNCa            %uM
125.290,    %HTRPNCa            %uM
0.0,        %RyR Modulation factor
14237.1,    %Myoplasmic Na      %uM
143720,     %Myoplasmic K       %uM
0.265563e-2,%a_to_f
0.999977,   %i_to_f
0.262753e-3,%n_Ks
0.417069e-3,%a_to_s
0.998543,   %i_to_s
0.417069e-3,%a_ur
0.998543,   %i_ur
0.417069e-3,%a_Kss
1.0,        %i_Kss
S_LCC_0,
S_RyR_0,
S_Na_0,
S_IKr_0
];

state_vector_labels = [
"Membrane potential [mV]",
"Myoplasmic calcium [uM]",
"Subspace calcium [uM]",
"JSR calcium [uM]",
"NSR calcium [uM]",
"Calcium-bound low affinity troponin [uM]",
"Calcium=bound high affinity troponin [uM]",
"RyR modulation factor",
"Myoplasmic sodium [uM]",
"Myoplasmic potassium [uM]",
"Gating variable a - I_Kto,f",
"Gating variable i - I_Kto,f",
"Gating variable n - I_Ks",
"Gating variable a - I_Kto,s",
"Gating variable i - I_Kto,s",
"Gating variable a - I_Kur",
"Gating variable i - I_Kur",
"Gating variable a - I_Kss",
"Gating variable i - I_Kss",
"LCC MSM - O",
"LCC MSM - C1",
"LCC MSM - C2",
"LCC MSM - C3",
"LCC MSM - C4",
"LCC MSM - I1",
"LCC MSM - I2",
"LCC MSM - I3",
"RyR MSM - O1",
"RyR MSM - O2",
"RyR MSM - C1",
"RyR MSM - C2",
"Na MSM - O",
"Na MSM - C1",
"Na MSM - C2",
"Na MSM - C3",
"Na MSM - IF",
"Na MSM - I1",
"Na MSM - I2",
"Na MSM - IC2",
"Na MSM - IC3",
"IKr MSM - O",
"IKr MSM - C1",
"IKr MSM - C2",
"IKr MSM - C3",
"IKr MSM - I"
];

%Stimulation protocols:
p.protocol = "spike_smooth";
p.L = 0;      %B_AR ligand concentration [uM]
% none: No stimulation
% spike: single stim spike. Needed: start, dur, amp
% spike_smooth: sigmoidal spike. Needed: start, dur, amp, k
% two_spikes_smooth: two single spikes. Needed: start, 2nd_start, dur, amp, k
% train: regular spikes. Needed: start, period, dur, amp
% train_smooth: regular smoothed spikes. Needed: start, period, dur, amp, k
%
%
% k: sharpness of the smoothed curve
% Warning: Non-smooth protocols may cause integration failures.
% Warning: Very short smoothed spikes may not reach full amplitude
% Warning: protocol misspell leades to "value on the right hand side of assignment is undefined".

p.stim_start = 100;
p.stim_2nd_start = 130;
p.stim_period = 25;
p.stim_dur = 0.5;
p.stim_amp = -80;
p.stim_k = 5000;

%Solver
tspan = 0:0.1:300.0;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1e-2);
[t,X] = ode15s(@(t,x) odes.odes(t,x,p), tspan, X0, options);

%Plotting
plotting.line_plot(t, X(:,1));     %a:b, includes a but not b






