;1
%{
This is the main file of a mouse cardiomyocite model.
It should be in a directory with other directories such as +odes, +plotting and
+rates. This allows it to call scripts/function from those files
Refer to the original paper for details on parameters. I tried to stay as consistent
with the notation as possible, expect when it was too verbose.
Philippe Aumont
Bondarenko VE (2014) A Compartmentalized Mathematical Model of the
B1-Adrenergic Signaling System in Mouse Ventricular Myocytes. PLoS
ONE 9(2): e89113. https://doi.org/10.1371/journal.pone.0089113

Units
- Time: s
- Concentration: uM
- Volume: uL
- Voltage: mV
-...

%}

%Parameters
% %Cell Parameters

p.A_cap = 1.534e-4;   %cm2
p.V_cell = 38.00e-4;  %uL
p.V_cyt = 25.84e-6;   %uL
p.V_JSR = 0.12e-6;    %uL
p.V_NSR = 2.098e-6;   %uL
p.V_ss = 1.485e-9;    %uL
p.V_cav = 7.600e-05;  %uL (2% of cell volume)
p.V_ecav = 1.520e-04; %uL (4% of cell volume)

% %Extracellular Ion concentrations
p.K_o = 5400;    %uM
p.Na_o = 140000; %uM
p.Ca_o = 1800;   %uM

% %SR Parameters
p.v1 = 4500;      %1/s %From RyR Module Parameters
p.v2 = 1.74e-2;   %s-1
p.v3 = 306.0;     %1/s %From PLB module
p.t_tr = 0.02;    %s
p.t_xfer = 0.008; %s

% %Calmodulin and Calsequestrin
p.CMDN_tot = 50.0;    %uM
p.CSQN_tot = 15000.0; %uM
p.Km_CMDN = 0.238;     %uM
p.Km_CSQN = 800.0;     %uM

% %Membrane current parameters
p.C_m = 1.0;          %uF/cm2
p.F = 96.5;           %C/mmol
p.T = 298;            %K
p.R = 8.314;          %J/mol*K
p.k_NaCa = 275;       %pA/pF
p.Km_Na = 87500;      %uM
p.Km_Ca = 1380;       %uM
p.k_sat = 0.27;
p.n = 0.35;
p.I_max_pCa = 0.051;
p.Km_pCa = 0.5;       %uM
p.G_Cab = 0.000284;   %mS/uF
p.G_Nab = 0.0063;     %mS/uF
p.G_Kss = 0.0611;     %mS/uF
p.G_Ks = 0.00575;     %mS/uF
p.G_Kr = 0.078;       %mS/uF
p.k_f = 23.761;       %1/s
p.k_b = 36.778;       %1/s
p.GG_ClCa = 10.0;     %mS/uF
p.Km_Cl = 10.0;       %uM
p.E_Cl = -40;         %mV



%Initial Markov State conditions
S_LCC_cav_0 = [
0.320206e-11,   %O
0.973685,       %C1
0.524483e-2,    %C2
0.105944e-4,    %C3
0.951124e-8,    %C4
0.320207e-11,   %Cp
0.308577e-11,   %I1
0.217536e-7,    %I2
0.209641e-7,    %I3
0.562222e-10,   %O-p
0.206347e-1,    %C1-p
0.421668e-3,    %C2-p
0.323128e-5,    %C3-p
0.110051e-7,    %C4-p
0.140555e-10,   %Cp-p
0.541817e-10,   %I1-p
0.100783e-6,    %I2-p
0.970287e-7     %I3-p
];

S_LCC_ecav_0 = [
0.286851e-11,   %O
0.872261,       %C1
0.469850e-2,    %C2
0.949082e-5,    %C3
0.852050e-8,    %C4
0.286852e-11,   %Cp
0.276420e-11,   %I1
0.194870e-7,    %I2
0.187798e-7     %I3
0.328449e-9,    %O-p
0.120548,       %C1-p
0.246338e-2,    %C2-p
0.188771e-4,    %C3-p
0.642918e-7,    %C4-p
0.821123e-10,   %Cp-p
0.316528e-9,    %I1-p
0.588189e-6,    %I2-p
0.566840e-6     %I3-p
];

S_RyR_0 = [
0.854737e-5,    %O1
0.360412e-10,   %O2
0.996216,       %C1
0.961561e-4,    %C2
0.526065e-7,    %O1-p
0.369705e-12,   %O2-p
0.367832e-2,    %C1-p
0.986431e-6     %C2-p
];

S_Na_0 = [
0.367777e-6,    %O
0.161178e-3,    %C1
0.132248e-1,    %C2
0.436222,       %C3
0.153271e-3,    %IF
0.146044e-4,    %I1
0.545874e-7,    %I2
0.125760e-1,    %IC2
0.414822,       %IC3
0.515006e-7,    %O-p
0.225696e-4,    %C1-p
0.185179e-2,    %C2-p
0.610809e-1,    %C3-p
0.214630e-4,    %IF-p
0.217162e-5,    %I1-p
0.301835e-7,    %I2-p
0.176099e-2,    %IC2-p
0.580859e-1     %IC3-p
];

S_IKr_0 = [
0.332600e-3,    %O
0.997365,       %C0
0.135218e-2,    %C1
0.873596e-3,    %C2
0.763767e-4     %I
];


%Initial Conditions
X0 = [
-78.2787,       %membrane potential %mV
0.100157,       %myoplasmic Ca (Ca_i)%uM
0.100157,       %subspace Ca (Ca_ss)%uM
1081.23,        %JSR Ca (Ca_JSR)    %uM
1081.23,        %NSR Ca (Ca_NSR)    %uM
8.66981,        %LTRPNCa            %uM
123.369,        %HTRPNCa            %uM
10508.5,        %Myoplasmic Na      %uM
145400,         %Myoplasmic K       %uM
0.533799e-2,    %a_to_f
0.999945,       %i_to_f
0.713943e-3,    %a_ur
0.996991,       %i_ur
0.713943e-3,    %a_Kss
0.225905,       %f_cav_PLM_p
0.908852,       %f_ecav_IKur
0.713943e-3,    %a_urp
0.996991,       %i_urp
0.252661,       %f_ecav_IKto,f
0.111499e-2,    %a_to_fp
0.999983,       %i_to_fp
0.186637,       %f_cyt_PLB_p
0.364102,       %f_cyt_Tnl_p
0.799452e-3,    %R_cav_PKA
0.626341e-27,   %R_cav_GRK2
0.132189e-2,    %Gs_cav_aGTP
0.180824e-2,    %Gs_cav_By
0.487356e-3,    %Gs_cav_aGDP
0.478002e-1,    %R_ecav_PKA
0.626341e-27,   %R_ecav_GRK2
0.230801e-1,    %Gs_ecav_aGTP
0.237276e-1,    %Gs_ecav_By
0.648475e-3,    %Gs_ecav_aGDP
0.155949e-2,    %R_cyt_PKA
0.626341e-27,   %R_cyt_GRK2
0.331511e-3,    %Gs_cyt_aGTP
0.663570e-3,    %Gs_cyt_By
0.333058e-3,    %Gs_cyt_aGDP
0.000000,       %cAMP_cav_AC56
0.000000,       %cAMP_ecav_AC47
0.000000,       %cAMP_cyt_AC56
0.000000,       %cAMP_cyt_AC47
0.125103e-1,    %PDE3_cav_p
0.580798e-2,    %PDE4_cav_p
0.000000,       %cAMP_cav_PDE2
0.000000,       %cAMP_cav_PDE3
0.000000,       %cAMP_cav_PDE4
0.158226e-1,    %PDE4_ecav_p
0.000000,       %cAMP_ecav_PDE2
0.000000,       %cAMP_ecav_PDE4
0.120998e-2,    %PDE3_cyt_p
0.373102e-2,    %PDE4_cyt_p
0.000000,       %cAMP_cyt_PDE2
0.000000,       %cAMP_cyt_PDE3
0.000000,       %cAMP_cyt_PDE4
7.92317,        %cAMP_cav_PKA
0.299288,       %ARC_cav
0.303358,       %A2RC_cav
0.858440,       %A2R_cav
0.459397e-2,    %C_cav
0.823499,       %PKIC_cav
6.74029,        %cAMP_ecav_PKA
0.653988,       %ARC_ecav
0.132861,       %A2RC_ecav
1.17000,        %A2R_ecav
0.147623,       %C_ecav
1.03338,        %PKIC_ecav
9.32461,        %cAMP_cyt_PKA
0.996350e-1,    %ARC_cyt
0.140099e-1,    %A2RC_cyt
0.273868,       %A2R_cyt
0.665022e-1,    %C_cyt
0.218365,       %PKIC_cyt
0.213571e-1,    %Inhib1_cyt_p
0.253399,       %cAMP_cav
0.507889,       %cAMP_ecav
0.407775,       %cAMP_cyt
0.254152e-11,   %RyR Modulation factor
S_LCC_cav_0,
S_LCC_ecav_0,
S_RyR_0,
S_Na_0,
S_IKr_0
];

% Below needs update
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
p.protocol = "none";
p.L = 0;      %B_AR ligand concentration [uM]
p.IBMX = 0;   %PDE inhibitor concentration [uM]

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
% Warning: protocol misspell leads to "value on the right hand side of assignment is undefined".

p.stim_start = 100;
p.stim_2nd_start = 130;
p.stim_period = 25;
p.stim_dur = 0.5;
p.stim_amp = -80;
p.stim_k = 5000;

%Solver
tspan = 0:0.1:100.0; %THIS MODEL USES S AND NOT MS FOR TIME UNITS
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 1e-2);
%[t,X] = ode15s(@(t,x) odes.odes(t,x,p), tspan, X0, options);
odes.odes(0,X0,p);

%Plotting
%plotting.line_plot(t, X(:,1));     %a:b, includes a but not b






