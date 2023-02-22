%% Solution system Fluid properties
%clear
Fluid_density_kgm3 = 997;%Kg/m^3
Viscocity_cSt = 0.4733; % cSt
Bulk_Modulus_Pa=2.289e9;%Pascal
Trapped_air=0.005;

Sample_Time=0.001;

%% Hydraulic Motor Properties

AN403992StandardFlowMotorMaxDisp_ccperrev = 9.8*2.1; % CC/Rev
AN403992StandardFlowMotorRatedSpeed_rpm = 6000; %RPM

% AN403993HighFlowMotorMaxDisp_ccperrev = 19; % CC/Rev
% AN403993HighFlowMotorRatedSpeed_rpm = 4000; %RPM

MotorNominalDelP_pa = 100e5*0.01; %Pa, this affects initial regions of results. 
MotorHydOilVolEff = 0.95;

%% PV control Valve Properties

PV_CurrentCmd_mA = [0 300 1500]; % mA base on Pump Calibration tests
PV_FlowArea_m2 = [0 0 0.056]*0.00064516;% m^2
PV_Cd = 0.7; % Discharge Coefficient - 

%% Centrifugal Pump Properties

% AN402280HighFlowPumpRefVelocity_rpm=4700;     %Based on 14GPM Worst Case %6668 6200
% AN402280HighFlowPumpRefDensity_kgperm3=1000;      %kg/m^3
% AN402280HighFlowPumpDelievery_gpm=[0 40 80 120 160 200 240 280 290];% gpm 1lpm =0.2641728747 gpm
% AN402280HighFlowPDiff_psi = [155 142 132 121 108 92 78 62 20];
% 
% AN402280HighFlowPumpBrakePowerNQ_HP = 1.5*(AN402280HighFlowPDiff_psi.*AN402280HighFlowPumpDelievery_gpm)/(1714*5);
% 
% 
% AN402280StdFlowPumpRefVelocity_rpm=6000;     %Based on 14GPM Worst Case %6668 6200
% AN402280StdFlowPumpRefDensity_kgperm3=1000;      %kg/m^3
% AN402280StdFlowPumpDelievery_gpm=[0 40 80 120 160 200 240 280 290];% gpm 1lpm =0.2641728747 gpm
% AN402280StdFlowPDiff_psi = [155 142 132 121 108 92 78 62 20];
% 
% AN402280StdFlowPumpBrakePowerNQ_HP = (AN402280StdFlowPDiff_psi.*AN402280StdFlowPumpDelievery_gpm)/(1714*5);

% 9306C-HM1C 9306CHM1C

% HM1CPumpRefVelocity_rpm=4700;     %Based on 14GPM Worst Case %6668 6200
HM1CPumpRefVelocity_rpm=4500;     %Based on tuning for flow and dP agreement 
HM1CPumpRefDensity_kgperm3=997;      %kg/m^3
% HM1CPumpDelievery_gpm=[0 50 100 150 200 210 211 212];% gpm 1lpm =0.2641728747 gpm
% HM1CPumpPDiff_psi = [130 118 108 88 70 60 40 20];
HM1CPumpDelievery_gpm = xlsread('PumpPerformancePoints.xlsx','TunableTable','A1:A120');
HM1CPumpPDiff_psi = xlsread('PumpPerformancePoints.xlsx','TunableTable','B1:B120');

% load('HM1CPumpDelievery_gpm.mat');
% load('HM1CPumpPDiff_psi.mat');

% HM1CPumpBrakePowerNQ_HP = (HM1CPumpPDiff_psi.*HM1CPumpDelievery_gpm)/(1714*5);
HM1CPumpBrakePowerNQ_HP = (HM1CPumpPDiff_psi.*HM1CPumpDelievery_gpm)/(1714*500); % this seems to affect the latest section of the plots.


%% P-Q-W and N-Q-W tables for pump
% AN402280HighFlowPumpDelievery2D_gpm = [0 40 80 120 160 200 240 280 290];
% AN402280HighFlowAngulerVelocity2D_rpm = [2788 2987 3186 3385 3584 ];
% 
% AN402280HighFlowPDiff2D_psi = [155 142 132 121 108 92 78 62 20;
%                                141 131 120 110 98 85 70 54 20;
%                                128 119 110 98 86 74 62 46 20;
%                                115 106 98 87 76 64 50 35 20;
%                                100 94 87 76 65 54 42 30 20]';
% AN402280HighFlowPumpBrakePowerNQ2D_HP = (AN402280HighFlowPDiff2D_psi.*AN402280HighFlowPumpDelievery2D_gpm)/(1714*0.85);

%% Motor/Pump Shared Parameters

MotorPumpInertia_kgm2 = 0.001929258317282; % kg*m^2; 
MotorPumpDampingCoef = 2.630416384778607e-04; % N*m/(rad/s)

%% SolutionTankSumpValve

SolutionTankSumpValveCommand =[0 0.5 1];
SolutionTankSumpValveOrificeDiameter = 3; %in
SolutionTankSumpValveOpeningArea =[0 (pi*(SolutionTankSumpValveOrificeDiameter^2)/4)/2 (pi*(SolutionTankSumpValveOrificeDiameter^2)/4)]; %in^2

%% Pressure Compensator Valve

PreloadPressure_psi=160; % psi
PistonArea_in2= pi*((0.4372^2)/4); %in^2
SpringRate_lbfperin= (160/150)*28.7;
ValveOpeningVector_in =[0,0.0123,0.0246,0.0369,0.0492,0.0615,0.0738,0.0861,0.0984,0.1107,0.123] ;%[0,0.004572,0.009144,0.013716,0.018288,0.02286,0.027432,0.032004,0.036576,0.041148,0.04572,0.22225,0.238252,0.254254,0.270256,0.286258,0.30226,0.318262,0.334264,0.350266,0.366268,0.38227]*0.393701 ; %in
ValvePassageAreaVector_in2 =[0.094386049,0.084947157,0.074566564,0.063621532,0.052431398,0.041294188,0.030512361,0.020420639,0.011431788,0.004149875,0] ;%[0.046109286,0.048060932,0.04994279,0.051743444,0.053449845,0.055046626,0.05651496,0.057830474,0.058958786,0.059843262,0.06033366,0.06033366,0.057193653,0.051743444,0.045110525,0.037797225,0.03016683,0.022536434,0.015223135,0.008590215,0.003140006,0]*0.155;%in^2
Pascal_To_psi=0.0001450377;

