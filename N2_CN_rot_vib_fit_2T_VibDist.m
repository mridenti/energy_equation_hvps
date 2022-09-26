function [beta, Res] = N2_CN_rot_vib_fit_2T_VibDist

% correct data 

xq = 376:0.01:393;

file_id_absoluto = fopen('absoluto.txt','r');
file_id_resposta = fopen('resposta.txt','r');

DataAbs = fscanf(file_id_absoluto,'%f %f',[2 inf]);
DataAbs = DataAbs';
XDataAbs = DataAbs(1:end,1);
YDataAbs = DataAbs(1:end,2)/max(DataAbs(1:end,2));

DataResp = fscanf(file_id_resposta,'%f %f',[2 inf]);
DataResp = DataResp';
XDataResp = DataResp(1:end,1);
YDataResp = DataResp(1:end,2)/max(DataResp(1:end,2));

pAbs = polyfit(XDataAbs,YDataAbs,2);

pResp = polyfit(XDataResp,YDataResp,8);
% [AjusteResp] = mmq([0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18], XDataResp, YDataResp, 1);
% pResp = zeros(1,19);
% for i = 0:18
%     pResp(i+1) = AjusteResp.A(19-i); 
% end

YAbsInterp = polyval(pAbs,xq);
YRespInterp = polyval(pResp,xq);

calib_q =  YAbsInterp./YRespInterp;

figure(1)
plot(xq, YAbsInterp,'-', XDataAbs, YDataAbs, 'o');
hold on;
plot(xq, YRespInterp,'-', XDataResp, YDataResp, 'o');
hold off;

figure(2)
plot(xq, calib_q,'-');

figure(3);
file_id_entrada = fopen('medida15mm.dat','r');
nol = 17;
for i=1:nol
  fgetl(file_id_entrada);
end
ExpData = fscanf(file_id_entrada,'%f %f',[2 inf]);
ExpData = ExpData';
XData = ExpData(1:end,1);
YData = ExpData(1:end,2);
YDataSel = YData(XData<393 & XData>376);
XDataSel = XData(XData<393 & XData>376);
calib_q_sel = interp1(xq,calib_q,XDataSel); 
YDataSelCorr = YDataSel.*calib_q_sel;
plot(XDataSel,YDataSelCorr,'-b');

close all;

XData = XDataSel;
YData = YDataSelCorr;

% XData correction
XData = XData + 0.08;
%YData = YData - mean(YData(XData>391.7));
%YData = YData - mean(YData(XData<382));
YData = YData - min(YData);
YData = YData/max(YData);
plot(XData,YData);
axis([376 393 0 1]);
hold on;

global N_J;
global h_p c kb e_0 e m_e;
global cte;
global En_B1 En_B2 En_X1 En_X2  DEn_R1 DEn_R2 DEn_P1 DEn_P2 DEn_RQ21 DEn_PQ12 DT_1;
global En_B1_1 En_B2_1 En_X1_1 En_X2_1  DEn_R1_1 DEn_R2_1 DEn_P1_1 DEn_P2_1 DEn_RQ21_1 DEn_PQ12_1;
global CN_En_B1 CN_En_B2 CN_En_X1 CN_En_X2  CN_DEn_R1 CN_DEn_R2 CN_DEn_P1 CN_DEn_P2 CN_DEn_RQ21 CN_DEn_PQ12;
global CN_En_B1_1 CN_En_B2_1 CN_En_X1_1 CN_En_X2_1 CN_DEn_R1_1 CN_DEn_R2_1 CN_DEn_P1_1 CN_DEn_P2_1 CN_DEn_RQ21_1 CN_DEn_PQ12_1;
global CN_En_B1_2 CN_En_B2_2 CN_En_X1_2 CN_En_X2_2 CN_DEn_R1_2 CN_DEn_R2_2 CN_DEn_P1_2 CN_DEn_P2_2 CN_DEn_RQ21_2 CN_DEn_PQ12_2;
global CN_En_B1_3 CN_En_B2_3 CN_En_X1_3 CN_En_X2_3 CN_DEn_R1_3 CN_DEn_R2_3 CN_DEn_P1_3 CN_DEn_P2_3 CN_DEn_RQ21_3 CN_DEn_PQ12_3;
global CN_En_B1_4 CN_En_B2_4 CN_En_X1_4 CN_En_X2_4 CN_DEn_R1_4 CN_DEn_R2_4 CN_DEn_P1_4 CN_DEn_P2_4 CN_DEn_RQ21_4 CN_DEn_PQ12_4;
global CN_En_B1_5 CN_En_B2_5 CN_En_X1_5 CN_En_X2_5 CN_DEn_R1_5 CN_DEn_R2_5 CN_DEn_P1_5 CN_DEn_P2_5 CN_DEn_RQ21_5 CN_DEn_PQ12_5;
global CN_En_B1_6 CN_En_B2_6 CN_En_X1_6 CN_En_X2_6 CN_DEn_R1_6 CN_DEn_R2_6 CN_DEn_P1_6 CN_DEn_P2_6 CN_DEn_RQ21_6 CN_DEn_PQ12_6;
global CN_En_B1_7 CN_En_B2_7 CN_En_X1_7 CN_En_X2_7 CN_DEn_R1_7 CN_DEn_R2_7 CN_DEn_P1_7 CN_DEn_P2_7 CN_DEn_RQ21_7 CN_DEn_PQ12_7;
global CN_En_B1_8 CN_En_B2_8 CN_En_X1_8 CN_En_X2_8 CN_DEn_R1_8 CN_DEn_R2_8 CN_DEn_P1_8 CN_DEn_P2_8 CN_DEn_RQ21_8 CN_DEn_PQ12_8;
global CN_En_B1_9 CN_En_B2_9 CN_En_X1_9 CN_En_X2_9 CN_DEn_R1_9 CN_DEn_R2_9 CN_DEn_P1_9 CN_DEn_P2_9 CN_DEn_RQ21_9 CN_DEn_PQ12_9;
global CN_En_B1_10 CN_En_B2_10 CN_En_X1_10 CN_En_X2_10 CN_DEn_R1_10 CN_DEn_R2_10 CN_DEn_P1_10 CN_DEn_P2_10 CN_DEn_RQ21_10 CN_DEn_PQ12_10;
global CN_DT1 CN_DT2 CN_DT3 CN_DT4 CN_DT5 CN_DT6 CN_DT7 CN_DT8 CN_DT9 CN_DT10;
global delta_lambda_inst;
global T T_1 g_B g_B_1;
global CN_T CN_g_B; 
global N2_Y_B N2_Y_B_4 N2_Y_B_5 N2_Y_B_7 N2_Y_C N2_Y_C_2 N2_Y_C_4 N2_Y_C_1 ...
    N2_Y_B_3 N2_Y_B_2 N2_Y_C_0;
global N2_En_B1_6e N2_En_B2_6e N2_En_B3_6e N2_En_B1_6f N2_En_B2_6f N2_En_B3_6f ...
    N2_En_C1_3e N2_En_C2_3e N2_En_C3_3e N2_En_C1_3f N2_En_C2_3f N2_En_C3_3f... 
    N2_DEn_R1_36e N2_DEn_R2_36e N2_DEn_R3_36e N2_DEn_P1_36e N2_DEn_P2_36e N2_DEn_P3_36e ...
    N2_DEn_Q1_36e N2_DEn_Q2_36e N2_DEn_Q3_36e ...
    N2_DEn_R1_36f N2_DEn_R2_36f N2_DEn_R3_36f N2_DEn_P1_36f N2_DEn_P2_36f N2_DEn_P3_36f ...
    N2_DEn_Q1_36f N2_DEn_Q2_36f N2_DEn_Q3_36f;
global N2_En_B1_5e N2_En_B2_5e N2_En_B3_5e N2_En_C1_2e N2_En_C2_2e N2_En_C3_2e ...
    N2_En_B1_5f N2_En_B2_5f N2_En_B3_5f N2_En_C1_2f N2_En_C2_2f N2_En_C3_2f ...
    N2_DEn_R1_25e N2_DEn_R2_25e N2_DEn_R3_25e N2_DEn_P1_25e N2_DEn_P2_25e N2_DEn_P3_25e ...
    N2_DEn_Q1_25e N2_DEn_Q2_25e N2_DEn_Q3_25e ...
    N2_DEn_R1_25f N2_DEn_R2_25f N2_DEn_R3_25f N2_DEn_P1_25f N2_DEn_P2_25f N2_DEn_P3_25f ...
    N2_DEn_Q1_25f N2_DEn_Q2_25f N2_DEn_Q3_25f;
global N2_En_B1_7e N2_En_B2_7e N2_En_B3_7e N2_En_C1_4e N2_En_C2_4e N2_En_C3_4e ... 
    N2_En_B1_7f N2_En_B2_7f N2_En_B3_7f N2_En_C1_4f N2_En_C2_4f N2_En_C3_4f ...
    N2_DEn_R1_47e N2_DEn_R2_47e N2_DEn_R3_47e N2_DEn_P1_47e N2_DEn_P2_47e N2_DEn_P3_47e ...
    N2_DEn_Q1_47e N2_DEn_Q2_47e N2_DEn_Q3_47e ...
    N2_DEn_R1_47f N2_DEn_R2_47f N2_DEn_R3_47f N2_DEn_P1_47f N2_DEn_P2_47f N2_DEn_P3_47f ...
    N2_DEn_Q1_47f N2_DEn_Q2_47f N2_DEn_Q3_47f;
global N2_En_B1_4e N2_En_B2_4e N2_En_B3_4e N2_En_C1_1e N2_En_C2_1e N2_En_C3_1e ... 
    N2_En_B1_4f N2_En_B2_4f N2_En_B3_4f N2_En_C1_1f N2_En_C2_1f N2_En_C3_1f ... 
    N2_DEn_R1_14e N2_DEn_R2_14e N2_DEn_R3_14e N2_DEn_P1_14e N2_DEn_P2_14e N2_DEn_P3_14e ...
    N2_DEn_Q1_14e N2_DEn_Q2_14e N2_DEn_Q3_14e ...
    N2_DEn_R1_14f N2_DEn_R2_14f N2_DEn_R3_14f N2_DEn_P1_14f N2_DEn_P2_14f N2_DEn_P3_14f ...
    N2_DEn_Q1_14f N2_DEn_Q2_14f N2_DEn_Q3_14f;
global N2_En_B1_3e N2_En_B2_3e N2_En_B3_3e N2_En_C1_0e N2_En_C2_0e N2_En_C3_0e ... 
    N2_En_B1_3f N2_En_B2_3f N2_En_B3_3f N2_En_C1_0f N2_En_C2_0f N2_En_C3_0f ... 
    N2_DEn_R1_03e N2_DEn_R2_03e N2_DEn_R3_03e N2_DEn_P1_03e N2_DEn_P2_03e N2_DEn_P3_03e ...
    N2_DEn_Q1_03e N2_DEn_Q2_03e N2_DEn_Q3_03e ...
    N2_DEn_R1_03f N2_DEn_R2_03f N2_DEn_R3_03f N2_DEn_P1_03f N2_DEn_P2_03f N2_DEn_P3_03f ...
    N2_DEn_Q1_03f N2_DEn_Q2_03f N2_DEn_Q3_03f;
global N2_En_B1_2e N2_En_B2_2e N2_En_B3_2e N2_En_B1_2f N2_En_B2_2f N2_En_B3_2f ... 
    N2_DEn_R1_02e N2_DEn_R2_02e N2_DEn_R3_02e N2_DEn_P1_02e N2_DEn_P2_02e N2_DEn_P3_02e ...
    N2_DEn_Q1_02e N2_DEn_Q2_02e N2_DEn_Q3_02e ...
    N2_DEn_R1_02f N2_DEn_R2_02f N2_DEn_R3_02f N2_DEn_P1_02f N2_DEn_P2_02f N2_DEn_P3_02f ...
    N2_DEn_Q1_02f N2_DEn_Q2_02f N2_DEn_Q3_02f;
global  N2_DEn_R1_13e N2_DEn_R2_13e N2_DEn_R3_13e N2_DEn_P1_13e N2_DEn_P2_13e N2_DEn_P3_13e ...
    N2_DEn_Q1_13e N2_DEn_Q2_13e N2_DEn_Q3_13e ...
    N2_DEn_R1_13f N2_DEn_R2_13f N2_DEn_R3_13f N2_DEn_P1_13f N2_DEn_P2_13f N2_DEn_P3_13f ...
    N2_DEn_Q1_13f N2_DEn_Q2_13f N2_DEn_Q3_13f;
global  N2_DEn_R1_24e N2_DEn_R2_24e N2_DEn_R3_24e N2_DEn_P1_24e N2_DEn_P2_24e N2_DEn_P3_24e ...
    N2_DEn_Q1_24e N2_DEn_Q2_24e N2_DEn_Q3_24e ...
    N2_DEn_R1_24f N2_DEn_R2_24f N2_DEn_R3_24f N2_DEn_P1_24f N2_DEn_P2_24f N2_DEn_P3_24f ...
    N2_DEn_Q1_24f N2_DEn_Q2_24f N2_DEn_Q3_24f;
global V_C_3 V_C_2 V_C_4 V_C_1 V_C_0;
global A_00 A_11;
global A_CN_00 A_CN_11 A_CN_22 A_CN_33 A_CN_44 A_CN_55 A_CN_66 ...
    A_CN_77 A_CN_88 A_CN_99 A_CN_1010;
global A_N2_36 A_N2_25 A_N2_47 A_N2_14 A_N2_03 A_N2_02 A_N2_13 A_N2_24; 
global f_00 f_11 f_22 f_33 f_44 f_55 f_66 f_77 f_88 f_99 f_1010;
global nu_00 nu_11 nu_22 nu_33 nu_44 nu_55 nu_66 nu_77 nu_88 nu_99 nu_1010;

delta_lambda_inst = 5.6137e-02;
N_J = 200;

format long;

% physical constants
e_0 = 8.854187817E-12; %vacuum dielectric permitivity
m_e = 9.10938356E-31; %electron mass
c = 299792458; %light speed
e = 1.6021766208E-19;%electric charge
h_p = 6.626070040E-34;% Planck constant
kb = 1.38064852E-23;

% pre-factor in absorption coefficients
cte = e^2/(4*m_e*c*e_0);

% Spectroscopic constants and parameters
V_C_3 = 35480.326; %pure vibrational energy of v'' = 3, CPi state 
V_C_2 = 33606.182;
V_C_4 = 37261.611;
V_C_1 = 31655.392;
V_C_0 = 29670.942;
% Vibrational Einstein-coefficients 
A_00 = 1.214e7;
A_11 = 4.257E6;
A_CN_00 = 1.4785E7;
A_CN_11 = 1.2285E7;
A_CN_22 = 1.0271E07;
A_CN_33 = 8.6838E06;
A_CN_44 = 7.4763E06;
A_CN_55 = 6.6093E06;
A_CN_66 = 6.0437E06;
A_CN_77 = 5.7402E06;
A_CN_88 = 5.6556E06;
A_CN_99 = 5.7353E06;
A_CN_1010 = 5.9031E06;
A_N2_36 = 2.972E6;
A_N2_25 = 3.081E6;
A_N2_47 = 2.326E6; 
A_N2_14 = 2.375E6;
A_N2_03 = 1.086E6;
A_N2_02 = 3.532E6;
A_N2_13 = 4.885E6;
A_N2_24 = 4.045E6;

%oscillator strengths
f_00 = 3.45E-02;
f_11 = 2.87E-02;
f_22 = 2.41E-02;
f_33 = 2.05E-02;
f_44 = 1.75E-02;
f_55 = 1.53E-02;
f_66 = 1.38E-02;
f_77 = 1.29E-02;
f_88 = 1.27E-02;
f_99 = 1.30E-02;
f_1010 = 1.38E-02;

%band origins
nu_00 = c/387.7E-09;
nu_11 = c/386.5E-09;
nu_22 = c/385.5E-09;
nu_33 = c/384.7E-09;
nu_44 = c/384.2E-09;
nu_55 = c/384.0E-09;
nu_66 = c/384.3E-09;
nu_77 = c/384.9E-09;
nu_88 = c/386.0E-09;
nu_99 = c/387.5E-09;
nu_1010 = c/389.6E-09;

% molecular constants

% N2(+) constants (taken from Journal of Molecular Spectroscopy 203, 1-8 (2000))

% excited state B sigma

% (nu = 0)

T = 25566.062;
B_B = 2.07456;
D_B = 0.6273E-05;
H_B = -0.573E-11;
g_B = 0.0240;

% (nu = 1)

T_1 = 27937.681;
DT_1 = 2371.619;
B_B_1 = 2.05145;
D_B_1 = 0.6459E-05;
H_B_1 = -0.808E-11;
g_B_1 = 0.0202;


% ground state X

% (nu = 0)

B_X = 1.922316;
D_X = 0.5919E-5;
g_X = 0.918E-2;

% (nu = 1)

T_X1 = 2174.746;
B_X_1 = 1.903380;
D_X_1 = 0.5958E-5;
g_X_1 = 0.920E-2;

En_B1 = zeros(N_J,1);
En_B2 = zeros(N_J,1);
En_X1 = zeros(N_J,1);
En_X2 = zeros(N_J,1);
DEn_R1 = zeros(N_J,1);
DEn_R2 = zeros(N_J,1);
DEn_P1 = zeros(N_J,1);
DEn_P2 = zeros(N_J,1);
DEn_RQ21 = zeros(N_J,1);
DEn_PQ12 = zeros(N_J,1);

En_B1_1 = zeros(N_J,1);
En_B2_1 = zeros(N_J,1);
En_X1_1 = zeros(N_J,1);
En_X2_1 = zeros(N_J,1);
DEn_R1_1 = zeros(N_J,1);
DEn_R2_1 = zeros(N_J,1);
DEn_P1_1 = zeros(N_J,1);
DEn_P2_1 = zeros(N_J,1);
DEn_RQ21_1 = zeros(N_J,1);
DEn_PQ12_1 = zeros(N_J,1);

% CN constants (taken from JOURNAL OF MOLECULAR SPECTROSCOPY 156, 327-340 (1992) )

% excited state B sigma

% (nu = 0)

CN_T =  26828.963;%25797.8693;
CN_B_B = 1.95874;%1.958728;
CN_D_B = 6.580e-6;%0.6606E-05;
CN_H_B = 0.0;
CN_L_B = 0.0;
CN_M_B = 0.0;
CN_N_B = 0.0;
CN_g_B = 0.017155;

% (nu = 1)

CN_T_1 = 28952.539; %27921.4652;
CN_B_B_1 = 1.93792;%1.937995;
CN_D_B_1 = 6.690e-6;%0.6654E-05;
CN_H_B_1 = 0.0;
CN_L_B_1 = 0.0;
CN_M_B_1 = 0.0;
CN_N_B_1 = 0.0;
CN_g_B_1 = 0.02;%0.01803;

% (nu = 2)

CN_T_2 = 31035.958; %30004.8966;
CN_B_B_2 = 1.91625;%1.916649;
CN_D_B_2 = 6.820e-6;%0.768E-05;
CN_H_B_2 = 0.0;
CN_L_B_2 = 0.0;
CN_M_B_2 = 0.0;
CN_N_B_2 = 0.0;
CN_g_B_2 = 0.02;%0.01829;

% (nu = 3)

CN_T_3 = 33076.969; %32045.9383;
CN_B_B_3 = 1.893825;%1.894004;
CN_D_B_3 = 6.980e-6;%0.542E-05;
CN_H_B_3 = 0.0;
CN_L_B_3 = 0.0;
CN_M_B_3 = 0.0;
CN_N_B_3 = 0.0;
CN_g_B_3 = 0.02;%0.02406;

% (nu = 4)

CN_T_4 = 35072.973; %34041.9520;
CN_B_B_4 = 1.87045;%1.870330;
CN_D_B_4 = 7.190e-6;%0.638E-05;
CN_H_B_4 = 0.0;
CN_L_B_4 = 0.0;
CN_M_B_4 = 0.0;
CN_N_B_4 = 0.0;
CN_g_B_4 = 0.02;%0.01892;

% (nu = 5)

CN_T_5 = 37021.124; %5990.0220;
CN_B_B_5 = 1.845075;%1.85158;
CN_D_B_5 = 7.430e-6;%9.41E-05;
CN_H_B_5 = 0.0;
CN_L_B_5 = 0.0;
CN_M_B_5 = 0.0;
CN_N_B_5 = 0.0;
CN_g_B_5 = 0.02;%0.0111;

% (nu = 6)

CN_T_6 = 38918.365; %37887.3796;
CN_B_B_6 = 1.81933;%1.819169;
CN_D_B_6 = 7.720e-6;%.606E-05;
CN_H_B_6 = 0.0;
CN_L_B_6 = 0.0;
CN_M_B_6 = 0.0;
CN_N_B_6 = 0.0;
CN_g_B_6 = 0.02;%0.02420;

% (nu = 7)

CN_T_7 = 40761.457; %39730.4754;
CN_B_B_7 = 1.79091;%1.790946;
CN_D_B_7 = 7.927e-6;%1.217E-05;
CN_H_B_7 = -26.9e-12;
CN_L_B_7 = 4.2e-16;
CN_M_B_7 = -6.3e-20;
CN_N_B_7 = 4.6e-25;
CN_g_B_7 = 0.02;%0.00657;

% (nu = 8)

CN_T_8 = 42547.834; %41516.5744;
CN_B_B_8 = 1.76175;%1.762140;
CN_D_B_8 = 8.661e-6;%0.862E-05;
CN_H_B_8 = -18.0e-12;
CN_L_B_8 = 4.2e-16;
CN_M_B_8 = -6.3e-20;
CN_N_B_8 = 4.6e-25;
CN_g_B_8 = 0.02;%0.03192;

% (nu = 9)

CN_T_9 = 44273.969; %43242.9051;
CN_B_B_9 = 1.73047;%1.730417;
CN_D_B_9 = 9.145e-6;%10.086E-06;
CN_H_B_9 = -23.9e-12;
CN_L_B_9 = 4.2e-16;
CN_M_B_9 = -6.3e-20;
CN_N_B_9 = 4.6e-25;
CN_g_B_9 = 0.02;%0.00626;

% (nu = 10)

CN_T_10 = 46511.3138;
CN_B_B_10 = 1.66477;
CN_D_B_10 = 0.81E-05;
CN_H_B_10 = 0.0;
CN_L_B_10 = 0.0;
CN_M_B_10 = 0.0;
CN_N_B_10 = 0.0;
CN_g_B_10 = 0.0244;

CN_DT1 = CN_T_1 - CN_T;
CN_DT2 = CN_T_2 - CN_T;
CN_DT3 = CN_T_3 - CN_T;
CN_DT4 = CN_T_4 - CN_T;
CN_DT5 = CN_T_5 - CN_T;
CN_DT6 = CN_T_6 - CN_T;
CN_DT7 = CN_T_7 - CN_T;
CN_DT8 = CN_T_8 - CN_T;
CN_DT9 = CN_T_9 - CN_T;
CN_DT10 = CN_T_10 - CN_T;

% ground state X

% (nu = 0)
CN_T_X = 1031.124;%0.0;
CN_B_X = 1.891025;%1.89109067;
CN_D_X = 6.410e-6;%0.64094E-5;
CN_H_X = 0.0;%4.9E-12;
CN_g_X = 0.6E-2;%0.725517E-2;

% (nu = 1)

CN_T_X1 = 3073.527;%2042.41851;
CN_B_X_1 = 1.87355;%1.87366616;
CN_D_X_1 = 6.420e-6;%0.64157E-5;
CN_H_X_1 = 0.0;%4.1E-12;
CN_g_X_1 = 0.6E-2;%0.717398E-2;

% (nu = 2)

CN_T_X2 = 5089.623;%4058.54467;
CN_B_X_2 = 1.855725;%1.85618743;
CN_D_X_2 = 6.430e-6 ;%0.64239E-5;
CN_H_X_2 = 0.0;%3.3E-12;
CN_g_X_2 = 0.6E-2;%0.70833E-2;

% (nu = 3)

CN_T_X3 = 7079.4;%6048.33680;
CN_B_X_3 = 1.838675;%1.83865304;
CN_D_X_3 = 6.440e-6;%0.64355E-5;
CN_H_X_3 = 0.0;%3.4E-12;
CN_g_X_3 = 0.6E-2;%0.69819E-2;

% (nu = 4)

CN_T_X4 = 9042.802;%8011.7554;
CN_B_X_4 = 1.820955;%1.82105956;
CN_D_X_4 = 6.450e-6;%0.64403E-5;
CN_H_X_4 = 0.0;%0;
CN_g_X_4 = 0.6E-2;%0.68645E-2;

% (nu = 5)

CN_T_X5 = 10979.769 ;%9948.7499;
CN_B_X_5 = 1.80331;%1.80340517;
CN_D_X_5 = 6.470e-6;%0.658E-5;
CN_H_X_5 = 0.0;%0;
CN_g_X_5 = 0.6E-2;%0.67208E-2;

% (nu = 6)

CN_T_X6 = 12890.32;%11859.2870;
CN_B_X_6 = 1.78584;%1.78568568;
CN_D_X_6 = 6.480e-6;%0.6545E-5;
CN_H_X_6 = 0.0;%0;
CN_g_X_6 = 0.6E-2;%0.65436E-2;

% (nu = 7)

CN_T_X7 = 14774.3715;%13743.3235;
CN_B_X_7 = 1.767898777;%1.76789664;
CN_D_X_7 = 6.500e-6;%0.6207E-5;
CN_H_X_7 = 0.0;%0;
CN_g_X_7 = 0.6E-2;%0.63137E-2;

% (nu = 8)

CN_T_X8 = 16631.8586;%15600.8106;
CN_B_X_8 = 1.750040523;%1.75004005;
CN_D_X_8 = 6.510e-6;%0.6479E-5;
CN_H_X_8 = 0.0;%0;
CN_g_X_8 = 0.6E-2;%0.60126E-2;

% (nu = 9)

CN_T_X9 = 18462.7308;%17431.6828;
CN_B_X_9 = 1.73326872;%1.7321008;
CN_D_X_9 = 6.530e-6;%0.6494E-5;
CN_H_X_9 = 0.0;%0;
CN_g_X_9 = 0.6E-2;%0.56139E-2;

% (nu = 10)

CN_T_X10 = 19235.8808;
CN_B_X_10 = 1.714045;
CN_D_X_10 = 0.64E-5;
CN_H_X_10 = 0;
CN_g_X_10 = 0.52322E-2;

CN_En_B1 = zeros(N_J,1);
CN_En_B2 = zeros(N_J,1);
CN_En_X1 = zeros(N_J,1);
CN_En_X2 = zeros(N_J,1);
CN_DEn_R1 = zeros(N_J,1);
CN_DEn_R2 = zeros(N_J,1);
CN_DEn_P1 = zeros(N_J,1);
CN_DEn_P2 = zeros(N_J,1);
CN_DEn_RQ21 = zeros(N_J,1);
CN_DEn_PQ12 = zeros(N_J,1);

CN_En_B1_1 = zeros(N_J,1);
CN_En_B2_1 = zeros(N_J,1);
CN_En_X1_1 = zeros(N_J,1);
CN_En_X2_1 = zeros(N_J,1);
CN_DEn_R1_1 = zeros(N_J,1);
CN_DEn_R2_1 = zeros(N_J,1);
CN_DEn_P1_1 = zeros(N_J,1);
CN_DEn_P2_1 = zeros(N_J,1);
CN_DEn_RQ21_1 = zeros(N_J,1);
CN_DEn_PQ12_1 = zeros(N_J,1);

CN_En_B1_2 = zeros(N_J,1);
CN_En_B2_2 = zeros(N_J,1);
CN_En_X1_2 = zeros(N_J,1);
CN_En_X2_2 = zeros(N_J,1);
CN_DEn_R1_2 = zeros(N_J,1);
CN_DEn_R2_2 = zeros(N_J,1);
CN_DEn_P1_2 = zeros(N_J,1);
CN_DEn_P2_2 = zeros(N_J,1);
CN_DEn_RQ21_2 = zeros(N_J,1);
CN_DEn_PQ12_2 = zeros(N_J,1);

CN_En_B1_3 = zeros(N_J,1);
CN_En_B2_3 = zeros(N_J,1);
CN_En_X1_3 = zeros(N_J,1);
CN_En_X2_3 = zeros(N_J,1);
CN_DEn_R1_3 = zeros(N_J,1);
CN_DEn_R2_3 = zeros(N_J,1);
CN_DEn_P1_3 = zeros(N_J,1);
CN_DEn_P2_3 = zeros(N_J,1);
CN_DEn_RQ21_3 = zeros(N_J,1);
CN_DEn_PQ12_3 = zeros(N_J,1);

CN_En_B1_4 = zeros(N_J,1);
CN_En_B2_4 = zeros(N_J,1);
CN_En_X1_4 = zeros(N_J,1);
CN_En_X2_4 = zeros(N_J,1);
CN_DEn_R1_4 = zeros(N_J,1);
CN_DEn_R2_4 = zeros(N_J,1);
CN_DEn_P1_4 = zeros(N_J,1);
CN_DEn_P2_4 = zeros(N_J,1);
CN_DEn_RQ21_4 = zeros(N_J,1);
CN_DEn_PQ12_4 = zeros(N_J,1);

CN_En_B1_5 = zeros(N_J,1);
CN_En_B2_5 = zeros(N_J,1);
CN_En_X1_5 = zeros(N_J,1);
CN_En_X2_5 = zeros(N_J,1);
CN_DEn_R1_5 = zeros(N_J,1);
CN_DEn_R2_5 = zeros(N_J,1);
CN_DEn_P1_5 = zeros(N_J,1);
CN_DEn_P2_5 = zeros(N_J,1);
CN_DEn_RQ21_5 = zeros(N_J,1);
CN_DEn_PQ12_5 = zeros(N_J,1);

CN_En_B1_6 = zeros(N_J,1);
CN_En_B2_6 = zeros(N_J,1);
CN_En_X1_6 = zeros(N_J,1);
CN_En_X2_6 = zeros(N_J,1);
CN_DEn_R1_6 = zeros(N_J,1);
CN_DEn_R2_6 = zeros(N_J,1);
CN_DEn_P1_6 = zeros(N_J,1);
CN_DEn_P2_6 = zeros(N_J,1);
CN_DEn_RQ21_6 = zeros(N_J,1);
CN_DEn_PQ12_6 = zeros(N_J,1);

CN_En_B1_7 = zeros(N_J,1);
CN_En_B2_7 = zeros(N_J,1);
CN_En_X1_7 = zeros(N_J,1);
CN_En_X2_7 = zeros(N_J,1);
CN_DEn_R1_7 = zeros(N_J,1);
CN_DEn_R2_7 = zeros(N_J,1);
CN_DEn_P1_7 = zeros(N_J,1);
CN_DEn_P2_7 = zeros(N_J,1);
CN_DEn_RQ21_7 = zeros(N_J,1);
CN_DEn_PQ12_7 = zeros(N_J,1);

CN_En_B1_8 = zeros(N_J,1);
CN_En_B2_8 = zeros(N_J,1);
CN_En_X1_8 = zeros(N_J,1);
CN_En_X2_8 = zeros(N_J,1);
CN_DEn_R1_8 = zeros(N_J,1);
CN_DEn_R2_8 = zeros(N_J,1);
CN_DEn_P1_8 = zeros(N_J,1);
CN_DEn_P2_8 = zeros(N_J,1);
CN_DEn_RQ21_8 = zeros(N_J,1);
CN_DEn_PQ12_8 = zeros(N_J,1);

CN_En_B1_9 = zeros(N_J,1);
CN_En_B2_9 = zeros(N_J,1);
CN_En_X1_9 = zeros(N_J,1);
CN_En_X2_9 = zeros(N_J,1);
CN_DEn_R1_9 = zeros(N_J,1);
CN_DEn_R2_9 = zeros(N_J,1);
CN_DEn_P1_9 = zeros(N_J,1);
CN_DEn_P2_9 = zeros(N_J,1);
CN_DEn_RQ21_9 = zeros(N_J,1);
CN_DEn_PQ12_9 = zeros(N_J,1);

CN_En_B1_10 = zeros(N_J,1);
CN_En_B2_10 = zeros(N_J,1);
CN_En_X1_10 = zeros(N_J,1);
CN_En_X2_10 = zeros(N_J,1);
CN_DEn_R1_10 = zeros(N_J,1);
CN_DEn_R2_10 = zeros(N_J,1);
CN_DEn_P1_10 = zeros(N_J,1);
CN_DEn_P2_10 = zeros(N_J,1);
CN_DEn_RQ21_10 = zeros(N_J,1);
CN_DEn_PQ12_10 = zeros(N_J,1);

% Second positive system, vibrational transition v'-v'' (3 - 6)
% Data Roux 1993

% Lower State B 3Pi


% (nu = 2)

N2_T_B_2 = 3381.501;
N2_A_B_2 = 42.1357;
N2_B_B_2 = 1.59227;
N2_D_B_2 = 0.596E-5;
N2_H_B_2 = 0.65E-11;
N2_Y_B_2 = N2_A_B_2/N2_B_B_2;
y1_B_2 = N2_Y_B_2*(N2_Y_B_2-4)+4/3;
y2_B_2 = N2_Y_B_2*(N2_Y_B_2-1)-4/9;
N2_l_B_2 = -0.2104;
N2_g_B_2 = 0.17E-2;
N2_o_B_2 = 1.1471;
N2_p_B_2 = 0.445E-2;
N2_q_B_2 = 0.834E-4;

% (nu = 3)

N2_T_B_3 = 5028.896;
N2_A_B_3 = 42.023;
N2_B_B_3 = 1.57380;
N2_D_B_3 = 0.598E-5;
N2_Y_B_3 = N2_A_B_3/N2_B_B_3;
y1_B_3 = N2_Y_B_3*(N2_Y_B_3-4)+4/3;
y2_B_3 = N2_Y_B_3*(N2_Y_B_3-1)-4/9;
N2_l_B_3 = -0.2101;
N2_g_B_3 = 0.18E-2;
N2_o_B_3 = 1.1447;
N2_p_B_3 = 0.442E-2;
N2_q_B_3 = 0.819E-4;

% (nu = 4)

N2_T_B_4 = 6647.314;
N2_A_B_4 = 42.023;
N2_B_B_4 = 1.55520;
N2_D_B_4 = 0.598E-5;
N2_Y_B_4 = N2_A_B_4/N2_B_B_4;
y1_B_4 = N2_Y_B_4*(N2_Y_B_4-4)+4/3;
y2_B_4 = N2_Y_B_4*(N2_Y_B_4-1)-4/9;
N2_l_B_4 = -0.2154;
N2_g_B_4 = 0.15E-2;
N2_o_B_4 = 1.1422;
N2_p_B_4 = 0.44E-2;
N2_q_B_4 = 0.789E-4;

% (nu = 5)

N2_T_B_5 = 8236.686;
N2_A_B_5 = 41.955;
N2_B_B_5 = 1.53644;
N2_D_B_5 = 5.83E-06;
N2_Y_B_5 = N2_A_B_5/N2_B_B_5;
y1_B_5 = N2_Y_B_5*(N2_Y_B_5-4)+4/3;
y2_B_5 = N2_Y_B_5*(N2_Y_B_5-1)-4/9;
N2_l_B_5 = -0.2187;
N2_g_B_5 = 0.19E-2;
N2_o_B_5 = 1.1397;
N2_p_B_5 = 0.439E-2;
N2_q_B_5 = 0.768E-4;

% (nu = 6)

N2_T_B = 9796.934;
N2_A_B = 41.883;
N2_B_B = 1.51756;
N2_D_B = 5.89E-06;
N2_Y_B = N2_A_B/N2_B_B;
y1_B = N2_Y_B*(N2_Y_B-4)+4/3;
y2_B = N2_Y_B*(N2_Y_B-1)-4/9;
N2_l_B = -0.2188;
N2_g_B = 0.19E-2;
N2_o_B = 1.1372;
N2_p_B = 0.437E-2;
N2_q_B = 0.762E-4;

% (nu = 7)

N2_T_B_7 = 11327.967;
N2_A_B_7 = 41.809;
N2_B_B_7 = 1.49863;
N2_D_B_7 = 0.593E-5;
N2_Y_B_7 = N2_A_B_7/N2_B_B_7;
y1_B_7 = N2_Y_B_7*(N2_Y_B_7-4)+4/3;
y2_B_7 = N2_Y_B_7*(N2_Y_B_7-1)-4/9;
N2_l_B_7 = -0.2231;
N2_g_B_7 = 0.15E-2;
N2_o_B_7 = 1.1338;
N2_p_B_7 = 0.433E-2;
N2_q_B_7 = 0.724E-4;

% upper state C 3Pi

% (nu = 0)

N2_T_C_0 = 29670.942;
N2_A_C_0 = 39.134;
N2_B_C_0 = 1.81530;
N2_D_C_0 = 0.595E-5;
N2_Y_C_0 = N2_A_C_0/N2_B_C_0;
y1_C_0 = N2_Y_C_0*(N2_Y_C_0-4)+4/3;
y2_C_0 = N2_Y_C_0*(N2_Y_C_0-1)-4/9;
N2_l_C_0 = 0.23;
N2_g_C_0 = 0.166E-2;
N2_o_C_0 = 1.106;
N2_p_C_0 = 0.503E-2;
N2_q_C_0 = -1.52E-4;

% (nu = 1)

N2_T_C_1 = 31665.392;
N2_A_C_1 = 38.544;
N2_B_C_1 = 1.79364;
N2_D_C_1 = 6.14E-06;
N2_Y_C_1 = N2_A_C_1/N2_B_C_1;
y1_C_1 = N2_Y_C_1*(N2_Y_C_1-4)+4/3;
y2_C_1 = N2_Y_C_1*(N2_Y_C_1-1)-4/9;
N2_l_C_1 = 0.234;
N2_g_C_1 = 0.159E-2;
N2_o_C_1 = 1.099;
N2_p_C_1 = 0.527E-2;
N2_q_C_1 = -0.36E-4;

% (nu = 2)

N2_T_C_2 = 33606.182;
N2_A_C_2 = 37.757;
N2_B_C_2 = 1.76934;
N2_D_C_2 = 0.633E-5;
N2_Y_C_2 = N2_A_C_2/N2_B_C_2;
y1_C_2 = N2_Y_C_2*(N2_Y_C_2-4)+4/3;
y2_C_2 = N2_Y_C_2*(N2_Y_C_2-1)-4/9;
N2_l_C_2 = 0.247;
N2_g_C_2 = 0.299E-2;
N2_o_C_2 = 1.081;
N2_p_C_2 = 0.496E-2;
N2_q_C_2 = -1.2E-4;

% (nu = 3)

N2_T_C = 35480.326;
N2_A_C = 36.634;
N2_B_C = 1.74028;
N2_D_C = 6.9500E-06;
N2_Y_C = N2_A_C/N2_B_C;
y1_C = N2_Y_C*(N2_Y_C-4)+4/3;
y2_C = N2_Y_C*(N2_Y_C-1)-4/9;
N2_l_C = 0.265;
N2_g_C = 0.495E-2;
N2_o_C = 1.054;
N2_p_C = 0.55E-2;
N2_q_C = -0.82E-4;

% (nu = 4)

N2_T_C_4 = 37261.611;
N2_A_C_4 = 34.826;
N2_B_C_4 = 1.70097;
N2_D_C_4 = 1.051E-5;
N2_Y_C_4 = N2_A_C_4/N2_B_C_4;
y1_C_4 = N2_Y_C_4*(N2_Y_C_4-4)+4/3;
y2_C_4 = N2_Y_C_4*(N2_Y_C_4-1)-4/9;
N2_l_C_4 = 0.290;
N2_g_C_4 = 0.706E-2;
N2_o_C_4 = 0.971;
N2_p_C_4 = 0.544E-2;
N2_q_C_4 = -0.97E-4;

%(3-6)
N2_En_B1_6e = zeros(N_J,1);
N2_En_B2_6e = zeros(N_J,1);
N2_En_B3_6e = zeros(N_J,1);
N2_En_B1_6f = zeros(N_J,1);
N2_En_B2_6f = zeros(N_J,1);
N2_En_B3_6f = zeros(N_J,1);
N2_En_C1_3e = zeros(N_J,1);
N2_En_C2_3e = zeros(N_J,1);
N2_En_C3_3e = zeros(N_J,1);
N2_En_C1_3f = zeros(N_J,1);
N2_En_C2_3f = zeros(N_J,1);
N2_En_C3_3f = zeros(N_J,1);
N2_DEn_R1_36f = zeros(N_J,1);
N2_DEn_R2_36f = zeros(N_J,1);
N2_DEn_R3_36f = zeros(N_J,1);
N2_DEn_P1_36f = zeros(N_J,1);
N2_DEn_P2_36f = zeros(N_J,1);
N2_DEn_P3_36f = zeros(N_J,1);
N2_DEn_Q1_36f = zeros(N_J,1);
N2_DEn_Q2_36f = zeros(N_J,1);
N2_DEn_Q3_36f = zeros(N_J,1);
N2_DEn_R1_36e = zeros(N_J,1);
N2_DEn_R2_36e = zeros(N_J,1);
N2_DEn_R3_36e = zeros(N_J,1);
N2_DEn_P1_36e = zeros(N_J,1);
N2_DEn_P2_36e = zeros(N_J,1);
N2_DEn_P3_36e = zeros(N_J,1);
N2_DEn_Q1_36e = zeros(N_J,1);
N2_DEn_Q2_36e = zeros(N_J,1);
N2_DEn_Q3_36e = zeros(N_J,1);
% (2-5)
N2_En_B1_5e = zeros(N_J,1);
N2_En_B2_5e = zeros(N_J,1);
N2_En_B3_5e = zeros(N_J,1);
N2_En_C1_2e = zeros(N_J,1);
N2_En_C2_2e = zeros(N_J,1);
N2_En_C3_2e = zeros(N_J,1);
N2_En_B1_5f = zeros(N_J,1);
N2_En_B2_5f = zeros(N_J,1);
N2_En_B3_5f = zeros(N_J,1);
N2_En_C1_2f = zeros(N_J,1);
N2_En_C2_2f = zeros(N_J,1);
N2_En_C3_2f = zeros(N_J,1);
N2_DEn_R1_25e = zeros(N_J,1);
N2_DEn_R2_25e = zeros(N_J,1);
N2_DEn_R3_25e = zeros(N_J,1);
N2_DEn_P1_25e = zeros(N_J,1);
N2_DEn_P2_25e = zeros(N_J,1);
N2_DEn_P3_25e = zeros(N_J,1);
N2_DEn_Q1_25e = zeros(N_J,1);
N2_DEn_Q2_25e = zeros(N_J,1);
N2_DEn_Q3_25e = zeros(N_J,1);
N2_DEn_R1_25f = zeros(N_J,1);
N2_DEn_R2_25f = zeros(N_J,1);
N2_DEn_R3_25f = zeros(N_J,1);
N2_DEn_P1_25f = zeros(N_J,1);
N2_DEn_P2_25f = zeros(N_J,1);
N2_DEn_P3_25f = zeros(N_J,1);
N2_DEn_Q1_25f = zeros(N_J,1);
N2_DEn_Q2_25f = zeros(N_J,1);
N2_DEn_Q3_25f = zeros(N_J,1);
% (4-7)
N2_En_B1_7e = zeros(N_J,1);
N2_En_B2_7e = zeros(N_J,1);
N2_En_B3_7e = zeros(N_J,1);
N2_En_C1_4e = zeros(N_J,1);
N2_En_C2_4e = zeros(N_J,1);
N2_En_C3_4e = zeros(N_J,1);
N2_En_B1_7f = zeros(N_J,1);
N2_En_B2_7f = zeros(N_J,1);
N2_En_B3_7f = zeros(N_J,1);
N2_En_C1_4f = zeros(N_J,1);
N2_En_C2_4f = zeros(N_J,1);
N2_En_C3_4f = zeros(N_J,1);
N2_DEn_R1_47e = zeros(N_J,1);
N2_DEn_R2_47e = zeros(N_J,1);
N2_DEn_R3_47e = zeros(N_J,1);
N2_DEn_P1_47e = zeros(N_J,1);
N2_DEn_P2_47e = zeros(N_J,1);
N2_DEn_P3_47e = zeros(N_J,1);
N2_DEn_Q1_47e = zeros(N_J,1);
N2_DEn_Q2_47e = zeros(N_J,1);
N2_DEn_Q3_47e = zeros(N_J,1);
N2_DEn_R1_47f = zeros(N_J,1);
N2_DEn_R2_47f = zeros(N_J,1);
N2_DEn_R3_47f = zeros(N_J,1);
N2_DEn_P1_47f = zeros(N_J,1);
N2_DEn_P2_47f = zeros(N_J,1);
N2_DEn_P3_47f = zeros(N_J,1);
N2_DEn_Q1_47f = zeros(N_J,1);
N2_DEn_Q2_47f = zeros(N_J,1);
N2_DEn_Q3_47f = zeros(N_J,1);
% (1-4)
N2_En_B1_4e = zeros(N_J,1);
N2_En_B2_4e = zeros(N_J,1);
N2_En_B3_4e = zeros(N_J,1);
N2_En_C1_1e = zeros(N_J,1);
N2_En_C2_1e = zeros(N_J,1);
N2_En_C3_1e = zeros(N_J,1);
N2_En_B1_4f = zeros(N_J,1);
N2_En_B2_4f = zeros(N_J,1);
N2_En_B3_4f = zeros(N_J,1);
N2_En_C1_1f = zeros(N_J,1);
N2_En_C2_1f = zeros(N_J,1);
N2_En_C3_1f = zeros(N_J,1);
N2_DEn_R1_14e = zeros(N_J,1);
N2_DEn_R2_14e = zeros(N_J,1);
N2_DEn_R3_14e = zeros(N_J,1);
N2_DEn_P1_14e = zeros(N_J,1);
N2_DEn_P2_14e = zeros(N_J,1);
N2_DEn_P3_14e = zeros(N_J,1);
N2_DEn_Q1_14e = zeros(N_J,1);
N2_DEn_Q2_14e = zeros(N_J,1);
N2_DEn_Q3_14e = zeros(N_J,1);
N2_DEn_R1_14f = zeros(N_J,1);
N2_DEn_R2_14f = zeros(N_J,1);
N2_DEn_R3_14f = zeros(N_J,1);
N2_DEn_P1_14f = zeros(N_J,1);
N2_DEn_P2_14f = zeros(N_J,1);
N2_DEn_P3_14f = zeros(N_J,1);
N2_DEn_Q1_14f = zeros(N_J,1);
N2_DEn_Q2_14f = zeros(N_J,1);
N2_DEn_Q3_14f = zeros(N_J,1);
% (0-3)
N2_En_B1_3e = zeros(N_J,1);
N2_En_B2_3e = zeros(N_J,1);
N2_En_B3_3e = zeros(N_J,1);
N2_En_C1_0e = zeros(N_J,1);
N2_En_C2_0e = zeros(N_J,1);
N2_En_C3_0e = zeros(N_J,1);
N2_En_B1_3f = zeros(N_J,1);
N2_En_B2_3f = zeros(N_J,1);
N2_En_B3_3f = zeros(N_J,1);
N2_En_C1_0f = zeros(N_J,1);
N2_En_C2_0f = zeros(N_J,1);
N2_En_C3_0f = zeros(N_J,1);
N2_DEn_R1_03e = zeros(N_J,1);
N2_DEn_R2_03e = zeros(N_J,1);
N2_DEn_R3_03e = zeros(N_J,1);
N2_DEn_P1_03e = zeros(N_J,1);
N2_DEn_P2_03e = zeros(N_J,1);
N2_DEn_P3_03e = zeros(N_J,1);
N2_DEn_Q1_03e = zeros(N_J,1);
N2_DEn_Q2_03e = zeros(N_J,1);
N2_DEn_Q3_03e = zeros(N_J,1);
N2_DEn_R1_03f = zeros(N_J,1);
N2_DEn_R2_03f = zeros(N_J,1);
N2_DEn_R3_03f = zeros(N_J,1);
N2_DEn_P1_03f = zeros(N_J,1);
N2_DEn_P2_03f = zeros(N_J,1);
N2_DEn_P3_03f = zeros(N_J,1);
N2_DEn_Q1_03f = zeros(N_J,1);
N2_DEn_Q2_03f = zeros(N_J,1);
N2_DEn_Q3_03f = zeros(N_J,1);
% (0-2)
N2_En_B1_2e = zeros(N_J,1);
N2_En_B2_2e = zeros(N_J,1);
N2_En_B3_2e = zeros(N_J,1);
N2_En_B1_2f = zeros(N_J,1);
N2_En_B2_2f = zeros(N_J,1);
N2_En_B3_2f = zeros(N_J,1);
N2_DEn_R1_02e = zeros(N_J,1);
N2_DEn_R2_02e = zeros(N_J,1);
N2_DEn_R3_02e = zeros(N_J,1);
N2_DEn_P1_02e = zeros(N_J,1);
N2_DEn_P2_02e = zeros(N_J,1);
N2_DEn_P3_02e = zeros(N_J,1);
N2_DEn_Q1_02e = zeros(N_J,1);
N2_DEn_Q2_02e = zeros(N_J,1);
N2_DEn_Q3_02e = zeros(N_J,1);
N2_DEn_R1_02f = zeros(N_J,1);
N2_DEn_R2_02f = zeros(N_J,1);
N2_DEn_R3_02f = zeros(N_J,1);
N2_DEn_P1_02f = zeros(N_J,1);
N2_DEn_P2_02f = zeros(N_J,1);
N2_DEn_P3_02f = zeros(N_J,1);
N2_DEn_Q1_02f = zeros(N_J,1);
N2_DEn_Q2_02f = zeros(N_J,1);
N2_DEn_Q3_02f = zeros(N_J,1);
% (1-3)
N2_DEn_R1_13e = zeros(N_J,1);
N2_DEn_R2_13e = zeros(N_J,1);
N2_DEn_R3_13e = zeros(N_J,1);
N2_DEn_P1_13e = zeros(N_J,1);
N2_DEn_P2_13e = zeros(N_J,1);
N2_DEn_P3_13e = zeros(N_J,1);
N2_DEn_Q1_13e = zeros(N_J,1);
N2_DEn_Q2_13e = zeros(N_J,1);
N2_DEn_Q3_13e = zeros(N_J,1);
N2_DEn_R1_13f = zeros(N_J,1);
N2_DEn_R2_13f = zeros(N_J,1);
N2_DEn_R3_13f = zeros(N_J,1);
N2_DEn_P1_13f = zeros(N_J,1);
N2_DEn_P2_13f = zeros(N_J,1);
N2_DEn_P3_13f = zeros(N_J,1);
N2_DEn_Q1_13f = zeros(N_J,1);
N2_DEn_Q2_13f = zeros(N_J,1);
N2_DEn_Q3_13f = zeros(N_J,1);
% (2-4)
N2_DEn_R1_24e = zeros(N_J,1);
N2_DEn_R2_24e = zeros(N_J,1);
N2_DEn_R3_24e = zeros(N_J,1);
N2_DEn_P1_24e = zeros(N_J,1);
N2_DEn_P2_24e = zeros(N_J,1);
N2_DEn_P3_24e = zeros(N_J,1);
N2_DEn_Q1_24e = zeros(N_J,1);
N2_DEn_Q2_24e = zeros(N_J,1);
N2_DEn_Q3_24e = zeros(N_J,1);
N2_DEn_R1_24f = zeros(N_J,1);
N2_DEn_R2_24f = zeros(N_J,1);
N2_DEn_R3_24f = zeros(N_J,1);
N2_DEn_P1_24f = zeros(N_J,1);
N2_DEn_P2_24f = zeros(N_J,1);
N2_DEn_P3_24f = zeros(N_J,1);
N2_DEn_Q1_24f = zeros(N_J,1);
N2_DEn_Q2_24f = zeros(N_J,1);
N2_DEn_Q3_24f = zeros(N_J,1);

    for J=1:N_J
        
        % N2(+) energies 
        En_B1(J) = T + B_B*J*(J+1) - D_B*J^2*(J+1)^2+H_B*J^3*(J+1)^3 +...
            0.5*g_B*J; 
        En_B1_1(J) = T_1 + B_B_1*J*(J+1) - D_B_1*J^2*(J+1)^2+H_B_1*J^3*(J+1)^3 +...
            0.5*g_B_1*J; 
        En_B2(J) = T + B_B*J*(J+1) - D_B*J^2*(J+1)^2+H_B*J^3*(J+1)^3 -...
            0.5*g_B*(J+1); 
        En_B2_1(J) = T_1 + B_B_1*J*(J+1) - D_B_1*J^2*(J+1)^2+H_B_1*J^3*(J+1)^3 -...
            0.5*g_B_1*(J+1);
        En_X1(J) =     B_X*J*(J+1) - D_X*J^2*(J+1)^2 + 0.5*g_X*J;
        En_X1_1(J) = T_X1 + B_X_1*J*(J+1) - D_X_1*J^2*(J+1)^2 + 0.5*g_X_1*J;
        En_X2(J) =    B_X*J*(J+1) - D_X*J^2*(J+1)^2 - 0.5*g_X*(J+1);
        En_X2_1(J) = T_X1 + B_X_1*J*(J+1) - D_X_1*J^2*(J+1)^2 - 0.5*g_X_1*(J+1);
        
        % CN energies
        CN_En_B1(J) = CN_T + CN_B_B*J*(J+1) - CN_D_B*J^2*(J+1)^2 + ...
             CN_H_B*J^3*(J+1)^3 + CN_L_B*J^4*(J+1)^4 + CN_M_B*J^5*(J+1)^5 ...
             + CN_N_B*J^6*(J+1)^6 + 0.5*CN_g_B*J; 
        CN_En_B1_1(J) = CN_T_1 + CN_B_B_1*J*(J+1) - CN_D_B_1*J^2*(J+1)^2 + ...
             CN_H_B_1*J^3*(J+1)^3 + CN_L_B_1*J^4*(J+1)^4 + CN_M_B_1*J^5*(J+1)^5 ...
             + CN_N_B_1*J^6*(J+1)^6 + 0.5*CN_g_B_1*J;
        CN_En_B1_2(J) = CN_T_2 + CN_B_B_2*J*(J+1) - CN_D_B_2*J^2*(J+1)^2 +...
             CN_H_B_2*J^3*(J+1)^3 + CN_L_B_2*J^4*(J+1)^4 + CN_M_B_2*J^5*(J+1)^5 ...
             + CN_N_B_2*J^6*(J+1)^6 + 0.5*CN_g_B_2*J;
        CN_En_B1_3(J) = CN_T_3 + CN_B_B_3*J*(J+1) - CN_D_B_3*J^2*(J+1)^2 +...
             CN_H_B_3*J^3*(J+1)^3 + CN_L_B_3*J^4*(J+1)^4 + CN_M_B_3*J^5*(J+1)^5 ...
             + CN_N_B_3*J^6*(J+1)^6 + 0.5*CN_g_B_3*J;
        CN_En_B1_4(J) = CN_T_4 + CN_B_B_4*J*(J+1) - CN_D_B_4*J^2*(J+1)^2 +...
             CN_H_B_4*J^3*(J+1)^3 + CN_L_B_4*J^4*(J+1)^4 + CN_M_B_4*J^5*(J+1)^5 ...
             + CN_N_B_4*J^6*(J+1)^6 + 0.5*CN_g_B_4*J;
        CN_En_B1_5(J) = CN_T_5 + CN_B_B_5*J*(J+1) - CN_D_B_5*J^2*(J+1)^2 +...
             CN_H_B_5*J^3*(J+1)^3 + CN_L_B_5*J^4*(J+1)^4 + CN_M_B_5*J^5*(J+1)^5 ...
             + CN_N_B_5*J^6*(J+1)^6 + 0.5*CN_g_B_5*J;
        CN_En_B1_6(J) = CN_T_6 + CN_B_B_6*J*(J+1) - CN_D_B_6*J^2*(J+1)^2 +...
             CN_H_B_6*J^3*(J+1)^3 + CN_L_B_6*J^4*(J+1)^4 + CN_M_B_6*J^5*(J+1)^5 ...
             + CN_N_B_6*J^6*(J+1)^6 + 0.5*CN_g_B_6*J;     
        CN_En_B1_7(J) = CN_T_7 + CN_B_B_7*J*(J+1) - CN_D_B_7*J^2*(J+1)^2 +...
             CN_H_B_7*J^3*(J+1)^3 + CN_L_B_7*J^4*(J+1)^4 + CN_M_B_7*J^5*(J+1)^5 ...
             + CN_N_B_7*J^6*(J+1)^6 + 0.5*CN_g_B_7*J;
        CN_En_B1_8(J) = CN_T_8 + CN_B_B_8*J*(J+1) - CN_D_B_8*J^2*(J+1)^2 +...
             CN_H_B_8*J^3*(J+1)^3 + CN_L_B_8*J^4*(J+1)^4 + CN_M_B_8*J^5*(J+1)^5 ...
             + CN_N_B_8*J^6*(J+1)^6 + 0.5*CN_g_B_8*J;
        CN_En_B1_9(J) = CN_T_9 + CN_B_B_9*J*(J+1) - CN_D_B_9*J^2*(J+1)^2 +...
             CN_H_B_9*J^3*(J+1)^3 + CN_L_B_9*J^4*(J+1)^4 + CN_M_B_9*J^5*(J+1)^5 ...
             + CN_N_B_9*J^6*(J+1)^6 + 0.5*CN_g_B_9*J;
        CN_En_B1_10(J) = CN_T_10 + CN_B_B_10*J*(J+1) - CN_D_B_10*J^2*(J+1)^2 +...
             CN_H_B_10*J^3*(J+1)^3 + CN_L_B_10*J^4*(J+1)^4 + CN_M_B_10*J^5*(J+1)^5 ...
             + CN_N_B_10*J^6*(J+1)^6 + 0.5*CN_g_B_10*J;
        CN_En_B2(J) = CN_T + CN_B_B*J*(J+1) - CN_D_B*J^2*(J+1)^2 + ...
             + CN_H_B*J^3*(J+1)^3 + CN_L_B*J^4*(J+1)^4 + CN_M_B*J^5*(J+1)^5 ...
             + CN_N_B*J^6*(J+1)^6 - 0.5*CN_g_B*(J+1); 
        CN_En_B2_1(J) = CN_T_1 + CN_B_B_1*J*(J+1) - CN_D_B_1*J^2*(J+1)^2 + ...
             + CN_H_B_1*J^3*(J+1)^3 + CN_L_B_1*J^4*(J+1)^4 + CN_M_B_1*J^5*(J+1)^5 ...
             + CN_N_B_1*J^6*(J+1)^6 - 0.5*CN_g_B_1*(J+1); 
        CN_En_B2_2(J) = CN_T_2 + CN_B_B_2*J*(J+1) - CN_D_B_2*J^2*(J+1)^2 + ...
             + CN_H_B_2*J^3*(J+1)^3 + CN_L_B_2*J^4*(J+1)^4 + CN_M_B_2*J^5*(J+1)^5 ...
             + CN_N_B_2*J^6*(J+1)^6 - 0.5*CN_g_B_2*(J+1);    
        CN_En_B2_3(J) = CN_T_3 + CN_B_B_3*J*(J+1) - CN_D_B_3*J^2*(J+1)^2 + ...
             + CN_H_B_3*J^3*(J+1)^3 + CN_L_B_3*J^4*(J+1)^4 + CN_M_B_3*J^5*(J+1)^5 ...
             + CN_N_B_3*J^6*(J+1)^6 - 0.5*CN_g_B_3*(J+1);
        CN_En_B2_4(J) = CN_T_4 + CN_B_B_4*J*(J+1) - CN_D_B_4*J^2*(J+1)^2 + ...
             + CN_H_B_4*J^3*(J+1)^3 + CN_L_B_4*J^4*(J+1)^4 + CN_M_B_4*J^5*(J+1)^5 ...
             + CN_N_B_4*J^6*(J+1)^6 - 0.5*CN_g_B_4*(J+1);
        CN_En_B2_5(J) = CN_T_5 + CN_B_B_5*J*(J+1) - CN_D_B_5*J^2*(J+1)^2 + ...
             + CN_H_B_5*J^3*(J+1)^3 + CN_L_B_5*J^4*(J+1)^4 + CN_M_B_5*J^5*(J+1)^5 ...
             + CN_N_B_5*J^6*(J+1)^6 - 0.5*CN_g_B_5*(J+1); 
        CN_En_B2_6(J) = CN_T_6 + CN_B_B_6*J*(J+1) - CN_D_B_6*J^2*(J+1)^2 + ...
             + CN_H_B_6*J^3*(J+1)^3 + CN_L_B_6*J^4*(J+1)^4 + CN_M_B_6*J^5*(J+1)^5 ...
             + CN_N_B_6*J^6*(J+1)^6 - 0.5*CN_g_B_6*(J+1);  
        CN_En_B2_7(J) = CN_T_7 + CN_B_B_7*J*(J+1) - CN_D_B_7*J^2*(J+1)^2 + ...
             + CN_H_B_7*J^3*(J+1)^3 + CN_L_B_7*J^4*(J+1)^4 + CN_M_B_7*J^5*(J+1)^5 ...
             + CN_N_B_7*J^6*(J+1)^6 - 0.5*CN_g_B_7*(J+1);
        CN_En_B2_8(J) = CN_T_8 + CN_B_B_8*J*(J+1) - CN_D_B_8*J^2*(J+1)^2 + ...
             + CN_H_B_8*J^3*(J+1)^3 + CN_L_B_8*J^4*(J+1)^4 + CN_M_B_8*J^5*(J+1)^5 ...
             + CN_N_B_8*J^6*(J+1)^6 - 0.5*CN_g_B_8*(J+1);
        CN_En_B2_9(J) = CN_T_9 + CN_B_B_9*J*(J+1) - CN_D_B_9*J^2*(J+1)^2 + ...
             + CN_H_B_9*J^3*(J+1)^3 + CN_L_B_9*J^4*(J+1)^4 + CN_M_B_9*J^5*(J+1)^5 ...
             + CN_N_B_9*J^6*(J+1)^6 - 0.5*CN_g_B_9*(J+1);
        CN_En_B2_10(J) = CN_T_10 + CN_B_B_10*J*(J+1) - CN_D_B_10*J^2*(J+1)^2 + ...
             + CN_H_B_10*J^3*(J+1)^3 + CN_L_B_10*J^4*(J+1)^4 + CN_M_B_10*J^5*(J+1)^5 ...
             + CN_N_B_10*J^6*(J+1)^6 - 0.5*CN_g_B_10*(J+1);
        CN_En_X1(J) = CN_T_X + CN_B_X*J*(J+1) - CN_D_X*J^2*(J+1)^2 + ... 
            CN_H_X*J^3*(J+1)^3 + 0.5*CN_g_X*J;
        CN_En_X1_1(J) = CN_T_X1 + CN_B_X_1*J*(J+1) - CN_D_X_1*J^2*(J+1)^2 ...
            + CN_H_X_1*J^3*(J+1)^3 + 0.5*CN_g_X_1*J;
        CN_En_X1_2(J) = CN_T_X2 + CN_B_X_2*J*(J+1) - CN_D_X_2*J^2*(J+1)^2 ...
            + CN_H_X_2*J^3*(J+1)^3 + 0.5*CN_g_X_2*J;
        CN_En_X1_3(J) = CN_T_X3 + CN_B_X_3*J*(J+1) - CN_D_X_3*J^2*(J+1)^2 ...
            + CN_H_X_3*J^3*(J+1)^3 + 0.5*CN_g_X_3*J;  
        CN_En_X1_4(J) = CN_T_X4 + CN_B_X_4*J*(J+1) - CN_D_X_4*J^2*(J+1)^2 ...
            + CN_H_X_4*J^3*(J+1)^3 + 0.5*CN_g_X_4*J; 
        CN_En_X1_5(J) = CN_T_X5 + CN_B_X_5*J*(J+1) - CN_D_X_5*J^2*(J+1)^2 ...
            + CN_H_X_5*J^3*(J+1)^3 + 0.5*CN_g_X_5*J;   
        CN_En_X1_6(J) = CN_T_X6 + CN_B_X_6*J*(J+1) - CN_D_X_6*J^2*(J+1)^2 ...
            + CN_H_X_6*J^3*(J+1)^3 + 0.5*CN_g_X_6*J;   
        CN_En_X1_7(J) = CN_T_X7 + CN_B_X_7*J*(J+1) - CN_D_X_7*J^2*(J+1)^2 ...
            + CN_H_X_7*J^3*(J+1)^3 + 0.5*CN_g_X_7*J;
        CN_En_X1_8(J) = CN_T_X8 + CN_B_X_8*J*(J+1) - CN_D_X_8*J^2*(J+1)^2 ...
            + CN_H_X_8*J^3*(J+1)^3 + 0.5*CN_g_X_8*J;
        CN_En_X1_9(J) = CN_T_X9 + CN_B_X_9*J*(J+1) - CN_D_X_9*J^2*(J+1)^2 ...
            + CN_H_X_9*J^3*(J+1)^3 + 0.5*CN_g_X_9*J;
        CN_En_X1_10(J) = CN_T_X10 + CN_B_X_10*J*(J+1) - CN_D_X_10*J^2*(J+1)^2 ...
            + CN_H_X_10*J^3*(J+1)^3 + 0.5*CN_g_X_10*J;
        CN_En_X2(J) =  CN_T_X + CN_B_X*J*(J+1) - CN_D_X*J^2*(J+1)^2 ...
            + CN_H_X*J^3*(J+1)^3 - 0.5*CN_g_X*(J+1);
        CN_En_X2_1(J) = CN_T_X1 + CN_B_X_1*J*(J+1) - CN_D_X_1*J^2*(J+1)^2 ...
           + CN_H_X_1*J^3*(J+1)^3 - 0.5*CN_g_X_1*(J+1);
        CN_En_X2_2(J) = CN_T_X2 + CN_B_X_2*J*(J+1) - CN_D_X_2*J^2*(J+1)^2 ...
           + CN_H_X_2*J^3*(J+1)^3 - 0.5*CN_g_X_2*(J+1);  
        CN_En_X2_3(J) = CN_T_X3 + CN_B_X_3*J*(J+1) - CN_D_X_3*J^2*(J+1)^2 ...
           + CN_H_X_3*J^3*(J+1)^3 - 0.5*CN_g_X_3*(J+1); 
        CN_En_X2_4(J) = CN_T_X4 + CN_B_X_4*J*(J+1) - CN_D_X_4*J^2*(J+1)^2 ...
           + CN_H_X_4*J^3*(J+1)^3 - 0.5*CN_g_X_4*(J+1); 
        CN_En_X2_5(J) = CN_T_X5 + CN_B_X_5*J*(J+1) - CN_D_X_5*J^2*(J+1)^2 ...
           + CN_H_X_5*J^3*(J+1)^3 - 0.5*CN_g_X_5*(J+1);    
        CN_En_X2_6(J) = CN_T_X6 + CN_B_X_6*J*(J+1) - CN_D_X_6*J^2*(J+1)^2 ...
           + CN_H_X_6*J^3*(J+1)^3 - 0.5*CN_g_X_6*(J+1);  
        CN_En_X2_7(J) = CN_T_X7 + CN_B_X_7*J*(J+1) - CN_D_X_7*J^2*(J+1)^2 ...
           + CN_H_X_7*J^3*(J+1)^3 - 0.5*CN_g_X_7*(J+1);
        CN_En_X2_8(J) = CN_T_X8 + CN_B_X_8*J*(J+1) - CN_D_X_8*J^2*(J+1)^2 ...
           + CN_H_X_8*J^3*(J+1)^3 - 0.5*CN_g_X_8*(J+1);
        CN_En_X2_9(J) = CN_T_X9 + CN_B_X_9*J*(J+1) - CN_D_X_9*J^2*(J+1)^2 ...
           + CN_H_X_9*J^3*(J+1)^3 - 0.5*CN_g_X_9*(J+1);
        CN_En_X2_10(J) = CN_T_X10 + CN_B_X_10*J*(J+1) - CN_D_X_10*J^2*(J+1)^2 ...
           + CN_H_X_10*J^3*(J+1)^3 - 0.5*CN_g_X_10*(J+1); 
       
       % N2 energies (solution of the secular equation)

        [N2_En_B1_7e(J), N2_En_B2_7e(J), N2_En_B3_7e(J)] = N2_Energy_e(N2_T_B_7,...
            N2_A_B_7,N2_B_B_7,N2_D_B_7,N2_l_B_7,N2_g_B_7,N2_o_B_7,N2_p_B_7,N2_q_B_7,J);
        [N2_En_B1_7f(J), N2_En_B2_7f(J), N2_En_B3_7f(J)] = N2_Energy_f(N2_T_B_7,...
            N2_A_B_7,N2_B_B_7,N2_D_B_7,N2_l_B_7,N2_g_B_7,N2_o_B_7,N2_p_B_7,N2_q_B_7,J);             
        [N2_En_B1_6e(J), N2_En_B2_6e(J), N2_En_B3_6e(J)] = N2_Energy_e(N2_T_B,...
            N2_A_B,N2_B_B,N2_D_B,N2_l_B,N2_g_B,N2_o_B,N2_p_B,N2_q_B,J);
        [N2_En_B1_6f(J), N2_En_B2_6f(J), N2_En_B3_6f(J)] = N2_Energy_f(N2_T_B,...
            N2_A_B,N2_B_B,N2_D_B,N2_l_B,N2_g_B,N2_o_B,N2_p_B,N2_q_B,J);
        [N2_En_B1_5e(J), N2_En_B2_5e(J), N2_En_B3_5e(J)] = N2_Energy_e(N2_T_B_5,...
            N2_A_B_5,N2_B_B_5,N2_D_B_5,N2_l_B_5,N2_g_B_5,N2_o_B_5,N2_p_B_5,N2_q_B_5,J);
        [N2_En_B1_5f(J), N2_En_B2_5f(J), N2_En_B3_5f(J)] = N2_Energy_f(N2_T_B_5,...
            N2_A_B_5,N2_B_B_5,N2_D_B_5,N2_l_B_5,N2_g_B_5,N2_o_B_5,N2_p_B_5,N2_q_B_5,J);       
        [N2_En_B1_4e(J), N2_En_B2_4e(J), N2_En_B3_4e(J)] = N2_Energy_e(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,J);
        [N2_En_B1_4f(J), N2_En_B2_4f(J), N2_En_B3_4f(J)] = N2_Energy_f(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,J); 
        [N2_En_B1_3e(J), N2_En_B2_3e(J), N2_En_B3_3e(J)] = N2_Energy_e(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,J);
        [N2_En_B1_3f(J), N2_En_B2_3f(J), N2_En_B3_3f(J)] = N2_Energy_f(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,J);   
        [N2_En_B1_2e(J), N2_En_B2_2e(J), N2_En_B3_2e(J)] = N2_Energy_e(N2_T_B_2,...
            N2_A_B_2,N2_B_B_2,N2_D_B_2,N2_l_B_2,N2_g_B_2,N2_o_B_2,N2_p_B_2,N2_q_B_2,J);        
        [N2_En_B1_2f(J), N2_En_B2_2f(J), N2_En_B3_2f(J)] = N2_Energy_f(N2_T_B_2,...
            N2_A_B_2,N2_B_B_2,N2_D_B_2,N2_l_B_2,N2_g_B_2,N2_o_B_2,N2_p_B_2,N2_q_B_2,J);  
        [N2_En_C1_3e(J), N2_En_C2_3e(J), N2_En_C3_3e(J)] = N2_Energy_e(N2_T_C,...
            N2_A_C,N2_B_C,N2_D_C,N2_l_C,N2_g_C,N2_o_C,N2_p_C,N2_q_C,J);
        [N2_En_C1_3f(J), N2_En_C2_3f(J), N2_En_C3_3f(J)] = N2_Energy_f(N2_T_C,...
            N2_A_C,N2_B_C,N2_D_C,N2_l_C,N2_g_C,N2_o_C,N2_p_C,N2_q_C,J); 
        [N2_En_C1_2e(J), N2_En_C2_2e(J), N2_En_C3_2e(J)] = N2_Energy_e(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,J);
        [N2_En_C1_2f(J), N2_En_C2_2f(J), N2_En_C3_2f(J)] = N2_Energy_f(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,J);    
        [N2_En_C1_4e(J), N2_En_C2_4e(J), N2_En_C3_4e(J)] = N2_Energy_e(N2_T_C_4,...
            N2_A_C_4,N2_B_C_4,N2_D_C_4,N2_l_C_4,N2_g_C_4,N2_o_C_4,N2_p_C_4,N2_q_C_4,J);
        [N2_En_C1_4f(J), N2_En_C2_4f(J), N2_En_C3_4f(J)] = N2_Energy_f(N2_T_C_4,...
            N2_A_C_4,N2_B_C_4,N2_D_C_4,N2_l_C_4,N2_g_C_4,N2_o_C_4,N2_p_C_4,N2_q_C_4,J);  
        [N2_En_C1_1e(J), N2_En_C2_1e(J), N2_En_C3_1e(J)] = N2_Energy_e(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,J);
        [N2_En_C1_1f(J), N2_En_C2_1f(J), N2_En_C3_1f(J)] = N2_Energy_f(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,J);           
        [N2_En_C1_0e(J), N2_En_C2_0e(J), N2_En_C3_0e(J)] = N2_Energy_e(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,J);
        [N2_En_C1_0f(J), N2_En_C2_0f(J), N2_En_C3_0f(J)] = N2_Energy_f(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,J); 
    end

% energy corrections due to A-Delta crossing in N2(+)        
    for N=1:N_J
        En_B1(N) = En_B1(N) + 0.5*5.45/(N-39+0.5);  
        En_B1(N) = En_B1(N) + 0.5*1.4/(N-43+0.5); 
        En_B1(N) = En_B1(N) + 0.5*2.96/(N-66+0.5);
        En_B1(N) = En_B1(N) + 0.5*0.93/(N-69+0.5);
        En_B1(N) = En_B1(N) + 0.5*2.5/(N-83+0.5);
        En_B1(N) = En_B1(N) + 0.5*0.76/(N-86+0.5);
        En_B1(N) = En_B1(N) + 0.5*1.49/(N-96+0.5);
        En_B1(N) = En_B1(N) + 0.5*0.85/(N-98+0.5);
        En_B2(N) = En_B2(N) + 0.5*1.58/(N-36+0.5);
        En_B2(N) = En_B2(N) + 0.5*1.82/(N-39+0.5);
        En_B2(N) = En_B2(N) + 0.5*0.81/(N-63+0.5);
        En_B2(N) = En_B2(N) + 0.5*3.37/(N-66+0.5);
        En_B2(N) = En_B2(N) + 0.5*0.06/(N-80+0.5);
        En_B2(N) = En_B2(N) + 0.5*2.67/(N-83+0.5);
        En_B2(N) = En_B2(N) + 0.5*0.14/(N-93+0.5);
        En_B2(N) = En_B2(N) + 0.5*1.28/(N-96+0.5);
    end

    for J=1:N_J
        if(J~=1)
            
            % N2(+) transition energies
            DEn_P1(J) = En_B1(J) - En_X1(J-1);
            DEn_P1_1(J) = En_B1_1(J) - En_X1_1(J-1);            
            DEn_P2(J) = En_B2(J) - En_X2(J-1);
            DEn_P2_1(J) = En_B2_1(J) - En_X2_1(J-1);
            DEn_R1(J) = En_B1(J-1) - En_X1(J);
            DEn_R1_1(J) = En_B1_1(J-1) - En_X1_1(J);
            DEn_R2(J) = En_B2(J-1) - En_X2(J);
            DEn_R2_1(J) = En_B2_1(J-1) - En_X2_1(J);
            DEn_PQ12(J) = En_B2(J) - En_X1(J-1);
            DEn_PQ12_1(J) = En_B2_1(J) - En_X1_1(J-1);
            DEn_RQ21(J) = En_B1(J-1) - En_X2(J);
            DEn_RQ21_1(J) = En_B1_1(J-1) - En_X2_1(J);
            
            % CN transition energies
            CN_DEn_P1(J) = CN_En_B1(J) - CN_En_X1(J-1);
            CN_DEn_P1_1(J) = CN_En_B1_1(J) - CN_En_X1_1(J-1); 
            CN_DEn_P1_2(J) = CN_En_B1_2(J) - CN_En_X1_2(J-1);
            CN_DEn_P1_3(J) = CN_En_B1_3(J) - CN_En_X1_3(J-1); 
            CN_DEn_P1_4(J) = CN_En_B1_4(J) - CN_En_X1_4(J-1); 
            CN_DEn_P1_5(J) = CN_En_B1_5(J) - CN_En_X1_5(J-1);
            CN_DEn_P1_6(J) = CN_En_B1_6(J) - CN_En_X1_6(J-1); 
            CN_DEn_P1_7(J) = CN_En_B1_7(J) - CN_En_X1_7(J-1);
            CN_DEn_P1_8(J) = CN_En_B1_8(J) - CN_En_X1_8(J-1);
            CN_DEn_P1_9(J) = CN_En_B1_9(J) - CN_En_X1_9(J-1);
            CN_DEn_P1_10(J) = CN_En_B1_10(J) - CN_En_X1_10(J-1);
            CN_DEn_P2(J) = CN_En_B2(J) - CN_En_X2(J-1);
            CN_DEn_P2_1(J) = CN_En_B2_1(J) - CN_En_X2_1(J-1);
            CN_DEn_P2_2(J) = CN_En_B2_2(J) - CN_En_X2_2(J-1);
            CN_DEn_P2_3(J) = CN_En_B2_3(J) - CN_En_X2_3(J-1);
            CN_DEn_P2_4(J) = CN_En_B2_4(J) - CN_En_X2_4(J-1);
            CN_DEn_P2_5(J) = CN_En_B2_5(J) - CN_En_X2_5(J-1);
            CN_DEn_P2_6(J) = CN_En_B2_6(J) - CN_En_X2_6(J-1);
            CN_DEn_P2_7(J) = CN_En_B2_7(J) - CN_En_X2_7(J-1);
            CN_DEn_P2_8(J) = CN_En_B2_8(J) - CN_En_X2_8(J-1);
            CN_DEn_P2_9(J) = CN_En_B2_9(J) - CN_En_X2_9(J-1);
            CN_DEn_P2_10(J) = CN_En_B2_10(J) - CN_En_X2_10(J-1);
            CN_DEn_R1(J) = CN_En_B1(J-1) - CN_En_X1(J);
            CN_DEn_R1_1(J) = CN_En_B1_1(J-1) - CN_En_X1_1(J);
            CN_DEn_R1_2(J) = CN_En_B1_2(J-1) - CN_En_X1_2(J);
            CN_DEn_R1_3(J) = CN_En_B1_3(J-1) - CN_En_X1_3(J);
            CN_DEn_R1_4(J) = CN_En_B1_4(J-1) - CN_En_X1_4(J);
            CN_DEn_R1_5(J) = CN_En_B1_5(J-1) - CN_En_X1_5(J);
            CN_DEn_R1_6(J) = CN_En_B1_6(J-1) - CN_En_X1_6(J);
            CN_DEn_R1_7(J) = CN_En_B1_7(J-1) - CN_En_X1_7(J);
            CN_DEn_R1_8(J) = CN_En_B1_8(J-1) - CN_En_X1_8(J);
            CN_DEn_R1_9(J) = CN_En_B1_9(J-1) - CN_En_X1_9(J);
            CN_DEn_R1_10(J) = CN_En_B1_10(J-1) - CN_En_X1_10(J);
            CN_DEn_R2(J) = CN_En_B2(J-1) - CN_En_X2(J);
            CN_DEn_R2_1(J) = CN_En_B2_1(J-1) - CN_En_X2_1(J);
            CN_DEn_R2_2(J) = CN_En_B2_2(J-1) - CN_En_X2_2(J);
            CN_DEn_R2_3(J) = CN_En_B2_3(J-1) - CN_En_X2_3(J);
            CN_DEn_R2_4(J) = CN_En_B2_4(J-1) - CN_En_X2_4(J);
            CN_DEn_R2_5(J) = CN_En_B2_5(J-1) - CN_En_X2_5(J);
            CN_DEn_R2_6(J) = CN_En_B2_6(J-1) - CN_En_X2_6(J);
            CN_DEn_R2_7(J) = CN_En_B2_7(J-1) - CN_En_X2_7(J);
            CN_DEn_R2_8(J) = CN_En_B2_8(J-1) - CN_En_X2_8(J);
            CN_DEn_R2_9(J) = CN_En_B2_9(J-1) - CN_En_X2_9(J);
            CN_DEn_R2_10(J) = CN_En_B2_10(J-1) - CN_En_X2_10(J);
            CN_DEn_PQ12(J) = CN_En_B2(J) - CN_En_X1(J-1);
            CN_DEn_PQ12_1(J) = CN_En_B2_1(J) - CN_En_X1_1(J-1);
            CN_DEn_PQ12_2(J) = CN_En_B2_2(J) - CN_En_X1_2(J-1);
            CN_DEn_PQ12_3(J) = CN_En_B2_3(J) - CN_En_X1_3(J-1);
            CN_DEn_PQ12_4(J) = CN_En_B2_4(J) - CN_En_X1_4(J-1);
            CN_DEn_PQ12_5(J) = CN_En_B2_5(J) - CN_En_X1_5(J-1);
            CN_DEn_PQ12_6(J) = CN_En_B2_6(J) - CN_En_X1_6(J-1);
            CN_DEn_PQ12_7(J) = CN_En_B2_7(J) - CN_En_X1_7(J-1);
            CN_DEn_PQ12_8(J) = CN_En_B2_8(J) - CN_En_X1_8(J-1);
            CN_DEn_PQ12_9(J) = CN_En_B2_9(J) - CN_En_X1_9(J-1);
            CN_DEn_PQ12_10(J) = CN_En_B2_10(J) - CN_En_X1_10(J-1);
            CN_DEn_RQ21(J) = CN_En_B1(J-1) - CN_En_X2(J);
            CN_DEn_RQ21_1(J) = CN_En_B1_1(J-1) - CN_En_X2_1(J);
            CN_DEn_RQ21_2(J) = CN_En_B1_2(J-1) - CN_En_X2_2(J);
            CN_DEn_RQ21_3(J) = CN_En_B1_3(J-1) - CN_En_X2_3(J);
            CN_DEn_RQ21_4(J) = CN_En_B1_4(J-1) - CN_En_X2_4(J);
            CN_DEn_RQ21_5(J) = CN_En_B1_5(J-1) - CN_En_X2_5(J);
            CN_DEn_RQ21_6(J) = CN_En_B1_6(J-1) - CN_En_X2_6(J);
            CN_DEn_RQ21_7(J) = CN_En_B1_7(J-1) - CN_En_X2_7(J);
            CN_DEn_RQ21_8(J) = CN_En_B1_8(J-1) - CN_En_X2_8(J);
            CN_DEn_RQ21_9(J) = CN_En_B1_9(J-1) - CN_En_X2_9(J);
            CN_DEn_RQ21_10(J) = CN_En_B1_10(J-1) - CN_En_X2_10(J);           
            
            % N2 transition energies (3-6), (2-5), (4-7), (1-4)
            % (0-2), (1-3), (2-4), (3-5)
            
            N2_DEn_R1_36e(J) = N2_En_C1_3e(J) - N2_En_B1_6e(J-1);
            N2_DEn_R2_36e(J) = N2_En_C2_3e(J) - N2_En_B2_6e(J-1);
            N2_DEn_R3_36e(J) = N2_En_C3_3e(J) - N2_En_B3_6e(J-1);
            N2_DEn_R1_36f(J) = N2_En_C1_3f(J) - N2_En_B1_6f(J-1);
            N2_DEn_R2_36f(J) = N2_En_C2_3f(J) - N2_En_B2_6f(J-1);
            N2_DEn_R3_36f(J) = N2_En_C3_3f(J) - N2_En_B3_6f(J-1);        
            N2_DEn_R1_25e(J) = N2_En_C1_2e(J) - N2_En_B1_5e(J-1);
            N2_DEn_R2_25e(J) = N2_En_C2_2e(J) - N2_En_B2_5e(J-1);
            N2_DEn_R3_25e(J) = N2_En_C3_2e(J) - N2_En_B3_5e(J-1);
            N2_DEn_R1_25f(J) = N2_En_C1_2f(J) - N2_En_B1_5f(J-1);
            N2_DEn_R2_25f(J) = N2_En_C2_2f(J) - N2_En_B2_5f(J-1);
            N2_DEn_R3_25f(J) = N2_En_C3_2f(J) - N2_En_B3_5f(J-1);           
            N2_DEn_R1_47e(J) = N2_En_C1_4e(J) - N2_En_B1_7e(J-1);
            N2_DEn_R2_47e(J) = N2_En_C2_4e(J) - N2_En_B2_7e(J-1);
            N2_DEn_R3_47e(J) = N2_En_C3_4e(J) - N2_En_B3_7e(J-1);
            N2_DEn_R1_47f(J) = N2_En_C1_4f(J) - N2_En_B1_7f(J-1);
            N2_DEn_R2_47f(J) = N2_En_C2_4f(J) - N2_En_B2_7f(J-1);
            N2_DEn_R3_47f(J) = N2_En_C3_4f(J) - N2_En_B3_7f(J-1);          
            N2_DEn_R1_14e(J) = N2_En_C1_1e(J) - N2_En_B1_4e(J-1);
            N2_DEn_R2_14e(J) = N2_En_C2_1e(J) - N2_En_B2_4e(J-1);
            N2_DEn_R3_14e(J) = N2_En_C3_1e(J) - N2_En_B3_4e(J-1); 
            N2_DEn_R1_14f(J) = N2_En_C1_1f(J) - N2_En_B1_4f(J-1);
            N2_DEn_R2_14f(J) = N2_En_C2_1f(J) - N2_En_B2_4f(J-1);
            N2_DEn_R3_14f(J) = N2_En_C3_1f(J) - N2_En_B3_4f(J-1);   
            N2_DEn_R1_03e(J) = N2_En_C1_0e(J) - N2_En_B1_3e(J-1);
            N2_DEn_R2_03e(J) = N2_En_C2_0e(J) - N2_En_B2_3e(J-1);
            N2_DEn_R3_03e(J) = N2_En_C3_0e(J) - N2_En_B3_3e(J-1); 
            N2_DEn_R1_03f(J) = N2_En_C1_0f(J) - N2_En_B1_3f(J-1);
            N2_DEn_R2_03f(J) = N2_En_C2_0f(J) - N2_En_B2_3f(J-1);
            N2_DEn_R3_03f(J) = N2_En_C3_0f(J) - N2_En_B3_3f(J-1);            
            N2_DEn_R1_02e(J) = N2_En_C1_0e(J) - N2_En_B1_2e(J-1);
            N2_DEn_R2_02e(J) = N2_En_C2_0e(J) - N2_En_B2_2e(J-1);
            N2_DEn_R3_02e(J) = N2_En_C3_0e(J) - N2_En_B3_2e(J-1); 
            N2_DEn_R1_02f(J) = N2_En_C1_0f(J) - N2_En_B1_2f(J-1);
            N2_DEn_R2_02f(J) = N2_En_C2_0f(J) - N2_En_B2_2f(J-1);
            N2_DEn_R3_02f(J) = N2_En_C3_0f(J) - N2_En_B3_2f(J-1);            
            N2_DEn_R1_13e(J) = N2_En_C1_1e(J) - N2_En_B1_3e(J-1);
            N2_DEn_R2_13e(J) = N2_En_C2_1e(J) - N2_En_B2_3e(J-1);
            N2_DEn_R3_13e(J) = N2_En_C3_1e(J) - N2_En_B3_3e(J-1); 
            N2_DEn_R1_13f(J) = N2_En_C1_1f(J) - N2_En_B1_3f(J-1);
            N2_DEn_R2_13f(J) = N2_En_C2_1f(J) - N2_En_B2_3f(J-1);
            N2_DEn_R3_13f(J) = N2_En_C3_1f(J) - N2_En_B3_3f(J-1);
            N2_DEn_R1_24e(J) = N2_En_C1_2e(J) - N2_En_B1_4e(J-1);
            N2_DEn_R2_24e(J) = N2_En_C2_2e(J) - N2_En_B2_4e(J-1);
            N2_DEn_R3_24e(J) = N2_En_C3_2e(J) - N2_En_B3_4e(J-1);
            N2_DEn_R1_24f(J) = N2_En_C1_2f(J) - N2_En_B1_4f(J-1);
            N2_DEn_R2_24f(J) = N2_En_C2_2f(J) - N2_En_B2_4f(J-1);
            N2_DEn_R3_24f(J) = N2_En_C3_2f(J) - N2_En_B3_4f(J-1);                  
            if(J >= 3)
                N2_DEn_P1_36e(J) = N2_En_C1_3e(J-1) - N2_En_B1_6e(J);
                N2_DEn_P1_25e(J) = N2_En_C1_2e(J-1) - N2_En_B1_5e(J);
                N2_DEn_P1_47e(J) = N2_En_C1_4e(J-1) - N2_En_B1_7e(J);
                N2_DEn_P1_14e(J) = N2_En_C1_1e(J-1) - N2_En_B1_4e(J);
                N2_DEn_P1_03e(J) = N2_En_C1_0e(J-1) - N2_En_B1_3e(J);
                N2_DEn_P1_36f(J) = N2_En_C1_3f(J-1) - N2_En_B1_6f(J);
                N2_DEn_P1_25f(J) = N2_En_C1_2f(J-1) - N2_En_B1_5f(J);
                N2_DEn_P1_47f(J) = N2_En_C1_4f(J-1) - N2_En_B1_7f(J);
                N2_DEn_P1_14f(J) = N2_En_C1_1f(J-1) - N2_En_B1_4f(J);
                N2_DEn_P1_03f(J) = N2_En_C1_0f(J-1) - N2_En_B1_3f(J);
                N2_DEn_P1_02f(J) = N2_En_C1_0f(J-1) - N2_En_B1_2f(J);
                N2_DEn_P1_13f(J) = N2_En_C1_1f(J-1) - N2_En_B1_3f(J);
                N2_DEn_P1_24f(J) = N2_En_C1_2f(J-1) - N2_En_B1_4f(J);
            end
            N2_DEn_P2_36e(J) = N2_En_C2_3e(J-1) - N2_En_B2_6e(J);
            N2_DEn_P3_36e(J) = N2_En_C3_3e(J-1) - N2_En_B3_6e(J);
            N2_DEn_P2_36f(J) = N2_En_C2_3f(J-1) - N2_En_B2_6f(J);
            N2_DEn_P3_36f(J) = N2_En_C3_3f(J-1) - N2_En_B3_6f(J);           
            N2_DEn_Q1_36e(J) = N2_En_C1_3f(J) - N2_En_B1_6e(J);
            N2_DEn_Q2_36e(J) = N2_En_C2_3f(J) - N2_En_B2_6e(J);
            N2_DEn_Q3_36e(J) = N2_En_C3_3f(J) - N2_En_B3_6e(J); 
            N2_DEn_Q1_36f(J) = N2_En_C1_3e(J) - N2_En_B1_6f(J);
            N2_DEn_Q2_36f(J) = N2_En_C2_3e(J) - N2_En_B2_6f(J);
            N2_DEn_Q3_36f(J) = N2_En_C3_3e(J) - N2_En_B3_6f(J);          
            N2_DEn_P2_25e(J) = N2_En_C2_2e(J-1) - N2_En_B2_5e(J);
            N2_DEn_P3_25e(J) = N2_En_C3_2e(J-1) - N2_En_B3_5e(J);
            N2_DEn_P2_25f(J) = N2_En_C2_2f(J-1) - N2_En_B2_5f(J);
            N2_DEn_P3_25f(J) = N2_En_C3_2f(J-1) - N2_En_B3_5f(J);         
            N2_DEn_Q1_25e(J) = N2_En_C1_2f(J) - N2_En_B1_5e(J);
            N2_DEn_Q2_25e(J) = N2_En_C2_2f(J) - N2_En_B2_5e(J);
            N2_DEn_Q3_25e(J) = N2_En_C3_2f(J) - N2_En_B3_5e(J);
            N2_DEn_Q1_25f(J) = N2_En_C1_2e(J) - N2_En_B1_5f(J);
            N2_DEn_Q2_25f(J) = N2_En_C2_2e(J) - N2_En_B2_5f(J);
            N2_DEn_Q3_25f(J) = N2_En_C3_2e(J) - N2_En_B3_5f(J);          
            N2_DEn_P2_47e(J) = N2_En_C2_4e(J-1) - N2_En_B2_7e(J);
            N2_DEn_P3_47e(J) = N2_En_C3_4e(J-1) - N2_En_B3_7e(J);
            N2_DEn_P2_47f(J) = N2_En_C2_4f(J-1) - N2_En_B2_7f(J);
            N2_DEn_P3_47f(J) = N2_En_C3_4f(J-1) - N2_En_B3_7f(J);          
            N2_DEn_Q1_47e(J) = N2_En_C1_4f(J) - N2_En_B1_7e(J);
            N2_DEn_Q2_47e(J) = N2_En_C2_4f(J) - N2_En_B2_7e(J);
            N2_DEn_Q3_47e(J) = N2_En_C3_4f(J) - N2_En_B3_7e(J);
            N2_DEn_Q1_47f(J) = N2_En_C1_4e(J) - N2_En_B1_7f(J);
            N2_DEn_Q2_47f(J) = N2_En_C2_4e(J) - N2_En_B2_7f(J);
            N2_DEn_Q3_47f(J) = N2_En_C3_4e(J) - N2_En_B3_7f(J);          
            N2_DEn_P2_14e(J) = N2_En_C2_1e(J-1) - N2_En_B2_4e(J);
            N2_DEn_P3_14e(J) = N2_En_C3_1e(J-1) - N2_En_B3_4e(J);
            N2_DEn_P2_14f(J) = N2_En_C2_1f(J-1) - N2_En_B2_4f(J);
            N2_DEn_P3_14f(J) = N2_En_C3_1f(J-1) - N2_En_B3_4f(J);          
            N2_DEn_Q1_14e(J) = N2_En_C1_1f(J) - N2_En_B1_4e(J);
            N2_DEn_Q2_14e(J) = N2_En_C2_1f(J) - N2_En_B2_4e(J);
            N2_DEn_Q3_14e(J) = N2_En_C3_1f(J) - N2_En_B3_4e(J);
            N2_DEn_Q1_14f(J) = N2_En_C1_1e(J) - N2_En_B1_4f(J);
            N2_DEn_Q2_14f(J) = N2_En_C2_1e(J) - N2_En_B2_4f(J);
            N2_DEn_Q3_14f(J) = N2_En_C3_1e(J) - N2_En_B3_4f(J);            
            N2_DEn_P2_03e(J) = N2_En_C2_0e(J-1) - N2_En_B2_3e(J);
            N2_DEn_P3_03e(J) = N2_En_C3_0e(J-1) - N2_En_B3_3e(J);
            N2_DEn_P2_03f(J) = N2_En_C2_0f(J-1) - N2_En_B2_3f(J);
            N2_DEn_P3_03f(J) = N2_En_C3_0f(J-1) - N2_En_B3_3f(J);          
            N2_DEn_Q1_03e(J) = N2_En_C1_0f(J) - N2_En_B1_3e(J);
            N2_DEn_Q2_03e(J) = N2_En_C2_0f(J) - N2_En_B2_3e(J);
            N2_DEn_Q3_03e(J) = N2_En_C3_0f(J) - N2_En_B3_3e(J);
            N2_DEn_Q1_03f(J) = N2_En_C1_0e(J) - N2_En_B1_3f(J);
            N2_DEn_Q2_03f(J) = N2_En_C2_0e(J) - N2_En_B2_3f(J);
            N2_DEn_Q3_03f(J) = N2_En_C3_0e(J) - N2_En_B3_3f(J);   
            N2_DEn_P2_02e(J) = N2_En_C2_0e(J-1) - N2_En_B2_2e(J);
            N2_DEn_P3_02e(J) = N2_En_C3_0e(J-1) - N2_En_B3_2e(J);
            N2_DEn_P2_02f(J) = N2_En_C2_0f(J-1) - N2_En_B2_2f(J);
            N2_DEn_P3_02f(J) = N2_En_C3_0f(J-1) - N2_En_B3_2f(J);          
            N2_DEn_Q1_02e(J) = N2_En_C1_0f(J) - N2_En_B1_2e(J);
            N2_DEn_Q2_02e(J) = N2_En_C2_0f(J) - N2_En_B2_2e(J);
            N2_DEn_Q3_02e(J) = N2_En_C3_0f(J) - N2_En_B3_2e(J);
            N2_DEn_Q1_02f(J) = N2_En_C1_0e(J) - N2_En_B1_2f(J);
            N2_DEn_Q2_02f(J) = N2_En_C2_0e(J) - N2_En_B2_2f(J);
            N2_DEn_Q3_02f(J) = N2_En_C3_0e(J) - N2_En_B3_2f(J); 
            N2_DEn_P2_13e(J) = N2_En_C2_1e(J-1) - N2_En_B2_3e(J);
            N2_DEn_P3_13e(J) = N2_En_C3_1e(J-1) - N2_En_B3_3e(J);
            N2_DEn_P2_13f(J) = N2_En_C2_1f(J-1) - N2_En_B2_3f(J);
            N2_DEn_P3_13f(J) = N2_En_C3_1f(J-1) - N2_En_B3_3f(J);          
            N2_DEn_Q1_13e(J) = N2_En_C1_1f(J) - N2_En_B1_3e(J);
            N2_DEn_Q2_13e(J) = N2_En_C2_1f(J) - N2_En_B2_3e(J);
            N2_DEn_Q3_13e(J) = N2_En_C3_1f(J) - N2_En_B3_3e(J);
            N2_DEn_Q1_13f(J) = N2_En_C1_1e(J) - N2_En_B1_3f(J);
            N2_DEn_Q2_13f(J) = N2_En_C2_1e(J) - N2_En_B2_3f(J);
            N2_DEn_Q3_13f(J) = N2_En_C3_1e(J) - N2_En_B3_3f(J);   
            N2_DEn_P2_24e(J) = N2_En_C2_2e(J-1) - N2_En_B2_4e(J);
            N2_DEn_P3_24e(J) = N2_En_C3_2e(J-1) - N2_En_B3_4e(J);
            N2_DEn_P2_24f(J) = N2_En_C2_2f(J-1) - N2_En_B2_4f(J);
            N2_DEn_P3_24f(J) = N2_En_C3_2f(J-1) - N2_En_B3_4f(J);         
            N2_DEn_Q1_24e(J) = N2_En_C1_2f(J) - N2_En_B1_4e(J);
            N2_DEn_Q2_24e(J) = N2_En_C2_2f(J) - N2_En_B2_4e(J);
            N2_DEn_Q3_24e(J) = N2_En_C3_2f(J) - N2_En_B3_4e(J);
            N2_DEn_Q1_24f(J) = N2_En_C1_2e(J) - N2_En_B1_4f(J);
            N2_DEn_Q2_24f(J) = N2_En_C2_2e(J) - N2_En_B2_4f(J);
            N2_DEn_Q3_24f(J) = N2_En_C3_2e(J) - N2_En_B3_4f(J);
            
        else
            
            % N2(+) transitions
            DEn_P1(1) = En_B1(1);
            DEn_P1_1(1) = En_B1_1(1);
            DEn_P2(1) = En_B2(1)+0.5*g_X;
            DEn_P2_1(1) = En_B2_1(1)-T_X1+0.5*g_X_1;
            DEn_R1(1) = T - En_X1(1);
            DEn_R1_1(1) = T_1 - En_X1_1(1);
            DEn_R2(1) = T - 0.5*g_B - En_X2(1);
            DEn_R2_1(1) = T_1 - 0.5*g_B_1 - En_X2(1);
            DEn_PQ12(1) = En_B2(1);
            DEn_PQ12_1(1) = En_B2_1(1)-T_X1;
            DEn_RQ21(1) = T - En_X2(1);
            DEn_RQ21_1(1) = T_1 - En_X2_1(1);
            
            % CN transitions
            CN_DEn_P1(1) = CN_En_B1(1);
            CN_DEn_P1_1(1) = CN_En_B1_1(1);
            CN_DEn_P1_2(1) = CN_En_B1_2(1);
            CN_DEn_P1_3(1) = CN_En_B1_3(1);
            CN_DEn_P1_4(1) = CN_En_B1_4(1);
            CN_DEn_P1_5(1) = CN_En_B1_5(1);
            CN_DEn_P1_6(1) = CN_En_B1_6(1);
            CN_DEn_P1_7(1) = CN_En_B1_7(1);
            CN_DEn_P1_8(1) = CN_En_B1_8(1);
            CN_DEn_P1_9(1) = CN_En_B1_9(1);
            CN_DEn_P1_10(1) = CN_En_B1_10(1);
            %CN_DEn_P2(1) = CN_En_B2(1)+0.5*CN_g_X;
            %CN_DEn_P2_1(1) = CN_En_B2_1(1)-CN_T_X1+0.5*CN_g_X_1;
            %CN_DEn_P2_2(1) = CN_En_B2_2(1)-CN_T_X2+0.5*CN_g_X_2;
            %CN_DEn_P2_3(1) = CN_En_B2_3(1)-CN_T_X3+0.5*CN_g_X_3;
            %CN_DEn_P2_4(1) = CN_En_B2_4(1)-CN_T_X4+0.5*CN_g_X_4;
            %N_DEn_P2_5(1) = CN_En_B2_5(1)-CN_T_X5+0.5*CN_g_X_5;
            %N_DEn_P2_6(1) = CN_En_B2_6(1)-CN_T_X6+0.5*CN_g_X_6;
            CN_DEn_PQ12(1) = CN_En_B2(1)-CN_T;
            CN_DEn_PQ12_1(1) = CN_En_B2_1(1)-CN_T_X1; 
            CN_DEn_PQ12_2(1) = CN_En_B2_2(1)-CN_T_X2; 
            CN_DEn_PQ12_3(1) = CN_En_B2_3(1)-CN_T_X3;
            CN_DEn_PQ12_4(1) = CN_En_B2_4(1)-CN_T_X4;
            CN_DEn_PQ12_5(1) = CN_En_B2_5(1)-CN_T_X5;
            CN_DEn_PQ12_6(1) = CN_En_B2_6(1)-CN_T_X6;
            CN_DEn_PQ12_7(1) = CN_En_B2_7(1)-CN_T_X7;
            CN_DEn_PQ12_8(1) = CN_En_B2_8(1)-CN_T_X8;
            CN_DEn_PQ12_9(1) = CN_En_B2_9(1)-CN_T_X9;
            CN_DEn_PQ12_10(1) = CN_En_B2_10(1)-CN_T_X10;
            
            CN_DEn_R1(1) = CN_T - CN_En_X1(1);
            CN_DEn_R1_1(1) = CN_T_1 - CN_En_X1_1(1);
            CN_DEn_R1_2(1) = CN_T_2 - CN_En_X1_2(1);
            CN_DEn_R1_3(1) = CN_T_3 - CN_En_X1_3(1);
            CN_DEn_R1_4(1) = CN_T_4 - CN_En_X1_4(1);
            CN_DEn_R1_5(1) = CN_T_5 - CN_En_X1_5(1);
            CN_DEn_R1_6(1) = CN_T_6 - CN_En_X1_6(1);
            CN_DEn_R1_7(1) = CN_T_7 - CN_En_X1_7(1);
            CN_DEn_R1_8(1) = CN_T_8 - CN_En_X1_8(1);
            CN_DEn_R1_9(1) = CN_T_9 - CN_En_X1_9(1);
            CN_DEn_R1_10(1) = CN_T_10 - CN_En_X1_10(1);
            %CN_DEn_R2(1) = CN_T-0.5*CN_g_B-CN_En_X2(1);
            %CN_DEn_R2_1(1) = CN_T_1-0.5*CN_g_B_1-CN_En_X2_1(1);
            %CN_DEn_R2_2(1) = CN_T_2-0.5*CN_g_B_2-CN_En_X2_2(1);
            %CN_DEn_R2_3(1) = CN_T_3-0.5*CN_g_B_3-CN_En_X2_3(1);
            %CN_DEn_R2_4(1) = CN_T_4-0.5*CN_g_B_4-CN_En_X2_4(1);
            %CN_DEn_R2_5(1) = CN_T_5-0.5*CN_g_B_5-CN_En_X2_5(1);
            %CN_DEn_R2_6(1) = CN_T_6-0.5*CN_g_B_6-CN_En_X2_6(1);
            CN_DEn_RQ21(1) = CN_T - CN_En_X2(1);
            CN_DEn_RQ21_1(1) = CN_T_1 - CN_En_X2_1(1); 
            CN_DEn_RQ21_2(1) = CN_T_2 - CN_En_X2_2(1); 
            CN_DEn_RQ21_3(1) = CN_T_3 - CN_En_X2_3(1);
            CN_DEn_RQ21_4(1) = CN_T_4 - CN_En_X2_4(1);
            CN_DEn_RQ21_5(1) = CN_T_5 - CN_En_X2_5(1);
            CN_DEn_RQ21_6(1) = CN_T_6 - CN_En_X2_6(1);
            CN_DEn_RQ21_7(1) = CN_T_7 - CN_En_X2_7(1);
            CN_DEn_RQ21_8(1) = CN_T_8 - CN_En_X2_8(1);
            CN_DEn_RQ21_9(1) = CN_T_9 - CN_En_X2_9(1);
            CN_DEn_RQ21_10(1) = CN_T_10 - CN_En_X2_10(1);
                       
            % N2 transitions
            %(3-6)
            
            [~ , ~ ,N2_En_B3_6_0e] = N2_Energy_e(N2_T_B,...
            N2_A_B,N2_B_B,N2_D_B,N2_l_B,N2_g_B,N2_o_B,N2_p_B,N2_q_B,0);
            [~ , ~ ,N2_En_B3_6_0f] = N2_Energy_f(N2_T_B,...
            N2_A_B,N2_B_B,N2_D_B,N2_l_B,N2_g_B,N2_o_B,N2_p_B,N2_q_B,0);
            
            [~ , ~ ,N2_En_C3_3_0e] = N2_Energy_e(N2_T_C,...
            N2_A_C,N2_B_C,N2_D_C,N2_l_C,N2_g_C,N2_o_C,N2_p_C,N2_q_C,0);
            [~ , ~ ,N2_En_C3_3_0f] = N2_Energy_f(N2_T_C,...
            N2_A_C,N2_B_C,N2_D_C,N2_l_C,N2_g_C,N2_o_C,N2_p_C,N2_q_C,0);    
        
            N2_DEn_R3_36e(1) = N2_En_C3_3e(1) - N2_En_B3_6_0e;
            N2_DEn_R3_36f(1) = N2_En_C3_3f(1) - N2_En_B3_6_0f;
            N2_DEn_P3_36e(1) = N2_En_C3_3_0e - N2_En_B3_6e(1);
            N2_DEn_P3_36f(1) = N2_En_C3_3_0f - N2_En_B3_6f(1);
            N2_DEn_Q2_36e(1) = N2_En_C2_3f(1) - N2_En_B2_6e(1);
            N2_DEn_Q2_36f(1) = N2_En_C2_3e(1) - N2_En_B2_6f(1);
            N2_DEn_Q3_36e(1) = N2_En_C3_3f(1) - N2_En_B3_6e(1);
            N2_DEn_Q3_36f(1) = N2_En_C3_3e(1) - N2_En_B3_6f(1);           
            
            %(2-5)
            
            [~ , ~ ,N2_En_B3_5_0e] = N2_Energy_e(N2_T_B_5,...
            N2_A_B_5,N2_B_B_5,N2_D_B_5,N2_l_B_5,N2_g_B_5,N2_o_B_5,N2_p_B_5,N2_q_B_5,0);
            [~ , ~ ,N2_En_B3_5_0f] = N2_Energy_f(N2_T_B_5,...
            N2_A_B_5,N2_B_B_5,N2_D_B_5,N2_l_B_5,N2_g_B_5,N2_o_B_5,N2_p_B_5,N2_q_B_5,0);
            
            [~ , ~ ,N2_En_C3_2_0e] = N2_Energy_e(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,0);
            [~ , ~ ,N2_En_C3_2_0f] = N2_Energy_f(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,0); 
            
            N2_DEn_R3_25e(1) = N2_En_C3_2e(1) - N2_En_B3_5_0e;
            N2_DEn_R3_25f(1) = N2_En_C3_2f(1) - N2_En_B3_5_0f;
            N2_DEn_P3_25e(1) = N2_En_C3_2_0e - N2_En_B3_5e(1);
            N2_DEn_P3_25f(1) = N2_En_C3_2_0f - N2_En_B3_5f(1);
            N2_DEn_Q2_25e(1) = N2_En_C2_2f(1) - N2_En_B2_5e(1);
            N2_DEn_Q2_25f(1) = N2_En_C2_2f(1) - N2_En_B2_5e(1);
            N2_DEn_Q3_25e(1) = N2_En_C3_2f(1) - N2_En_B3_5e(1);
            N2_DEn_Q3_25f(1) = N2_En_C3_2e(1) - N2_En_B3_5f(1);
            
            %(4-7)
            
            [~ , ~ ,N2_En_B3_7_0e] = N2_Energy_e(N2_T_B_7,...
            N2_A_B_7,N2_B_B_7,N2_D_B_7,N2_l_B_7,N2_g_B_7,N2_o_B_7,N2_p_B_7,N2_q_B_7,0);
            [~ , ~ ,N2_En_B3_7_0f] = N2_Energy_f(N2_T_B_7,...
            N2_A_B_7,N2_B_B_7,N2_D_B_7,N2_l_B_7,N2_g_B_7,N2_o_B_7,N2_p_B_7,N2_q_B_7,0);
            
            [~ , ~ ,N2_En_C3_4_0e] = N2_Energy_e(N2_T_C_4,...
            N2_A_C_4,N2_B_C_4,N2_D_C_4,N2_l_C_4,N2_g_C_4,N2_o_C_4,N2_p_C_4,N2_q_C_4,0);
            [~ , ~ ,N2_En_C3_4_0f] = N2_Energy_f(N2_T_C_4,...
            N2_A_C_4,N2_B_C_4,N2_D_C_4,N2_l_C_4,N2_g_C_4,N2_o_C_4,N2_p_C_4,N2_q_C_4,0);       
            
            N2_DEn_R3_47e(1) = N2_En_C3_4e(1) - N2_En_B3_7_0e;
            N2_DEn_R3_47f(1) = N2_En_C3_4f(1) - N2_En_B3_7_0f;
            N2_DEn_P3_47e(1) = N2_En_C3_4_0e - N2_En_B3_7e(1);
            N2_DEn_P3_47f(1) = N2_En_C3_4_0f - N2_En_B3_7f(1);
            N2_DEn_Q2_47e(1) = N2_En_C2_4f(1) - N2_En_B2_7e(1);
            N2_DEn_Q2_47f(1) = N2_En_C2_4e(1) - N2_En_B2_7f(1);
            N2_DEn_Q3_47e(1) = N2_En_C3_4f(1) - N2_En_B3_7e(1);
            N2_DEn_Q3_47f(1) = N2_En_C3_4e(1) - N2_En_B3_7f(1);
            
            %(1-4)
            
            [~ , ~ ,N2_En_B3_4_0e] = N2_Energy_e(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,0);
            [~ , ~ ,N2_En_B3_4_0f] = N2_Energy_f(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,0);
            
            [~ , ~ ,N2_En_C3_1_0e] = N2_Energy_e(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,0);
            [~ , ~ ,N2_En_C3_1_0f] = N2_Energy_f(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,0);  
            
            N2_DEn_R3_14e(1) = N2_En_C3_1e(1) - N2_En_B3_4_0e;
            N2_DEn_R3_14f(1) = N2_En_C3_1f(1) - N2_En_B3_4_0f;
            N2_DEn_P3_14e(1) = N2_En_C3_1_0e - N2_En_B3_4e(1);
            N2_DEn_P3_14f(1) = N2_En_C3_1_0f - N2_En_B3_4f(1);
            N2_DEn_Q2_14e(1) = N2_En_C2_1f(1) - N2_En_B2_4e(1);
            N2_DEn_Q2_14f(1) = N2_En_C2_1e(1) - N2_En_B2_4e(1);
            N2_DEn_Q3_14e(1) = N2_En_C3_1f(1) - N2_En_B3_4e(1); 
            N2_DEn_Q3_14f(1) = N2_En_C3_1e(1) - N2_En_B3_4f(1); 
            
             %(0-3)
            
            [~ , ~ ,N2_En_B3_3_0e] = N2_Energy_e(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,0);
            [~ , ~ ,N2_En_B3_3_0f] = N2_Energy_f(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,0);
            
            [~ , ~ ,N2_En_C3_0_0e] = N2_Energy_e(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,0);
            [~ , ~ ,N2_En_C3_0_0f] = N2_Energy_f(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,0);  
            
            N2_DEn_R3_03e(1) = N2_En_C3_0e(1) - N2_En_B3_3_0e;
            N2_DEn_R3_03f(1) = N2_En_C3_0f(1) - N2_En_B3_3_0f;
            N2_DEn_P3_03e(1) = N2_En_C3_0_0e - N2_En_B3_3e(1);
            N2_DEn_P3_03f(1) = N2_En_C3_0_0f - N2_En_B3_3f(1);
            N2_DEn_Q2_03e(1) = N2_En_C2_0f(1) - N2_En_B2_3e(1);
            N2_DEn_Q2_03f(1) = N2_En_C2_0e(1) - N2_En_B2_3e(1);
            N2_DEn_Q3_03e(1) = N2_En_C3_0f(1) - N2_En_B3_3e(1); 
            N2_DEn_Q3_03f(1) = N2_En_C3_0e(1) - N2_En_B3_3f(1);  
            
             %(0-2)
            
            [~ , ~ ,N2_En_B3_2_0e] = N2_Energy_e(N2_T_B_2,...
            N2_A_B_2,N2_B_B_2,N2_D_B_2,N2_l_B_2,N2_g_B_2,N2_o_B_2,N2_p_B_2,N2_q_B_2,0);
            [~ , ~ ,N2_En_B3_2_0f] = N2_Energy_f(N2_T_B_2,...
            N2_A_B_2,N2_B_B_2,N2_D_B_2,N2_l_B_2,N2_g_B_2,N2_o_B_2,N2_p_B_2,N2_q_B_2,0);
            
            [~ , ~ ,N2_En_C3_0_0e] = N2_Energy_e(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,0);
            [~ , ~ ,N2_En_C3_0_0f] = N2_Energy_f(N2_T_C_0,...
            N2_A_C_0,N2_B_C_0,N2_D_C_0,N2_l_C_0,N2_g_C_0,N2_o_C_0,N2_p_C_0,N2_q_C_0,0);  
            
            N2_DEn_R3_02e(1) = N2_En_C3_0e(1) - N2_En_B3_2_0e;
            N2_DEn_R3_02f(1) = N2_En_C3_0f(1) - N2_En_B3_2_0f;
            N2_DEn_P3_02e(1) = N2_En_C3_0_0e - N2_En_B3_2e(1);
            N2_DEn_P3_02f(1) = N2_En_C3_0_0f - N2_En_B3_2f(1);
            N2_DEn_Q2_02e(1) = N2_En_C2_0f(1) - N2_En_B2_2e(1);
            N2_DEn_Q2_02f(1) = N2_En_C2_0e(1) - N2_En_B2_2e(1);
            N2_DEn_Q3_02e(1) = N2_En_C3_0f(1) - N2_En_B3_2e(1); 
            N2_DEn_Q3_02f(1) = N2_En_C3_0e(1) - N2_En_B3_2f(1); 
            
             %(1-3)
            
            [~ , ~ ,N2_En_B3_3_0e] = N2_Energy_e(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,0);
            [~ , ~ ,N2_En_B3_3_0f] = N2_Energy_f(N2_T_B_3,...
            N2_A_B_3,N2_B_B_3,N2_D_B_3,N2_l_B_3,N2_g_B_3,N2_o_B_3,N2_p_B_3,N2_q_B_3,0);
            
            [~ , ~ ,N2_En_C3_1_0e] = N2_Energy_e(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,0);
            [~ , ~ ,N2_En_C3_1_0f] = N2_Energy_f(N2_T_C_1,...
            N2_A_C_1,N2_B_C_1,N2_D_C_1,N2_l_C_1,N2_g_C_1,N2_o_C_1,N2_p_C_1,N2_q_C_1,0);  
            
            N2_DEn_R3_13e(1) = N2_En_C3_1e(1) - N2_En_B3_3_0e;
            N2_DEn_R3_13f(1) = N2_En_C3_1f(1) - N2_En_B3_3_0f;
            N2_DEn_P3_13e(1) = N2_En_C3_1_0e - N2_En_B3_3e(1);
            N2_DEn_P3_13f(1) = N2_En_C3_1_0f - N2_En_B3_3f(1);
            N2_DEn_Q2_13e(1) = N2_En_C2_1f(1) - N2_En_B2_3e(1);
            N2_DEn_Q2_13f(1) = N2_En_C2_1e(1) - N2_En_B2_3e(1);
            N2_DEn_Q3_13e(1) = N2_En_C3_1f(1) - N2_En_B3_3e(1); 
            N2_DEn_Q3_13f(1) = N2_En_C3_1e(1) - N2_En_B3_3f(1);
           
             %(2-4)
            
            [~ , ~ ,N2_En_B3_4_0e] = N2_Energy_e(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,0);
            [~ , ~ ,N2_En_B3_4_0f] = N2_Energy_f(N2_T_B_4,...
            N2_A_B_4,N2_B_B_4,N2_D_B_4,N2_l_B_4,N2_g_B_4,N2_o_B_4,N2_p_B_4,N2_q_B_4,0);
            
            [~ , ~ ,N2_En_C3_2_0e] = N2_Energy_e(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,0);
            [~ , ~ ,N2_En_C3_2_0f] = N2_Energy_f(N2_T_C_2,...
            N2_A_C_2,N2_B_C_2,N2_D_C_2,N2_l_C_2,N2_g_C_2,N2_o_C_2,N2_p_C_2,N2_q_C_2,0); 
            
            N2_DEn_R3_24e(1) = N2_En_C3_2e(1) - N2_En_B3_4_0e;
            N2_DEn_R3_24f(1) = N2_En_C3_2f(1) - N2_En_B3_4_0f;
            N2_DEn_P3_24e(1) = N2_En_C3_2_0e - N2_En_B3_4e(1);
            N2_DEn_P3_24f(1) = N2_En_C3_2_0f - N2_En_B3_4f(1);
            N2_DEn_Q2_24e(1) = N2_En_C2_2f(1) - N2_En_B2_4e(1);
            N2_DEn_Q2_24f(1) = N2_En_C2_2f(1) - N2_En_B2_4e(1);
            N2_DEn_Q3_24e(1) = N2_En_C3_2f(1) - N2_En_B3_4e(1);
            N2_DEn_Q3_24f(1) = N2_En_C3_2e(1) - N2_En_B3_4f(1);
            
        end
    end 

%Tvib_ab = 1.2e+04;
modelo = @(beta,lambda) I_syntetic(beta,lambda);
beta0 = [4.00e+03 2.3e-01 delta_lambda_inst ...
    6.4 0.3 -0.5];
%beta0 = [4.00e+03 2.3e-01 delta_lambda_inst ...
%    6.4 0.96 -0.5]; 
plot(XData,feval(modelo,beta0,XData),'r');
hold off;
figure(2)
[I1, I2, I3, I4] =  I_syntetic_post(beta0, XData);
plot(XData, I1, '-b');
hold on;
plot(XData, I2, '-g');
hold on;
plot(XData, I3, '-k');
hold on;
plot(XData, I4, '-r');
hold on;
plot(XData,YData);
[beta, Res] = fitGM(XData, YData, 1, modelo, beta0, 30); 
Yerror = sqrt(Res.qui2/Res.ngl);
[beta, Res] = fitGM(XData, YData, Yerror, modelo, beta0, 30); 
I = I_syntetic(beta,XData);
figure(2)
plot(XData,I);
axis([376 395 0 1]);
figure(3)
plot(XData,I,'-r');
hold on;
plot(XData,YData);
hold on;
axis([376 395 0 1]);
[I1, I2, I3, ~] =  I_syntetic_post(beta, XData);
plot(XData, I1, '-b');
hold on;
plot(XData, I2, '-g');
hold on;
plot(XData, I3, '-k');
hold off;
save output
end


function [I1, I2, I3, I4] = I_syntetic_post(beta, lambda)
global En_B1 En_B2 DEn_R1 DEn_R2 DEn_P1 DEn_P2 DEn_RQ21 DEn_PQ12;
global En_B1_1 En_B2_1 DEn_R1_1 DEn_R2_1 DEn_P1_1 DEn_P2_1 DEn_RQ21_1 DEn_PQ12_1 DT_1;
global CN_En_B1 CN_En_B2 CN_DEn_R1 CN_DEn_R2 CN_DEn_P1 CN_DEn_P2 CN_DEn_RQ21 CN_DEn_PQ12;
global CN_En_B1_1 CN_En_B2_1 CN_DEn_R1_1 CN_DEn_R2_1 CN_DEn_P1_1 CN_DEn_P2_1 CN_DEn_RQ21_1 CN_DEn_PQ12_1;
global CN_En_B1_2 CN_En_B2_2 CN_DEn_R1_2 CN_DEn_R2_2 CN_DEn_P1_2 CN_DEn_P2_2 CN_DEn_RQ21_2 CN_DEn_PQ12_2;
global CN_En_B1_3 CN_En_B2_3 CN_DEn_R1_3 CN_DEn_R2_3 CN_DEn_P1_3 CN_DEn_P2_3 CN_DEn_RQ21_3 CN_DEn_PQ12_3;
global CN_En_B1_4 CN_En_B2_4 CN_DEn_R1_4 CN_DEn_R2_4 CN_DEn_P1_4 CN_DEn_P2_4 CN_DEn_RQ21_4 CN_DEn_PQ12_4;
global CN_En_B1_5 CN_En_B2_5 CN_DEn_R1_5 CN_DEn_R2_5 CN_DEn_P1_5 CN_DEn_P2_5 CN_DEn_RQ21_5 CN_DEn_PQ12_5;
global CN_En_B1_6 CN_En_B2_6 CN_DEn_R1_6 CN_DEn_R2_6 CN_DEn_P1_6 CN_DEn_P2_6 CN_DEn_RQ21_6 CN_DEn_PQ12_6;
global CN_En_B1_7 CN_En_B2_7 CN_DEn_R1_7 CN_DEn_R2_7 CN_DEn_P1_7 CN_DEn_P2_7 CN_DEn_RQ21_7 CN_DEn_PQ12_7;
global CN_En_B1_8 CN_En_B2_8 CN_DEn_R1_8 CN_DEn_R2_8 CN_DEn_P1_8 CN_DEn_P2_8 CN_DEn_RQ21_8 CN_DEn_PQ12_8;
global CN_En_B1_9 CN_En_B2_9 CN_DEn_R1_9 CN_DEn_R2_9 CN_DEn_P1_9 CN_DEn_P2_9 CN_DEn_RQ21_9 CN_DEn_PQ12_9;
global CN_En_B1_10 CN_En_B2_10 CN_DEn_R1_10 CN_DEn_R2_10 CN_DEn_P1_10 CN_DEn_P2_10 CN_DEn_RQ21_10 CN_DEn_PQ12_10;
global CN_DT1 CN_DT2 CN_DT3 CN_DT4 CN_DT5 CN_DT6 CN_DT7 CN_DT8 CN_DT9 CN_DT10; 
global delta_lambda_inst;
global N_J;
global h_p c kb;
global N2_Y_B N2_Y_B_4 N2_Y_B_5 N2_Y_B_7 N2_Y_C N2_Y_C_2 N2_Y_C_4 N2_Y_C_1 ...
    N2_Y_B_3 N2_Y_B_2 N2_Y_C_0;
global N2_En_C1_3e N2_En_C2_3e N2_En_C3_3e N2_En_C1_3f N2_En_C2_3f N2_En_C3_3f... 
    N2_DEn_R1_36e N2_DEn_R2_36e N2_DEn_R3_36e N2_DEn_P1_36e N2_DEn_P2_36e N2_DEn_P3_36e ...
    N2_DEn_Q1_36e N2_DEn_Q2_36e N2_DEn_Q3_36e ...
    N2_DEn_R1_36f N2_DEn_R2_36f N2_DEn_R3_36f N2_DEn_P1_36f N2_DEn_P2_36f N2_DEn_P3_36f ...
    N2_DEn_Q1_36f N2_DEn_Q2_36f N2_DEn_Q3_36f;
global N2_En_C1_2e N2_En_C2_2e N2_En_C3_2e N2_En_C1_2f N2_En_C2_2f N2_En_C3_2f ...
    N2_DEn_R1_25e N2_DEn_R2_25e N2_DEn_R3_25e N2_DEn_P1_25e N2_DEn_P2_25e N2_DEn_P3_25e ...
    N2_DEn_Q1_25e N2_DEn_Q2_25e N2_DEn_Q3_25e ...
    N2_DEn_R1_25f N2_DEn_R2_25f N2_DEn_R3_25f N2_DEn_P1_25f N2_DEn_P2_25f N2_DEn_P3_25f ...
    N2_DEn_Q1_25f N2_DEn_Q2_25f N2_DEn_Q3_25f;
global N2_En_C1_4e N2_En_C2_4e N2_En_C3_4e N2_En_C1_4f N2_En_C2_4f N2_En_C3_4f ...
    N2_DEn_R1_47e N2_DEn_R2_47e N2_DEn_R3_47e N2_DEn_P1_47e N2_DEn_P2_47e N2_DEn_P3_47e ...
    N2_DEn_Q1_47e N2_DEn_Q2_47e N2_DEn_Q3_47e ...
    N2_DEn_R1_47f N2_DEn_R2_47f N2_DEn_R3_47f N2_DEn_P1_47f N2_DEn_P2_47f N2_DEn_P3_47f ...
    N2_DEn_Q1_47f N2_DEn_Q2_47f N2_DEn_Q3_47f;
global N2_En_C1_1e N2_En_C2_1e N2_En_C3_1e N2_En_C1_1f N2_En_C2_1f N2_En_C3_1f ... 
    N2_DEn_R1_14e N2_DEn_R2_14e N2_DEn_R3_14e N2_DEn_P1_14e N2_DEn_P2_14e N2_DEn_P3_14e ...
    N2_DEn_Q1_14e N2_DEn_Q2_14e N2_DEn_Q3_14e ...
    N2_DEn_R1_14f N2_DEn_R2_14f N2_DEn_R3_14f N2_DEn_P1_14f N2_DEn_P2_14f N2_DEn_P3_14f ...
    N2_DEn_Q1_14f N2_DEn_Q2_14f N2_DEn_Q3_14f;
global N2_En_C1_0e N2_En_C2_0e N2_En_C3_0e N2_En_C1_0f N2_En_C2_0f N2_En_C3_0f ... 
    N2_DEn_R1_03e N2_DEn_R2_03e N2_DEn_R3_03e N2_DEn_P1_03e N2_DEn_P2_03e N2_DEn_P3_03e ...
    N2_DEn_Q1_03e N2_DEn_Q2_03e N2_DEn_Q3_03e ...
    N2_DEn_R1_03f N2_DEn_R2_03f N2_DEn_R3_03f N2_DEn_P1_03f N2_DEn_P2_03f N2_DEn_P3_03f ...
    N2_DEn_Q1_03f N2_DEn_Q2_03f N2_DEn_Q3_03f;
global N2_DEn_R1_02e N2_DEn_R2_02e N2_DEn_R3_02e N2_DEn_P1_02e N2_DEn_P2_02e N2_DEn_P3_02e ...
    N2_DEn_Q1_02e N2_DEn_Q2_02e N2_DEn_Q3_02e ...
    N2_DEn_R1_02f N2_DEn_R2_02f N2_DEn_R3_02f N2_DEn_P1_02f N2_DEn_P2_02f N2_DEn_P3_02f ...
    N2_DEn_Q1_02f N2_DEn_Q2_02f N2_DEn_Q3_02f;
global  N2_DEn_R1_13e N2_DEn_R2_13e N2_DEn_R3_13e N2_DEn_P1_13e N2_DEn_P2_13e N2_DEn_P3_13e ...
    N2_DEn_Q1_13e N2_DEn_Q2_13e N2_DEn_Q3_13e ...
    N2_DEn_R1_13f N2_DEn_R2_13f N2_DEn_R3_13f N2_DEn_P1_13f N2_DEn_P2_13f N2_DEn_P3_13f ...
    N2_DEn_Q1_13f N2_DEn_Q2_13f N2_DEn_Q3_13f;
global  N2_DEn_R1_24e N2_DEn_R2_24e N2_DEn_R3_24e N2_DEn_P1_24e N2_DEn_P2_24e N2_DEn_P3_24e ...
    N2_DEn_Q1_24e N2_DEn_Q2_24e N2_DEn_Q3_24e ...
    N2_DEn_R1_24f N2_DEn_R2_24f N2_DEn_R3_24f N2_DEn_P1_24f N2_DEn_P2_24f N2_DEn_P3_24f ...
    N2_DEn_Q1_24f N2_DEn_Q2_24f N2_DEn_Q3_24f;
global V_C_3 V_C_2 V_C_4 V_C_1 V_C_0;
global A_00 A_11;
global A_CN_00 A_CN_11 A_CN_22 A_CN_33 A_CN_44 A_CN_55 A_CN_66 ...
    A_CN_77 A_CN_88 A_CN_99 A_CN_1010;
global A_N2_36 A_N2_25 A_N2_47 A_N2_14 A_N2_03 A_N2_02 A_N2_13 A_N2_24; 

lambda = lambda(:);
Trot = beta(1);
Tvib = beta(1);
TrotCN = beta(1);
TrotN2 = beta(1);
TvibCN = beta(1);
TvibN2 = beta(1);
CN_ratio = 1;
N2p_ratio = beta(2);
N2_ratio = beta(4);
gamma = beta(5);
delta = beta(6);
eta = 0.0;
delta_lambda_inst = beta(3);
sum1 = zeros(length(lambda),1);
sum2 = zeros(length(lambda),1);
sum3 = zeros(length(lambda),1);
sum4 = zeros(length(lambda),1);
f_DT1 = -c*h_p*CN_DT1*1E2*(1-gamma + eta*CN_DT1)/kb/TvibCN+delta;
f_DT2 = -c*h_p*CN_DT2*1E2*(1-gamma + eta*CN_DT2)/kb/TvibCN+delta;
f_DT3 = -c*h_p*CN_DT3*1E2*(1-gamma + eta*CN_DT3)/kb/TvibCN+delta;
f_DT4 = -c*h_p*CN_DT4*1E2*(1-gamma + eta*CN_DT4)/kb/TvibCN+delta;
f_DT5 = -c*h_p*CN_DT5*1E2*(1-gamma + eta*CN_DT5)/kb/TvibCN+delta;
f_DT6 = -c*h_p*CN_DT6*1E2*(1-gamma + eta*CN_DT6)/kb/TvibCN+delta;
f_DT7 = -c*h_p*CN_DT7*1E2*(1-gamma + eta*CN_DT7)/kb/TvibCN+delta;
f_DT8 = -c*h_p*CN_DT8*1E2*(1-gamma + eta*CN_DT8)/kb/TvibCN+delta;
f_DT9 = -c*h_p*CN_DT9*1E2*(1-gamma + eta*CN_DT9)/kb/TvibCN+delta;
    for J=1:N_J
        if(mod(J,2)==1)
            gP = 2*(2/3)*(2*J+1);
            gR = 2*(1/3)*(2*(J-1)+1);
            gP_N2 = 3*(2/3)*(2*J+1);
            gR_N2 = 3*(1/3)*(2*(J-1)+1);
            gQ_N2 = 3*(2/3)*(2*J+1);
        else
            gP = 2*(1/3)*(2*J+1);
            gR = 2*(2/3)*(2*(J-1)+1);
            gP_N2 = 3*(1/3)*(2*J+1);
            gR_N2 = 3*(2/3)*(2*(J-1)+1);
            gQ_N2 = 3*(1/3)*(2*J+1);
        end
        gCN = 2*(2*J+1);
        
        % N2(+) frequencies and wavelengths
        nu_B1 = (DEn_R1(J))*1E2*c;
        nu_B2 = (DEn_R2(J))*1E2*c;
        nu_RQ21 = (DEn_RQ21(J))*1E2*c;
        nu_B1_1 = (DEn_R1_1(J))*1E2*c;
        nu_B2_1 = (DEn_R2_1(J))*1E2*c;
        nu_RQ21_1 = (DEn_RQ21_1(J))*1E2*c;
        lambda_B1 = 1E7/DEn_R1(J);
        lambda_B2 = 1E7/DEn_R2(J);
        lambda_RQ21 = 1E7/DEn_RQ21(J);
        lambda_B1_1 = 1E7/DEn_R1_1(J);
        lambda_B2_1 = 1E7/DEn_R2_1(J);
        lambda_RQ21_1 = 1E7/DEn_RQ21_1(J);
        
        % CN frequencies and wavelengths
        CN_nu_B1 = (CN_DEn_R1(J))*1E2*c;
        CN_nu_B2 = (CN_DEn_R2(J))*1E2*c;
        CN_nu_RQ21 = (CN_DEn_RQ21(J))*1E2*c;
        CN_nu_B1_1 = (CN_DEn_R1_1(J))*1E2*c;
        CN_nu_B1_2 = (CN_DEn_R1_2(J))*1E2*c;
        CN_nu_B1_3 = (CN_DEn_R1_3(J))*1E2*c;
        CN_nu_B1_4 = (CN_DEn_R1_4(J))*1E2*c;
        CN_nu_B1_5 = (CN_DEn_R1_5(J))*1E2*c;
        CN_nu_B1_6 = (CN_DEn_R1_6(J))*1E2*c;
        CN_nu_B1_7 = (CN_DEn_R1_7(J))*1E2*c;
        CN_nu_B1_8 = (CN_DEn_R1_8(J))*1E2*c;
        CN_nu_B1_9 = (CN_DEn_R1_9(J))*1E2*c;
        CN_nu_B1_10 = (CN_DEn_R1_10(J))*1E2*c;
        CN_nu_B2_1 = (CN_DEn_R2_1(J))*1E2*c;
        CN_nu_B2_2 = (CN_DEn_R2_2(J))*1E2*c;
        CN_nu_B2_3 = (CN_DEn_R2_3(J))*1E2*c;
        CN_nu_B2_4 = (CN_DEn_R2_4(J))*1E2*c;
        CN_nu_B2_5 = (CN_DEn_R2_5(J))*1E2*c;
        CN_nu_B2_6 = (CN_DEn_R2_6(J))*1E2*c;
        CN_nu_B2_7 = (CN_DEn_R2_7(J))*1E2*c;
        CN_nu_B2_8 = (CN_DEn_R2_8(J))*1E2*c;
        CN_nu_B2_9 = (CN_DEn_R2_9(J))*1E2*c;
        CN_nu_B2_10 = (CN_DEn_R2_10(J))*1E2*c;
        CN_nu_RQ21_1 = (CN_DEn_RQ21_1(J))*1E2*c;
        CN_nu_RQ21_2 = (CN_DEn_RQ21_2(J))*1E2*c;
        CN_nu_RQ21_3 = (CN_DEn_RQ21_3(J))*1E2*c;
        CN_nu_RQ21_4 = (CN_DEn_RQ21_4(J))*1E2*c;
        CN_nu_RQ21_5 = (CN_DEn_RQ21_5(J))*1E2*c;
        CN_nu_RQ21_6 = (CN_DEn_RQ21_6(J))*1E2*c;
        CN_nu_RQ21_7 = (CN_DEn_RQ21_7(J))*1E2*c;
        CN_nu_RQ21_8 = (CN_DEn_RQ21_8(J))*1E2*c;
        CN_nu_RQ21_9 = (CN_DEn_RQ21_9(J))*1E2*c;
        CN_nu_RQ21_10 = (CN_DEn_RQ21_10(J))*1E2*c;
        CN_lambda_B1 = 1E7/CN_DEn_R1(J);
        CN_lambda_B2 = 1E7/CN_DEn_R2(J);
        CN_lambda_RQ21 = 1E7/CN_DEn_RQ21(J);
        CN_lambda_B1_1 = 1E7/CN_DEn_R1_1(J);
        CN_lambda_B1_2 = 1E7/CN_DEn_R1_2(J);
        CN_lambda_B1_3 = 1E7/CN_DEn_R1_3(J);
        CN_lambda_B1_4 = 1E7/CN_DEn_R1_4(J);
        CN_lambda_B1_5 = 1E7/CN_DEn_R1_5(J);
        CN_lambda_B1_6 = 1E7/CN_DEn_R1_6(J);
        CN_lambda_B1_7 = 1E7/CN_DEn_R1_7(J);
        CN_lambda_B1_8 = 1E7/CN_DEn_R1_8(J);
        CN_lambda_B1_9 = 1E7/CN_DEn_R1_9(J);
        CN_lambda_B1_10 = 1E7/CN_DEn_R1_10(J);
        CN_lambda_B2_1 = 1E7/CN_DEn_R2_1(J);
        CN_lambda_B2_2 = 1E7/CN_DEn_R2_2(J);
        CN_lambda_B2_3 = 1E7/CN_DEn_R2_3(J);
        CN_lambda_B2_4 = 1E7/CN_DEn_R2_4(J);
        CN_lambda_B2_5 = 1E7/CN_DEn_R2_5(J);
        CN_lambda_B2_6 = 1E7/CN_DEn_R2_6(J);
        CN_lambda_B2_7 = 1E7/CN_DEn_R2_7(J);
        CN_lambda_B2_8 = 1E7/CN_DEn_R2_8(J);
        CN_lambda_B2_9 = 1E7/CN_DEn_R2_9(J);
        CN_lambda_B2_10 = 1E7/CN_DEn_R2_10(J);
        CN_lambda_RQ21_1 = 1E7/CN_DEn_RQ21_1(J);
        CN_lambda_RQ21_2 = 1E7/CN_DEn_RQ21_2(J);
        CN_lambda_RQ21_3 = 1E7/CN_DEn_RQ21_3(J);
        CN_lambda_RQ21_4 = 1E7/CN_DEn_RQ21_4(J);
        CN_lambda_RQ21_5 = 1E7/CN_DEn_RQ21_5(J);
        CN_lambda_RQ21_6 = 1E7/CN_DEn_RQ21_6(J);
        CN_lambda_RQ21_7 = 1E7/CN_DEn_RQ21_7(J);
        CN_lambda_RQ21_8 = 1E7/CN_DEn_RQ21_8(J);
        CN_lambda_RQ21_9 = 1E7/CN_DEn_RQ21_9(J);
        CN_lambda_RQ21_10 = 1E7/CN_DEn_RQ21_10(J);
      
            
        % N2(+) frequencies and wavelengths
        nu_P_B1 = (DEn_P1(J))*1E2*c;
        nu_P_B2 = (DEn_P2(J))*1E2*c;
        nu_PQ12 = (DEn_PQ12(J))*1E2*c;
        nu_P_B1_1 = (DEn_P1_1(J))*1E2*c;
        nu_P_B2_1 = (DEn_P2_1(J))*1E2*c;
        nu_PQ12_1 = (DEn_PQ12_1(J))*1E2*c;
        lambda_P_B1 = 1E7/DEn_P1(J);
        lambda_P_B2 = 1E7/DEn_P2(J);
        lambda_PQ12 = 1E7/DEn_PQ12(J);
        lambda_P_B1_1 = 1E7/DEn_P1_1(J);
        lambda_P_B2_1 = 1E7/DEn_P2_1(J);
        lambda_PQ12_1 = 1E7/DEn_PQ12_1(J);
            
        % CN frequencies and wavelengths
        CN_nu_P_B1 = (CN_DEn_P1(J))*1E2*c;
        CN_nu_P_B2 = (CN_DEn_P2(J))*1E2*c;
        CN_nu_PQ12 = (CN_DEn_PQ12(J))*1E2*c;
        CN_nu_P_B1_1 = (CN_DEn_P1_1(J))*1E2*c;
        CN_nu_P_B1_2 = (CN_DEn_P1_2(J))*1E2*c;
        CN_nu_P_B1_3 = (CN_DEn_P1_3(J))*1E2*c;
        CN_nu_P_B1_4 = (CN_DEn_P1_4(J))*1E2*c;
        CN_nu_P_B1_5 = (CN_DEn_P1_5(J))*1E2*c;
        CN_nu_P_B1_6 = (CN_DEn_P1_6(J))*1E2*c;
        CN_nu_P_B1_7 = (CN_DEn_P1_7(J))*1E2*c;
        CN_nu_P_B1_8 = (CN_DEn_P1_8(J))*1E2*c;
        CN_nu_P_B1_9 = (CN_DEn_P1_9(J))*1E2*c;
        CN_nu_P_B1_10 = (CN_DEn_P1_10(J))*1E2*c;
        CN_nu_P_B2_1 = (CN_DEn_P2_1(J))*1E2*c;
        CN_nu_P_B2_2 = (CN_DEn_P2_2(J))*1E2*c;
        CN_nu_P_B2_3 = (CN_DEn_P2_3(J))*1E2*c;
        CN_nu_P_B2_4 = (CN_DEn_P2_4(J))*1E2*c;
        CN_nu_P_B2_5 = (CN_DEn_P2_5(J))*1E2*c;
        CN_nu_P_B2_6 = (CN_DEn_P2_6(J))*1E2*c;
        CN_nu_P_B2_7 = (CN_DEn_P2_7(J))*1E2*c;
        CN_nu_P_B2_8 = (CN_DEn_P2_8(J))*1E2*c;
        CN_nu_P_B2_9 = (CN_DEn_P2_9(J))*1E2*c;
        CN_nu_P_B2_10 = (CN_DEn_P2_10(J))*1E2*c;
        CN_nu_PQ12_1 = (CN_DEn_PQ12_1(J))*1E2*c;
        CN_nu_PQ12_2 = (CN_DEn_PQ12_2(J))*1E2*c;
        CN_nu_PQ12_3 = (CN_DEn_PQ12_3(J))*1E2*c;
        CN_nu_PQ12_4 = (CN_DEn_PQ12_4(J))*1E2*c;
        CN_nu_PQ12_5 = (CN_DEn_PQ12_5(J))*1E2*c;
        CN_nu_PQ12_6 = (CN_DEn_PQ12_6(J))*1E2*c;
        CN_nu_PQ12_7 = (CN_DEn_PQ12_7(J))*1E2*c;
        CN_nu_PQ12_8 = (CN_DEn_PQ12_9(J))*1E2*c;
        CN_nu_PQ12_9 = (CN_DEn_PQ12_10(J))*1E2*c;
        CN_nu_PQ12_10 = (CN_DEn_PQ12_8(J))*1E2*c;
        CN_lambda_P_B1 = 1E7/CN_DEn_P1(J);
        CN_lambda_P_B2 = 1E7/CN_DEn_P2(J);
        CN_lambda_PQ12 = 1E7/CN_DEn_PQ12(J);
        CN_lambda_P_B1_1 = 1E7/CN_DEn_P1_1(J);
        CN_lambda_P_B1_2 = 1E7/CN_DEn_P1_2(J);
        CN_lambda_P_B1_3 = 1E7/CN_DEn_P1_3(J);
        CN_lambda_P_B1_4 = 1E7/CN_DEn_P1_4(J);
        CN_lambda_P_B1_5 = 1E7/CN_DEn_P1_5(J);
        CN_lambda_P_B1_6 = 1E7/CN_DEn_P1_6(J);
        CN_lambda_P_B1_7 = 1E7/CN_DEn_P1_7(J);
        CN_lambda_P_B1_8 = 1E7/CN_DEn_P1_8(J);
        CN_lambda_P_B1_9 = 1E7/CN_DEn_P1_9(J);
        CN_lambda_P_B1_10 = 1E7/CN_DEn_P1_10(J);
        CN_lambda_P_B2_1 = 1E7/CN_DEn_P2_1(J);
        CN_lambda_P_B2_2 = 1E7/CN_DEn_P2_2(J);
        CN_lambda_P_B2_3 = 1E7/CN_DEn_P2_3(J);
        CN_lambda_P_B2_4 = 1E7/CN_DEn_P2_4(J);
        CN_lambda_P_B2_5 = 1E7/CN_DEn_P2_5(J);
        CN_lambda_P_B2_6 = 1E7/CN_DEn_P2_6(J);
        CN_lambda_P_B2_7 = 1E7/CN_DEn_P2_7(J);
        CN_lambda_P_B2_8 = 1E7/CN_DEn_P2_8(J);
        CN_lambda_P_B2_9 = 1E7/CN_DEn_P2_9(J);
        CN_lambda_P_B2_10 = 1E7/CN_DEn_P2_10(J);
        CN_lambda_PQ12_1 = 1E7/CN_DEn_PQ12_1(J);
        CN_lambda_PQ12_2 = 1E7/CN_DEn_PQ12_2(J);  
        CN_lambda_PQ12_3 = 1E7/CN_DEn_PQ12_3(J);
        CN_lambda_PQ12_4 = 1E7/CN_DEn_PQ12_4(J);
        CN_lambda_PQ12_5 = 1E7/CN_DEn_PQ12_5(J);
        CN_lambda_PQ12_6 = 1E7/CN_DEn_PQ12_6(J);
        CN_lambda_PQ12_7 = 1E7/CN_DEn_PQ12_7(J); 
        CN_lambda_PQ12_8 = 1E7/CN_DEn_PQ12_8(J);
        CN_lambda_PQ12_9 = 1E7/CN_DEn_PQ12_9(J); 
        CN_lambda_PQ12_10 = 1E7/CN_DEn_PQ12_10(J); 
        
        % N2 frequencies and wavelengths
        % (3-6)
        N2_nu_R1e = (N2_DEn_R1_36e(J))*1E2*c;
        N2_nu_R2e = (N2_DEn_R2_36e(J))*1E2*c;
        N2_nu_R3e = (N2_DEn_R3_36e(J))*1E2*c;
        N2_nu_P1e = (N2_DEn_P1_36e(J))*1E2*c;
        N2_nu_P2e = (N2_DEn_P2_36e(J))*1E2*c;
        N2_nu_P3e = (N2_DEn_P3_36e(J))*1E2*c;
        N2_nu_Q1e = (N2_DEn_Q1_36e(J))*1E2*c;
        N2_nu_Q2e = (N2_DEn_Q2_36e(J))*1E2*c;
        N2_nu_Q3e = (N2_DEn_Q3_36e(J))*1E2*c;
        N2_nu_R1f = (N2_DEn_R1_36f(J))*1E2*c;
        N2_nu_R2f = (N2_DEn_R2_36f(J))*1E2*c;
        N2_nu_R3f = (N2_DEn_R3_36f(J))*1E2*c;
        N2_nu_P1f = (N2_DEn_P1_36f(J))*1E2*c;
        N2_nu_P2f = (N2_DEn_P2_36f(J))*1E2*c;
        N2_nu_P3f = (N2_DEn_P3_36f(J))*1E2*c;
        N2_nu_Q1f = (N2_DEn_Q1_36f(J))*1E2*c;
        N2_nu_Q2f = (N2_DEn_Q2_36f(J))*1E2*c;
        N2_nu_Q3f = (N2_DEn_Q3_36f(J))*1E2*c;     
        if(J>2)
            N2_lambda_Q1e = 1E7/N2_DEn_Q1_36e(J);
            N2_lambda_P1e = 1E7/N2_DEn_P1_36e(J);
            N2_lambda_R1e = 1E7/N2_DEn_R1_36e(J);
            N2_lambda_Q1f = 1E7/N2_DEn_Q1_36f(J);
            N2_lambda_P1f = 1E7/N2_DEn_P1_36f(J);
            N2_lambda_R1f = 1E7/N2_DEn_R1_36f(J);         
        end
        if(J>1)
            N2_lambda_Q2e = 1E7/N2_DEn_Q2_36e(J);
            N2_lambda_P2e = 1E7/N2_DEn_P2_36e(J);
            N2_lambda_R2e = 1E7/N2_DEn_R2_36e(J);
            N2_lambda_Q2f = 1E7/N2_DEn_Q2_36f(J);
            N2_lambda_P2f = 1E7/N2_DEn_P2_36f(J);
            N2_lambda_R2f = 1E7/N2_DEn_R2_36f(J);          
        end
        N2_lambda_Q3e = 1E7/N2_DEn_Q3_36e(J);            
        N2_lambda_P3e = 1E7/N2_DEn_P3_36e(J);
        N2_lambda_R3e = 1E7/N2_DEn_R3_36e(J);
        N2_lambda_Q3f = 1E7/N2_DEn_Q3_36f(J);            
        N2_lambda_P3f = 1E7/N2_DEn_P3_36f(J);
        N2_lambda_R3f = 1E7/N2_DEn_R3_36f(J);     
        
        % (2-5)
        N2_nu_R1_25e = (N2_DEn_R1_25e(J))*1E2*c;
        N2_nu_R2_25e = (N2_DEn_R2_25e(J))*1E2*c;
        N2_nu_R3_25e = (N2_DEn_R3_25e(J))*1E2*c;
        N2_nu_P1_25e = (N2_DEn_P1_25e(J))*1E2*c;
        N2_nu_P2_25e = (N2_DEn_P2_25e(J))*1E2*c;
        N2_nu_P3_25e = (N2_DEn_P3_25e(J))*1E2*c;
        N2_nu_Q1_25e = (N2_DEn_Q1_25e(J))*1E2*c;
        N2_nu_Q2_25e = (N2_DEn_Q2_25e(J))*1E2*c;
        N2_nu_Q3_25e = (N2_DEn_Q3_25e(J))*1E2*c;
        N2_nu_R1_25f = (N2_DEn_R1_25f(J))*1E2*c;
        N2_nu_R2_25f = (N2_DEn_R2_25f(J))*1E2*c;
        N2_nu_R3_25f = (N2_DEn_R3_25f(J))*1E2*c;
        N2_nu_P1_25f = (N2_DEn_P1_25f(J))*1E2*c;
        N2_nu_P2_25f = (N2_DEn_P2_25f(J))*1E2*c;
        N2_nu_P3_25f = (N2_DEn_P3_25f(J))*1E2*c;
        N2_nu_Q1_25f = (N2_DEn_Q1_25f(J))*1E2*c;
        N2_nu_Q2_25f = (N2_DEn_Q2_25f(J))*1E2*c;
        N2_nu_Q3_25f = (N2_DEn_Q3_25f(J))*1E2*c;  
        if(J>2)
            N2_lambda_Q1_25e = 1E7/N2_DEn_Q1_25e(J);
            N2_lambda_P1_25e = 1E7/N2_DEn_P1_25e(J);
            N2_lambda_R1_25e = 1E7/N2_DEn_R1_25e(J);
            N2_lambda_Q1_25f = 1E7/N2_DEn_Q1_25f(J);
            N2_lambda_P1_25f = 1E7/N2_DEn_P1_25f(J);
            N2_lambda_R1_25f = 1E7/N2_DEn_R1_25f(J);           
        end
        if(J>1)
            N2_lambda_Q2_25e = 1E7/N2_DEn_Q2_25e(J);
            N2_lambda_P2_25e = 1E7/N2_DEn_P2_25e(J);
            N2_lambda_R2_25e = 1E7/N2_DEn_R2_25e(J);
            N2_lambda_Q2_25f = 1E7/N2_DEn_Q2_25f(J);
            N2_lambda_P2_25f = 1E7/N2_DEn_P2_25f(J);
            N2_lambda_R2_25f = 1E7/N2_DEn_R2_25f(J);          
        end
        N2_lambda_Q3_25e = 1E7/N2_DEn_Q3_25e(J);            
        N2_lambda_P3_25e = 1E7/N2_DEn_P3_25e(J);
        N2_lambda_R3_25e = 1E7/N2_DEn_R3_25e(J);
        N2_lambda_Q3_25f = 1E7/N2_DEn_Q3_25f(J);            
        N2_lambda_P3_25f = 1E7/N2_DEn_P3_25f(J);
        N2_lambda_R3_25f = 1E7/N2_DEn_R3_25f(J);      
        
        % (4-7)
        N2_nu_R1_47e = (N2_DEn_R1_47e(J))*1E2*c;
        N2_nu_R2_47e = (N2_DEn_R2_47e(J))*1E2*c;
        N2_nu_R3_47e = (N2_DEn_R3_47e(J))*1E2*c;
        N2_nu_P1_47e = (N2_DEn_P1_47e(J))*1E2*c;
        N2_nu_P2_47e = (N2_DEn_P2_47e(J))*1E2*c;
        N2_nu_P3_47e = (N2_DEn_P3_47e(J))*1E2*c;
        N2_nu_Q1_47e = (N2_DEn_Q1_47e(J))*1E2*c;
        N2_nu_Q2_47e = (N2_DEn_Q2_47e(J))*1E2*c;
        N2_nu_Q3_47e = (N2_DEn_Q3_47e(J))*1E2*c;
        N2_nu_R1_47f = (N2_DEn_R1_47f(J))*1E2*c;
        N2_nu_R2_47f = (N2_DEn_R2_47f(J))*1E2*c;
        N2_nu_R3_47f = (N2_DEn_R3_47f(J))*1E2*c;
        N2_nu_P1_47f = (N2_DEn_P1_47f(J))*1E2*c;
        N2_nu_P2_47f = (N2_DEn_P2_47f(J))*1E2*c;
        N2_nu_P3_47f = (N2_DEn_P3_47f(J))*1E2*c;
        N2_nu_Q1_47f = (N2_DEn_Q1_47f(J))*1E2*c;
        N2_nu_Q2_47f = (N2_DEn_Q2_47f(J))*1E2*c;
        N2_nu_Q3_47f = (N2_DEn_Q3_47f(J))*1E2*c;    
        if(J>2)
            N2_lambda_Q1_47e = 1E7/N2_DEn_Q1_47e(J);
            N2_lambda_P1_47e = 1E7/N2_DEn_P1_47e(J);
            N2_lambda_R1_47e = 1E7/N2_DEn_R1_47e(J);
            N2_lambda_Q1_47f = 1E7/N2_DEn_Q1_47f(J);
            N2_lambda_P1_47f = 1E7/N2_DEn_P1_47f(J);
            N2_lambda_R1_47f = 1E7/N2_DEn_R1_47f(J);           
        end
        if(J>1)
            N2_lambda_Q2_47e = 1E7/N2_DEn_Q2_47e(J);
            N2_lambda_P2_47e = 1E7/N2_DEn_P2_47e(J);
            N2_lambda_R2_47e = 1E7/N2_DEn_R2_47e(J);
            N2_lambda_Q2_47f = 1E7/N2_DEn_Q2_47f(J);
            N2_lambda_P2_47f = 1E7/N2_DEn_P2_47f(J);
            N2_lambda_R2_47f = 1E7/N2_DEn_R2_47f(J);         
        end
        N2_lambda_Q3_47e = 1E7/N2_DEn_Q3_47e(J);            
        N2_lambda_P3_47e = 1E7/N2_DEn_P3_47e(J);
        N2_lambda_R3_47e = 1E7/N2_DEn_R3_47e(J);
        N2_lambda_Q3_47f = 1E7/N2_DEn_Q3_47f(J);            
        N2_lambda_P3_47f = 1E7/N2_DEn_P3_47f(J);
        N2_lambda_R3_47f = 1E7/N2_DEn_R3_47f(J);     

        % (1-4)
        N2_nu_R1_14e = (N2_DEn_R1_14e(J))*1E2*c;
        N2_nu_R2_14e = (N2_DEn_R2_14e(J))*1E2*c;
        N2_nu_R3_14e = (N2_DEn_R3_14e(J))*1E2*c;
        N2_nu_P1_14e = (N2_DEn_P1_14e(J))*1E2*c;
        N2_nu_P2_14e = (N2_DEn_P2_14e(J))*1E2*c;
        N2_nu_P3_14e = (N2_DEn_P3_14e(J))*1E2*c;
        N2_nu_Q1_14e = (N2_DEn_Q1_14e(J))*1E2*c;
        N2_nu_Q2_14e = (N2_DEn_Q2_14e(J))*1E2*c;
        N2_nu_Q3_14e = (N2_DEn_Q3_14e(J))*1E2*c;
        N2_nu_R1_14f = (N2_DEn_R1_14f(J))*1E2*c;
        N2_nu_R2_14f = (N2_DEn_R2_14f(J))*1E2*c;
        N2_nu_R3_14f = (N2_DEn_R3_14f(J))*1E2*c;
        N2_nu_P1_14f = (N2_DEn_P1_14f(J))*1E2*c;
        N2_nu_P2_14f = (N2_DEn_P2_14f(J))*1E2*c;
        N2_nu_P3_14f = (N2_DEn_P3_14f(J))*1E2*c;
        N2_nu_Q1_14f = (N2_DEn_Q1_14f(J))*1E2*c;
        N2_nu_Q2_14f = (N2_DEn_Q2_14f(J))*1E2*c;
        N2_nu_Q3_14f = (N2_DEn_Q3_14f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_14e = 1E7/N2_DEn_Q1_14e(J);
            N2_lambda_P1_14e = 1E7/N2_DEn_P1_14e(J);
            N2_lambda_R1_14e = 1E7/N2_DEn_R1_14e(J);
            N2_lambda_Q1_14f = 1E7/N2_DEn_Q1_14f(J);
            N2_lambda_P1_14f = 1E7/N2_DEn_P1_14f(J);
            N2_lambda_R1_14f = 1E7/N2_DEn_R1_14f(J);           
        end
        if(J>1)
            N2_lambda_Q2_14e = 1E7/N2_DEn_Q2_14e(J);
            N2_lambda_P2_14e = 1E7/N2_DEn_P2_14e(J);
            N2_lambda_R2_14e = 1E7/N2_DEn_R2_14e(J);
            N2_lambda_Q2_14f = 1E7/N2_DEn_Q2_14f(J);
            N2_lambda_P2_14f = 1E7/N2_DEn_P2_14f(J);
            N2_lambda_R2_14f = 1E7/N2_DEn_R2_14f(J);         
        end
        N2_lambda_Q3_14e = 1E7/N2_DEn_Q3_14e(J);            
        N2_lambda_P3_14e = 1E7/N2_DEn_P3_14e(J);
        N2_lambda_R3_14e = 1E7/N2_DEn_R3_14e(J);
        N2_lambda_Q3_14f = 1E7/N2_DEn_Q3_14f(J);            
        N2_lambda_P3_14f = 1E7/N2_DEn_P3_14f(J);
        N2_lambda_R3_14f = 1E7/N2_DEn_R3_14f(J);
        
        % (0-3)
        N2_nu_R1_03e = (N2_DEn_R1_03e(J))*1E2*c;
        N2_nu_R2_03e = (N2_DEn_R2_03e(J))*1E2*c;
        N2_nu_R3_03e = (N2_DEn_R3_03e(J))*1E2*c;
        N2_nu_P1_03e = (N2_DEn_P1_03e(J))*1E2*c;
        N2_nu_P2_03e = (N2_DEn_P2_03e(J))*1E2*c;
        N2_nu_P3_03e = (N2_DEn_P3_03e(J))*1E2*c;
        N2_nu_Q1_03e = (N2_DEn_Q1_03e(J))*1E2*c;
        N2_nu_Q2_03e = (N2_DEn_Q2_03e(J))*1E2*c;
        N2_nu_Q3_03e = (N2_DEn_Q3_03e(J))*1E2*c;
        N2_nu_R1_03f = (N2_DEn_R1_03f(J))*1E2*c;
        N2_nu_R2_03f = (N2_DEn_R2_03f(J))*1E2*c;
        N2_nu_R3_03f = (N2_DEn_R3_03f(J))*1E2*c;
        N2_nu_P1_03f = (N2_DEn_P1_03f(J))*1E2*c;
        N2_nu_P2_03f = (N2_DEn_P2_03f(J))*1E2*c;
        N2_nu_P3_03f = (N2_DEn_P3_03f(J))*1E2*c;
        N2_nu_Q1_03f = (N2_DEn_Q1_03f(J))*1E2*c;
        N2_nu_Q2_03f = (N2_DEn_Q2_03f(J))*1E2*c;
        N2_nu_Q3_03f = (N2_DEn_Q3_03f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_03e = 1E7/N2_DEn_Q1_03e(J);
            N2_lambda_P1_03e = 1E7/N2_DEn_P1_03e(J);
            N2_lambda_R1_03e = 1E7/N2_DEn_R1_03e(J);
            N2_lambda_Q1_03f = 1E7/N2_DEn_Q1_03f(J);
            N2_lambda_P1_03f = 1E7/N2_DEn_P1_03f(J);
            N2_lambda_R1_03f = 1E7/N2_DEn_R1_03f(J);           
        end
        if(J>1)
            N2_lambda_Q2_03e = 1E7/N2_DEn_Q2_03e(J);
            N2_lambda_P2_03e = 1E7/N2_DEn_P2_03e(J);
            N2_lambda_R2_03e = 1E7/N2_DEn_R2_03e(J);
            N2_lambda_Q2_03f = 1E7/N2_DEn_Q2_03f(J);
            N2_lambda_P2_03f = 1E7/N2_DEn_P2_03f(J);
            N2_lambda_R2_03f = 1E7/N2_DEn_R2_03f(J);         
        end
        N2_lambda_Q3_03e = 1E7/N2_DEn_Q3_03e(J);            
        N2_lambda_P3_03e = 1E7/N2_DEn_P3_03e(J);
        N2_lambda_R3_03e = 1E7/N2_DEn_R3_03e(J);
        N2_lambda_Q3_03f = 1E7/N2_DEn_Q3_03f(J);            
        N2_lambda_P3_03f = 1E7/N2_DEn_P3_03f(J);
        N2_lambda_R3_03f = 1E7/N2_DEn_R3_03f(J);      
        
        % (0-2)
        N2_nu_R1_02e = (N2_DEn_R1_02e(J))*1E2*c;
        N2_nu_R2_02e = (N2_DEn_R2_02e(J))*1E2*c;
        N2_nu_R3_02e = (N2_DEn_R3_02e(J))*1E2*c;
        N2_nu_P1_02e = (N2_DEn_P1_02e(J))*1E2*c;
        N2_nu_P2_02e = (N2_DEn_P2_02e(J))*1E2*c;
        N2_nu_P3_02e = (N2_DEn_P3_02e(J))*1E2*c;
        N2_nu_Q1_02e = (N2_DEn_Q1_02e(J))*1E2*c;
        N2_nu_Q2_02e = (N2_DEn_Q2_02e(J))*1E2*c;
        N2_nu_Q3_02e = (N2_DEn_Q3_02e(J))*1E2*c;
        N2_nu_R1_02f = (N2_DEn_R1_02f(J))*1E2*c;
        N2_nu_R2_02f = (N2_DEn_R2_02f(J))*1E2*c;
        N2_nu_R3_02f = (N2_DEn_R3_02f(J))*1E2*c;
        N2_nu_P1_02f = (N2_DEn_P1_02f(J))*1E2*c;
        N2_nu_P2_02f = (N2_DEn_P2_02f(J))*1E2*c;
        N2_nu_P3_02f = (N2_DEn_P3_02f(J))*1E2*c;
        N2_nu_Q1_02f = (N2_DEn_Q1_02f(J))*1E2*c;
        N2_nu_Q2_02f = (N2_DEn_Q2_02f(J))*1E2*c;
        N2_nu_Q3_02f = (N2_DEn_Q3_02f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_02e = 1E7/N2_DEn_Q1_02e(J);
            N2_lambda_P1_02e = 1E7/N2_DEn_P1_02e(J);
            N2_lambda_R1_02e = 1E7/N2_DEn_R1_02e(J);
            N2_lambda_Q1_02f = 1E7/N2_DEn_Q1_02f(J);
            N2_lambda_P1_02f = 1E7/N2_DEn_P1_02f(J);
            N2_lambda_R1_02f = 1E7/N2_DEn_R1_02f(J);           
        end
        if(J>1)
            N2_lambda_Q2_02e = 1E7/N2_DEn_Q2_02e(J);
            N2_lambda_P2_02e = 1E7/N2_DEn_P2_02e(J);
            N2_lambda_R2_02e = 1E7/N2_DEn_R2_02e(J);
            N2_lambda_Q2_02f = 1E7/N2_DEn_Q2_02f(J);
            N2_lambda_P2_02f = 1E7/N2_DEn_P2_02f(J);
            N2_lambda_R2_02f = 1E7/N2_DEn_R2_02f(J);         
        end
        N2_lambda_Q3_02e = 1E7/N2_DEn_Q3_02e(J);            
        N2_lambda_P3_02e = 1E7/N2_DEn_P3_02e(J);
        N2_lambda_R3_02e = 1E7/N2_DEn_R3_02e(J);
        N2_lambda_Q3_02f = 1E7/N2_DEn_Q3_02f(J);            
        N2_lambda_P3_02f = 1E7/N2_DEn_P3_02f(J);
        N2_lambda_R3_02f = 1E7/N2_DEn_R3_02f(J);         
        
        % (1-3)
        N2_nu_R1_13e = (N2_DEn_R1_13e(J))*1E2*c;
        N2_nu_R2_13e = (N2_DEn_R2_13e(J))*1E2*c;
        N2_nu_R3_13e = (N2_DEn_R3_13e(J))*1E2*c;
        N2_nu_P1_13e = (N2_DEn_P1_13e(J))*1E2*c;
        N2_nu_P2_13e = (N2_DEn_P2_13e(J))*1E2*c;
        N2_nu_P3_13e = (N2_DEn_P3_13e(J))*1E2*c;
        N2_nu_Q1_13e = (N2_DEn_Q1_13e(J))*1E2*c;
        N2_nu_Q2_13e = (N2_DEn_Q2_13e(J))*1E2*c;
        N2_nu_Q3_13e = (N2_DEn_Q3_13e(J))*1E2*c;
        N2_nu_R1_13f = (N2_DEn_R1_13f(J))*1E2*c;
        N2_nu_R2_13f = (N2_DEn_R2_13f(J))*1E2*c;
        N2_nu_R3_13f = (N2_DEn_R3_13f(J))*1E2*c;
        N2_nu_P1_13f = (N2_DEn_P1_13f(J))*1E2*c;
        N2_nu_P2_13f = (N2_DEn_P2_13f(J))*1E2*c;
        N2_nu_P3_13f = (N2_DEn_P3_13f(J))*1E2*c;
        N2_nu_Q1_13f = (N2_DEn_Q1_13f(J))*1E2*c;
        N2_nu_Q2_13f = (N2_DEn_Q2_13f(J))*1E2*c;
        N2_nu_Q3_13f = (N2_DEn_Q3_13f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_13e = 1E7/N2_DEn_Q1_13e(J);
            N2_lambda_P1_13e = 1E7/N2_DEn_P1_13e(J);
            N2_lambda_R1_13e = 1E7/N2_DEn_R1_13e(J);
            N2_lambda_Q1_13f = 1E7/N2_DEn_Q1_13f(J);
            N2_lambda_P1_13f = 1E7/N2_DEn_P1_13f(J);
            N2_lambda_R1_13f = 1E7/N2_DEn_R1_13f(J);           
        end
        if(J>1)
            N2_lambda_Q2_13e = 1E7/N2_DEn_Q2_13e(J);
            N2_lambda_P2_13e = 1E7/N2_DEn_P2_13e(J);
            N2_lambda_R2_13e = 1E7/N2_DEn_R2_13e(J);
            N2_lambda_Q2_13f = 1E7/N2_DEn_Q2_13f(J);
            N2_lambda_P2_13f = 1E7/N2_DEn_P2_13f(J);
            N2_lambda_R2_13f = 1E7/N2_DEn_R2_13f(J);         
        end
        N2_lambda_Q3_13e = 1E7/N2_DEn_Q3_13e(J);            
        N2_lambda_P3_13e = 1E7/N2_DEn_P3_13e(J);
        N2_lambda_R3_13e = 1E7/N2_DEn_R3_13e(J);
        N2_lambda_Q3_13f = 1E7/N2_DEn_Q3_13f(J);            
        N2_lambda_P3_13f = 1E7/N2_DEn_P3_13f(J);
        N2_lambda_R3_13f = 1E7/N2_DEn_R3_13f(J);        
        
        % (2-4)
        N2_nu_R1_24e = (N2_DEn_R1_24e(J))*1E2*c;
        N2_nu_R2_24e = (N2_DEn_R2_24e(J))*1E2*c;
        N2_nu_R3_24e = (N2_DEn_R3_24e(J))*1E2*c;
        N2_nu_P1_24e = (N2_DEn_P1_24e(J))*1E2*c;
        N2_nu_P2_24e = (N2_DEn_P2_24e(J))*1E2*c;
        N2_nu_P3_24e = (N2_DEn_P3_24e(J))*1E2*c;
        N2_nu_Q1_24e = (N2_DEn_Q1_24e(J))*1E2*c;
        N2_nu_Q2_24e = (N2_DEn_Q2_24e(J))*1E2*c;
        N2_nu_Q3_24e = (N2_DEn_Q3_24e(J))*1E2*c;
        N2_nu_R1_24f = (N2_DEn_R1_24f(J))*1E2*c;
        N2_nu_R2_24f = (N2_DEn_R2_24f(J))*1E2*c;
        N2_nu_R3_24f = (N2_DEn_R3_24f(J))*1E2*c;
        N2_nu_P1_24f = (N2_DEn_P1_24f(J))*1E2*c;
        N2_nu_P2_24f = (N2_DEn_P2_24f(J))*1E2*c;
        N2_nu_P3_24f = (N2_DEn_P3_24f(J))*1E2*c;
        N2_nu_Q1_24f = (N2_DEn_Q1_24f(J))*1E2*c;
        N2_nu_Q2_24f = (N2_DEn_Q2_24f(J))*1E2*c;
        N2_nu_Q3_24f = (N2_DEn_Q3_24f(J))*1E2*c;  
        if(J>2)
            N2_lambda_Q1_24e = 1E7/N2_DEn_Q1_24e(J);
            N2_lambda_P1_24e = 1E7/N2_DEn_P1_24e(J);
            N2_lambda_R1_24e = 1E7/N2_DEn_R1_24e(J);
            N2_lambda_Q1_24f = 1E7/N2_DEn_Q1_24f(J);
            N2_lambda_P1_24f = 1E7/N2_DEn_P1_24f(J);
            N2_lambda_R1_24f = 1E7/N2_DEn_R1_24f(J);           
        end
        if(J>1)
            N2_lambda_Q2_24e = 1E7/N2_DEn_Q2_24e(J);
            N2_lambda_P2_24e = 1E7/N2_DEn_P2_24e(J);
            N2_lambda_R2_24e = 1E7/N2_DEn_R2_24e(J);
            N2_lambda_Q2_24f = 1E7/N2_DEn_Q2_24f(J);
            N2_lambda_P2_24f = 1E7/N2_DEn_P2_24f(J);
            N2_lambda_R2_24f = 1E7/N2_DEn_R2_24f(J);          
        end
        N2_lambda_Q3_24e = 1E7/N2_DEn_Q3_24e(J);            
        N2_lambda_P3_24e = 1E7/N2_DEn_P3_24e(J);
        N2_lambda_R3_24e = 1E7/N2_DEn_R3_24e(J);
        N2_lambda_Q3_24f = 1E7/N2_DEn_Q3_24f(J);            
        N2_lambda_P3_24f = 1E7/N2_DEn_P3_24f(J);
        N2_lambda_R3_24f = 1E7/N2_DEn_R3_24f(J);
        
            % caso geral N2+ e CN
            SJR1 = SRDoublet1(J,0,0,0)/(2*(2*J+1));
            SJR2 = SRDoublet2(J,0,0,0)/(2*(2*J+1));
            %SJQ1 = SQDoublet1(J,0,0,0)/(2*(2*J+1));
            SJRQ21 = 2*((J+0.5)+1)/(4*(J-0.5)*(J+0.5))/(2*(2*J+1));
            
            % caso geral
            SJR1T = SRTriplet1(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJR2T = SRTriplet2(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJR3T = SRTriplet3(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ1T = SQTriplet1(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ2T = SQTriplet2(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ3T = SQTriplet3(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            
            SJR1T_25 = SRTriplet1(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJR2T_25 = SRTriplet2(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJR3T_25 = SRTriplet3(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ1T_25 = SQTriplet1(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ2T_25 = SQTriplet2(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ3T_25 = SQTriplet3(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));         
            
            SJR1T_47 = SRTriplet1(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJR2T_47 = SRTriplet2(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJR3T_47 = SRTriplet3(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ1T_47 = SQTriplet1(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ2T_47 = SQTriplet2(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ3T_47 = SQTriplet3(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1)); 
            
            SJR1T_14 = SRTriplet1(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJR2T_14 = SRTriplet2(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJR3T_14 = SRTriplet3(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ1T_14 = SQTriplet1(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ2T_14 = SQTriplet2(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ3T_14 = SQTriplet3(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1)); 
            
            SJR1T_03 = SRTriplet1(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJR2T_03 = SRTriplet2(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJR3T_03 = SRTriplet3(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ1T_03 = SQTriplet1(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ2T_03 = SQTriplet2(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ3T_03 = SQTriplet3(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));   
            
            SJR1T_02 = SRTriplet1(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJR2T_02 = SRTriplet2(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJR3T_02 = SRTriplet3(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ1T_02 = SQTriplet1(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ2T_02 = SQTriplet2(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ3T_02 = SQTriplet3(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1)); 
            
            SJR1T_13 = SRTriplet1(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJR2T_13 = SRTriplet2(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJR3T_13 = SRTriplet3(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ1T_13 = SQTriplet1(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ2T_13 = SQTriplet2(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ3T_13 = SQTriplet3(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1)); 
            
            SJR1T_24 = SRTriplet1(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJR2T_24 = SRTriplet2(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJR3T_24 = SRTriplet3(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ1T_24 = SQTriplet1(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ2T_24 = SQTriplet2(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ3T_24 = SQTriplet3(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            
            % (nu 0-0) transition in N2(+) R1, R2 and RQ21 branches
            
            sum1 = sum1 + N2p_ratio*A_00*gR*SJR1*nu_B1*exp(-c*h_p*En_B1(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B1).^2*log(2)/delta_lambda_inst^2);
            sum1 = sum1 + N2p_ratio*A_00*gR*SJR2*nu_B2*exp(-c*h_p*En_B2(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B2).^2*log(2)/delta_lambda_inst^2);
            sum1 = sum1 + N2p_ratio*A_00*gR*SJRQ21*nu_RQ21*exp(-c*h_p*En_B2(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_RQ21).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in N2(+) R1, R2 and RQ21 branches
            
            sum1 = sum1 + N2p_ratio*A_11*gP*SJR1*nu_B1_1*exp(-c*h_p*((En_B1_1(J)-DT_1)*1E2/kb/Trot)).*...
                exp(-(lambda-lambda_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            sum1 = sum1 + N2p_ratio*A_11*gP*SJR2*nu_B2_1*exp(-c*h_p*(En_B2_1(J)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            sum1 = sum1 + N2p_ratio*A_11*gP*SJRQ21*nu_RQ21_1*exp(-c*h_p*(En_B2_1(J)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_RQ21_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            
            % (nu 3-6) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR1T*N2_nu_R1e*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR1T*N2_nu_R1f*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ1T*N2_nu_Q1e*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ1T*N2_nu_Q1f*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));     
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR2T*N2_nu_R2e*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR2T*N2_nu_R2f*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));        
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T*N2_nu_Q2e*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T*N2_nu_Q2f*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            end
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR3T*N2_nu_R3e*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gR_N2*SJR3T*N2_nu_R3f*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ3T*N2_nu_Q3e*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ3T*N2_nu_Q3f*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            % (nu 2-5) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR1T_25*N2_nu_R1_25e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR1T_25*N2_nu_R1_25f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_25*gQ_N2*SJQ1T_25*N2_nu_Q1_25e*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_25*gQ_N2*SJQ1T_25*N2_nu_Q1_25f*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));        
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR2T_25*N2_nu_R2_25e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR2T_25*N2_nu_R2_25f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));            
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T_25*N2_nu_Q2_25e*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T_25*N2_nu_Q2_25f*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR3T_25*N2_nu_R3_25e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_25*gR_N2*SJR3T_25*N2_nu_R3_25f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));         
            sum3 = sum3 + N2_ratio*A_N2_25*gQ_N2*SJQ3T_25*N2_nu_Q3_25e*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_25*gQ_N2*SJQ3T_25*N2_nu_Q3_25f*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            
            % (nu 4-7) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR1T_47*N2_nu_R1_47e*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR1T_47*N2_nu_R1_47f*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ1T_47*N2_nu_Q1_47e*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ1T_47*N2_nu_Q1_47f*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR2T_47*N2_nu_R2_47e*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR2T_47*N2_nu_R2_47f*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ2T_47*N2_nu_Q2_47e*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ2T_47*N2_nu_Q2_47f*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            end
            sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR3T_47*N2_nu_R3_47e*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gR_N2*SJR3T_47*N2_nu_R3_47f*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ3T_47*N2_nu_Q3_47e*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_47*gQ_N2*SJQ3T_47*N2_nu_Q3_47f*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));             
            % (nu 1-4) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR1T_14*N2_nu_R1_14e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR1T_14*N2_nu_R1_14f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ1T_14*N2_nu_Q1_14f*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ1T_14*N2_nu_Q1_14e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR2T_14*N2_nu_R2_14e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR2T_14*N2_nu_R2_14f*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ2T_14*N2_nu_Q2_14e*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ2T_14*N2_nu_Q2_14f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));         
            end
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR3T_14*N2_nu_R3_14e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gR_N2*SJR3T_14*N2_nu_R3_14f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ3T_14*N2_nu_Q3_14f*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gQ_N2*SJQ3T_14*N2_nu_Q3_14e*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            
            % (nu 0-3) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR1T_03*N2_nu_R1_03e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR1T_03*N2_nu_R1_03f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ1T_03*N2_nu_Q1_03f*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ1T_03*N2_nu_Q1_03e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR2T_03*N2_nu_R2_03e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR2T_03*N2_nu_R2_03f*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ2T_03*N2_nu_Q2_03e*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ2T_03*N2_nu_Q2_03f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));         
            end
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR3T_03*N2_nu_R3_03e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gR_N2*SJR3T_03*N2_nu_R3_03f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ3T_03*N2_nu_Q3_03f*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gQ_N2*SJQ3T_03*N2_nu_Q3_03e*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2)); 
            
            % (nu 0-2) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR1T_02*N2_nu_R1_02e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR1T_02*N2_nu_R1_02f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ1T_02*N2_nu_Q1_02f*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ1T_02*N2_nu_Q1_02e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR2T_02*N2_nu_R2_02e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR2T_02*N2_nu_R2_02f*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ2T_02*N2_nu_Q2_02e*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ2T_02*N2_nu_Q2_02f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));         
            end
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR3T_02*N2_nu_R3_02e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gR_N2*SJR3T_02*N2_nu_R3_02f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ3T_02*N2_nu_Q3_02f*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gQ_N2*SJQ3T_02*N2_nu_Q3_02e*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));   
            
            % (nu 1-3) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR1T_13*N2_nu_R1_13e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR1T_13*N2_nu_R1_13f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ1T_13*N2_nu_Q1_13f*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ1T_13*N2_nu_Q1_13e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR2T_13*N2_nu_R2_13e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR2T_13*N2_nu_R2_13f*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ2T_13*N2_nu_Q2_13e*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ2T_13*N2_nu_Q2_13f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));         
            end
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR3T_13*N2_nu_R3_13e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gR_N2*SJR3T_13*N2_nu_R3_13f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ3T_13*N2_nu_Q3_13f*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gQ_N2*SJQ3T_13*N2_nu_Q3_13e*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));  
            
            % (nu 2-4) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR1T_24*N2_nu_R1_24e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR1T_24*N2_nu_R1_24f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_24*gQ_N2*SJQ1T_24*N2_nu_Q1_24e*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_24*gQ_N2*SJQ1T_24*N2_nu_Q1_24f*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));        
            end
            if(J>1)
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR2T_24*N2_nu_R2_24e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR2T_24*N2_nu_R2_24f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));            
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T_24*N2_nu_Q2_24e*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gQ_N2*SJQ2T_24*N2_nu_Q2_24f*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR3T_24*N2_nu_R3_24e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_24*gR_N2*SJR3T_24*N2_nu_R3_24f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));         
            sum3 = sum3 + N2_ratio*A_N2_24*gQ_N2*SJQ3T_24*N2_nu_Q3_24e*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_24*gQ_N2*SJQ3T_24*N2_nu_Q3_24f*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2)); 
            
            % (nu 0-0) transition in CN R1, R2 and RQ21 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJR1*CN_nu_B1*exp(-c*h_p*CN_En_B1(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1).^2*log(2)/delta_lambda_inst^2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJR2*CN_nu_B2*exp(-c*h_p*CN_En_B2(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2).^2*log(2)/delta_lambda_inst^2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJRQ21*CN_nu_RQ21*exp(-c*h_p*CN_En_B2(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in CN R1, R2 and RQ21 branches
        
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJR1*CN_nu_B1_1*exp(-c*h_p*(CN_En_B1_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJR2*CN_nu_B2_1*exp(-c*h_p*(CN_En_B2_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJRQ21*CN_nu_RQ21_1*exp(-c*h_p*(CN_En_B2_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
         
            % (nu 2-2) transition in CN R1, R2 and RQ21 branches
        
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJR1*CN_nu_B1_2*exp(-c*h_p*(CN_En_B1_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJR2*CN_nu_B2_2*exp(-c*h_p*(CN_En_B2_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJRQ21*CN_nu_RQ21_2*exp(-c*h_p*(CN_En_B2_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
         
            % (nu 3-3) transition in CN R1, R2 and RQ21 branches
        
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJR1*CN_nu_B1_3*exp(-c*h_p*(CN_En_B1_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJR2*CN_nu_B2_3*exp(-c*h_p*(CN_En_B2_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJRQ21*CN_nu_RQ21_3*exp(-c*h_p*(CN_En_B2_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);    
          
            % (nu 4-4) transition in CN R1, R2 and RQ21 branches
         
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJR1*CN_nu_B1_4*exp(-c*h_p*(CN_En_B1_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJR2*CN_nu_B2_4*exp(-c*h_p*(CN_En_B2_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJRQ21*CN_nu_RQ21_4*exp(-c*h_p*(CN_En_B2_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);            
            
            % (nu 5-5) transition in CN R1, R2 and RQ21 branches
        
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJR1*CN_nu_B1_5*exp(-c*h_p*(CN_En_B1_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJR2*CN_nu_B2_5*exp(-c*h_p*(CN_En_B2_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJRQ21*CN_nu_RQ21_5*exp(-c*h_p*(CN_En_B2_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
            
            % (nu 6-6) transition in CN R1, R2 and RQ21 branches
        
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJR1*CN_nu_B1_6*exp(-c*h_p*(CN_En_B1_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJR2*CN_nu_B2_6*exp(-c*h_p*(CN_En_B2_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJRQ21*CN_nu_RQ21_6*exp(-c*h_p*(CN_En_B2_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);  
            
            % (nu 7-7) transition in CN R1, R2 and RQ21 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJR1*CN_nu_B1_7*exp(-c*h_p*(CN_En_B1_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJR2*CN_nu_B2_7*exp(-c*h_p*(CN_En_B2_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJRQ21*CN_nu_RQ21_7*exp(-c*h_p*(CN_En_B2_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);
            
            % (nu 8-8) transition in CN R1, R2 and RQ21 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJR1*CN_nu_B1_8*exp(-c*h_p*(CN_En_B1_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJR2*CN_nu_B2_8*exp(-c*h_p*(CN_En_B2_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJRQ21*CN_nu_RQ21_8*exp(-c*h_p*(CN_En_B2_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            
            % (nu 9-9) transition in CN R1, R2 and RQ21 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJR1*CN_nu_B1_9*exp(-c*h_p*(CN_En_B1_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJR2*CN_nu_B2_9*exp(-c*h_p*(CN_En_B2_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJRQ21*CN_nu_RQ21_9*exp(-c*h_p*(CN_En_B2_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            
            % (nu 10-10) transition in CN R1, R2 and RQ21 branches
            
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJR1*CN_nu_B1_10*exp(-c*h_p*(CN_En_B1_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_B1_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJR2*CN_nu_B2_10*exp(-c*h_p*(CN_En_B2_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_B2_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJRQ21*CN_nu_RQ21_10*exp(-c*h_p*(CN_En_B2_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_RQ21_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
            
        if(J>1)
            SJP1 = SPDoublet1(J-1,0,0,0)/(2*(2*(J-1)+1));
            SJP2 = SPDoublet2(J-1,0,0,0)/(2*(2*(J-1)+1));
            %SJQ2 = SQDoublet2(J-1,0,0,0)/(2*(2*(J-1)+1));
            SJPQ12 = (2*(J+0.5)+1)/(4*(J-0.5)*(J+0.5))/(2*(2*(J-1)+1));
            
            SJP1T = SPTriplet1(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            SJP2T = SPTriplet2(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            SJP3T = SPTriplet3(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            
            SJP1T_25 = SPTriplet1(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP2T_25 = SPTriplet2(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP3T_25 = SPTriplet3(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));
            
            SJP1T_47 = SPTriplet1(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));
            SJP2T_47 = SPTriplet2(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));
            SJP3T_47 = SPTriplet3(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));         
            
            SJP1T_14 = SPTriplet1(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP2T_14 = SPTriplet2(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP3T_14 = SPTriplet3(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            
            SJP1T_03 = SPTriplet1(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP2T_03 = SPTriplet2(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP3T_03 = SPTriplet3(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1));
            
            SJP1T_02 = SPTriplet1(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP2T_02 = SPTriplet2(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP3T_02 = SPTriplet3(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            
            SJP1T_13 = SPTriplet1(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP2T_13 = SPTriplet2(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP3T_13 = SPTriplet3(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));

            SJP1T_24 = SPTriplet1(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP2T_24 = SPTriplet2(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP3T_24 = SPTriplet3(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            
            % (nu 0-0) transition in N2(+) P1, P2 and PQ12 branches
            sum1 = sum1 + N2p_ratio*A_00*gP*SJP1*nu_P_B1*exp(-c*h_p*En_B1(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B1).^2*log(2)/delta_lambda_inst^2); 
            sum1 = sum1 + N2p_ratio*A_00*gP*SJP2*nu_P_B2*exp(-c*h_p*En_B2(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B2).^2*log(2)/delta_lambda_inst^2);
             sum1 = sum1 + N2p_ratio*A_00*gP*SJPQ12*nu_PQ12*exp(-c*h_p*En_B1(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_PQ12).^2*log(2)/delta_lambda_inst^2); 
            
            % (nu 1-1) transition in N2(+) P1, P2 and PQ12 branches
            sum1 = sum1 + N2p_ratio*A_11*gR*SJP1*nu_P_B1_1*exp(-c*h_p*(En_B1_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib); 
            sum1 = sum1 + N2p_ratio*A_11*gR*SJP2*nu_P_B2_1*exp(-c*h_p*(En_B2_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib);
            sum1 = sum1 + N2p_ratio*A_11*gR*SJPQ12*nu_PQ12_1*exp(-c*h_p*(En_B1_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_PQ12_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib);  

            % (nu 3-6) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP1T*N2_nu_P1e*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP1T*N2_nu_P1f*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            end
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP2T*N2_nu_P2e*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP2T*N2_nu_P2f*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));            
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP3T*N2_nu_P3e*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2)); 
            sum3 = sum3 + N2_ratio*A_N2_36*gP_N2*SJP3T*N2_nu_P3f*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));          
            
            % (nu 2-5) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP1T_25*N2_nu_P1_25e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP1T_25*N2_nu_P1_25f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            end
            sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP2T_25*N2_nu_P2_25e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP2T_25*N2_nu_P2_25f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP3T_25*N2_nu_P3_25e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));   
            sum3 = sum3 + N2_ratio*A_N2_25*gP_N2*SJP3T_25*N2_nu_P3_25f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));             
            
            % (nu 4-7) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP1T_47*N2_nu_P1_47e*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP1T_47*N2_nu_P1_47f*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP2T_47*N2_nu_P2_47e*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP2T_47*N2_nu_P2_47f*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP3T_47*N2_nu_P3_47e*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2)); 
            sum3 = sum3 + N2_ratio*A_N2_47*gP_N2*SJP3T_47*N2_nu_P3_47f*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));            
            
            % (nu 1-4) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP1T_14*N2_nu_P1_14e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP1T_14*N2_nu_P1_14f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP2T_14*N2_nu_P2_14e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP2T_14*N2_nu_P2_14f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP3T_14*N2_nu_P3_14e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));       
            sum3 = sum3 + N2_ratio*A_N2_14*gP_N2*SJP3T_14*N2_nu_P3_14f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            
            % (nu 0-3) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP1T_03*N2_nu_P1_03e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP1T_03*N2_nu_P1_03f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP2T_03*N2_nu_P2_03e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP2T_03*N2_nu_P2_03f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP3T_03*N2_nu_P3_03e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));       
            sum3 = sum3 + N2_ratio*A_N2_03*gP_N2*SJP3T_03*N2_nu_P3_03f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));         
            
            % (nu 0-2) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP1T_02*N2_nu_P1_02e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP1T_02*N2_nu_P1_02f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP2T_02*N2_nu_P2_02e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP2T_02*N2_nu_P2_02f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP3T_02*N2_nu_P3_02e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));       
            sum3 = sum3 + N2_ratio*A_N2_02*gP_N2*SJP3T_02*N2_nu_P3_02f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));   
            
            % (nu 1-3) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP1T_13*N2_nu_P1_13e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP1T_13*N2_nu_P1_13f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP2T_13*N2_nu_P2_13e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP2T_13*N2_nu_P2_13f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP3T_13*N2_nu_P3_13e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));       
            sum3 = sum3 + N2_ratio*A_N2_13*gP_N2*SJP3T_13*N2_nu_P3_13f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            
            % (nu 2-4) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP1T_24*N2_nu_P1_24e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP1T_24*N2_nu_P1_24f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            end
            sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP2T_24*N2_nu_P2_24e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP2T_24*N2_nu_P2_24f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP3T_24*N2_nu_P3_24e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));   
            sum3 = sum3 + N2_ratio*A_N2_24*gP_N2*SJP3T_24*N2_nu_P3_24f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2)); 
            
            % (nu 0-0) transition in CN, P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJP1*CN_nu_P_B1*exp(-c*h_p*CN_En_B1(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1).^2*log(2)/delta_lambda_inst^2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJP2*CN_nu_P_B2*exp(-c*h_p*CN_En_B2(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2).^2*log(2)/delta_lambda_inst^2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_00*SJPQ12*CN_nu_PQ12*exp(-c*h_p*CN_En_B1(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in CN P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJP1*CN_nu_P_B1_1*exp(-c*h_p*(CN_En_B1_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJP2*CN_nu_P_B2_1*exp(-c*h_p*(CN_En_B2_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
            sum2 = sum2 + gCN*CN_ratio*A_CN_11*SJPQ12*CN_nu_PQ12_1*exp(-c*h_p*(CN_En_B1_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);  
        
            % (nu 2-2) transition in CN P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJP1*CN_nu_P_B1_2*exp(-c*h_p*(CN_En_B1_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJP2*CN_nu_P_B2_2*exp(-c*h_p*(CN_En_B2_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
            sum2 = sum2 + gCN*CN_ratio*A_CN_22*SJPQ12*CN_nu_PQ12_2*exp(-c*h_p*(CN_En_B1_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);  
        
            % (nu 3-3) transition in CN P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJP1*CN_nu_P_B1_3*exp(-c*h_p*(CN_En_B1_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJP2*CN_nu_P_B2_3*exp(-c*h_p*(CN_En_B2_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);
            sum2 = sum2 + gCN*CN_ratio*A_CN_33*SJPQ12*CN_nu_PQ12_3*exp(-c*h_p*(CN_En_B1_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);   
        
            % (nu 4-4) transition in CN P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJP1*CN_nu_P_B1_4*exp(-c*h_p*(CN_En_B1_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJP2*CN_nu_P_B2_4*exp(-c*h_p*(CN_En_B2_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);
            sum2 = sum2 + gCN*CN_ratio*A_CN_44*SJPQ12*CN_nu_PQ12_4*exp(-c*h_p*(CN_En_B1_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);   
        
            % (nu 5-5) transition in CN P1, P2 and PQ12 branches
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJP1*CN_nu_P_B1_5*exp(-c*h_p*(CN_En_B1_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJP2*CN_nu_P_B2_5*exp(-c*h_p*(CN_En_B2_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_55*SJPQ12*CN_nu_PQ12_5*exp(-c*h_p*(CN_En_B1_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
        
            % (nu 6-6) transition in CN P1, P2 and PQ12 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJP1*CN_nu_P_B1_6*exp(-c*h_p*(CN_En_B1_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJP2*CN_nu_P_B2_6*exp(-c*h_p*(CN_En_B2_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_66*SJPQ12*CN_nu_PQ12_6*exp(-c*h_p*(CN_En_B1_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6); 
        
            % (nu 7-7) transition in CN P1, P2 and PQ12 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJP1*CN_nu_P_B1_7*exp(-c*h_p*(CN_En_B1_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJP2*CN_nu_P_B2_7*exp(-c*h_p*(CN_En_B2_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_77*SJPQ12*CN_nu_PQ12_7*exp(-c*h_p*(CN_En_B1_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);   
        
            % (nu 8-8) transition in CN P1, P2 and PQ12 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJP1*CN_nu_P_B1_8*exp(-c*h_p*(CN_En_B1_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);  
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJP2*CN_nu_P_B2_8*exp(-c*h_p*(CN_En_B2_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            sum2 = sum2 + gCN*CN_ratio*A_CN_88*SJPQ12*CN_nu_PQ12_8*exp(-c*h_p*(CN_En_B1_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
        
            % (nu 9-9) transition in CN P1, P2 and PQ12 branches
            
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJP1*CN_nu_P_B1_9*exp(-c*h_p*(CN_En_B1_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9); 
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJP2*CN_nu_P_B2_9*exp(-c*h_p*(CN_En_B2_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            sum2 = sum2 + gCN*CN_ratio*A_CN_99*SJPQ12*CN_nu_PQ12_9*exp(-c*h_p*(CN_En_B1_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
        
            % (nu 10-10) transition in CN P1, P2 and PQ12 branches
            
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJP1*CN_nu_P_B1_10*exp(-c*h_p*(CN_En_B1_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_P_B1_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);  
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJP2*CN_nu_P_B2_10*exp(-c*h_p*(CN_En_B2_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_P_B2_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN); 
%             sum2 = sum2 + gCN*CN_ratio*A_CN_1010*SJPQ12*CN_nu_PQ12_10*exp(-c*h_p*(CN_En_B1_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_PQ12_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
            
        end
    end
    % add artificial gaussian shape profile to model the unknown line 
%    sum4 = sum4 + beta(9)*exp(-(lambda-395.028).^2*log(2)/delta_lambda_inst^2);
%    sum4 = sum4 + beta(10)*exp(-(lambda-394.89).^2*log(2)/delta_lambda_inst^2);
    %    beta(6)*exp(-(lambda-383.3).^2*log(2)/beta(7)^2);
    I1 = sum1/max(sum2);
    I2 = sum2/max(sum2);
    I3 = sum3/max(sum2);
    I4 = 0;
end

function I = I_syntetic(beta, lambda)
global En_B1 En_B2 DEn_R1 DEn_R2 DEn_P1 DEn_P2 DEn_RQ21 DEn_PQ12;
global En_B1_1 En_B2_1 DEn_R1_1 DEn_R2_1 DEn_P1_1 DEn_P2_1 DEn_RQ21_1 DEn_PQ12_1 DT_1;
global CN_En_B1 CN_En_B2 CN_DEn_R1 CN_DEn_R2 CN_DEn_P1 CN_DEn_P2 CN_DEn_RQ21 CN_DEn_PQ12;
global CN_En_B1_1 CN_En_B2_1 CN_DEn_R1_1 CN_DEn_R2_1 CN_DEn_P1_1 CN_DEn_P2_1 CN_DEn_RQ21_1 CN_DEn_PQ12_1;
global CN_En_B1_2 CN_En_B2_2 CN_DEn_R1_2 CN_DEn_R2_2 CN_DEn_P1_2 CN_DEn_P2_2 CN_DEn_RQ21_2 CN_DEn_PQ12_2;
global CN_En_B1_3 CN_En_B2_3 CN_DEn_R1_3 CN_DEn_R2_3 CN_DEn_P1_3 CN_DEn_P2_3 CN_DEn_RQ21_3 CN_DEn_PQ12_3;
global CN_En_B1_4 CN_En_B2_4 CN_DEn_R1_4 CN_DEn_R2_4 CN_DEn_P1_4 CN_DEn_P2_4 CN_DEn_RQ21_4 CN_DEn_PQ12_4;
global CN_En_B1_5 CN_En_B2_5 CN_DEn_R1_5 CN_DEn_R2_5 CN_DEn_P1_5 CN_DEn_P2_5 CN_DEn_RQ21_5 CN_DEn_PQ12_5;
global CN_En_B1_6 CN_En_B2_6 CN_DEn_R1_6 CN_DEn_R2_6 CN_DEn_P1_6 CN_DEn_P2_6 CN_DEn_RQ21_6 CN_DEn_PQ12_6;
global CN_En_B1_7 CN_En_B2_7 CN_DEn_R1_7 CN_DEn_R2_7 CN_DEn_P1_7 CN_DEn_P2_7 CN_DEn_RQ21_7 CN_DEn_PQ12_7;
global CN_En_B1_8 CN_En_B2_8 CN_DEn_R1_8 CN_DEn_R2_8 CN_DEn_P1_8 CN_DEn_P2_8 CN_DEn_RQ21_8 CN_DEn_PQ12_8;
global CN_En_B1_9 CN_En_B2_9 CN_DEn_R1_9 CN_DEn_R2_9 CN_DEn_P1_9 CN_DEn_P2_9 CN_DEn_RQ21_9 CN_DEn_PQ12_9;
global CN_En_B1_10 CN_En_B2_10 CN_DEn_R1_10 CN_DEn_R2_10 CN_DEn_P1_10 CN_DEn_P2_10 CN_DEn_RQ21_10 CN_DEn_PQ12_10;
global CN_DT1 CN_DT2 CN_DT3 CN_DT4 CN_DT5 CN_DT6 CN_DT7 CN_DT8 CN_DT9 CN_DT10;  
global delta_lambda_inst;
global N_J;
global h_p c kb;
global N2_Y_B N2_Y_B_4 N2_Y_B_5 N2_Y_B_7 N2_Y_C N2_Y_C_2 N2_Y_C_4 N2_Y_C_1 ...
    N2_Y_B_3 N2_Y_B_2 N2_Y_C_0;
global N2_En_C1_3e N2_En_C2_3e N2_En_C3_3e N2_En_C1_3f N2_En_C2_3f N2_En_C3_3f... 
    N2_DEn_R1_36e N2_DEn_R2_36e N2_DEn_R3_36e N2_DEn_P1_36e N2_DEn_P2_36e N2_DEn_P3_36e ...
    N2_DEn_Q1_36e N2_DEn_Q2_36e N2_DEn_Q3_36e ...
    N2_DEn_R1_36f N2_DEn_R2_36f N2_DEn_R3_36f N2_DEn_P1_36f N2_DEn_P2_36f N2_DEn_P3_36f ...
    N2_DEn_Q1_36f N2_DEn_Q2_36f N2_DEn_Q3_36f;
global N2_En_C1_2e N2_En_C2_2e N2_En_C3_2e N2_En_C1_2f N2_En_C2_2f N2_En_C3_2f ...
    N2_DEn_R1_25e N2_DEn_R2_25e N2_DEn_R3_25e N2_DEn_P1_25e N2_DEn_P2_25e N2_DEn_P3_25e ...
    N2_DEn_Q1_25e N2_DEn_Q2_25e N2_DEn_Q3_25e ...
    N2_DEn_R1_25f N2_DEn_R2_25f N2_DEn_R3_25f N2_DEn_P1_25f N2_DEn_P2_25f N2_DEn_P3_25f ...
    N2_DEn_Q1_25f N2_DEn_Q2_25f N2_DEn_Q3_25f;
global N2_En_C1_4e N2_En_C2_4e N2_En_C3_4e N2_En_C1_4f N2_En_C2_4f N2_En_C3_4f ...
    N2_DEn_R1_47e N2_DEn_R2_47e N2_DEn_R3_47e N2_DEn_P1_47e N2_DEn_P2_47e N2_DEn_P3_47e ...
    N2_DEn_Q1_47e N2_DEn_Q2_47e N2_DEn_Q3_47e ...
    N2_DEn_R1_47f N2_DEn_R2_47f N2_DEn_R3_47f N2_DEn_P1_47f N2_DEn_P2_47f N2_DEn_P3_47f ...
    N2_DEn_Q1_47f N2_DEn_Q2_47f N2_DEn_Q3_47f;
global N2_En_C1_1e N2_En_C2_1e N2_En_C3_1e N2_En_C1_1f N2_En_C2_1f N2_En_C3_1f ... 
    N2_DEn_R1_14e N2_DEn_R2_14e N2_DEn_R3_14e N2_DEn_P1_14e N2_DEn_P2_14e N2_DEn_P3_14e ...
    N2_DEn_Q1_14e N2_DEn_Q2_14e N2_DEn_Q3_14e ...
    N2_DEn_R1_14f N2_DEn_R2_14f N2_DEn_R3_14f N2_DEn_P1_14f N2_DEn_P2_14f N2_DEn_P3_14f ...
    N2_DEn_Q1_14f N2_DEn_Q2_14f N2_DEn_Q3_14f;
global N2_En_C1_0e N2_En_C2_0e N2_En_C3_0e N2_En_C1_0f N2_En_C2_0f N2_En_C3_0f ... 
    N2_DEn_R1_03e N2_DEn_R2_03e N2_DEn_R3_03e N2_DEn_P1_03e N2_DEn_P2_03e N2_DEn_P3_03e ...
    N2_DEn_Q1_03e N2_DEn_Q2_03e N2_DEn_Q3_03e ...
    N2_DEn_R1_03f N2_DEn_R2_03f N2_DEn_R3_03f N2_DEn_P1_03f N2_DEn_P2_03f N2_DEn_P3_03f ...
    N2_DEn_Q1_03f N2_DEn_Q2_03f N2_DEn_Q3_03f;
global N2_DEn_R1_02e N2_DEn_R2_02e N2_DEn_R3_02e N2_DEn_P1_02e N2_DEn_P2_02e N2_DEn_P3_02e ...
    N2_DEn_Q1_02e N2_DEn_Q2_02e N2_DEn_Q3_02e ...
    N2_DEn_R1_02f N2_DEn_R2_02f N2_DEn_R3_02f N2_DEn_P1_02f N2_DEn_P2_02f N2_DEn_P3_02f ...
    N2_DEn_Q1_02f N2_DEn_Q2_02f N2_DEn_Q3_02f;
global  N2_DEn_R1_13e N2_DEn_R2_13e N2_DEn_R3_13e N2_DEn_P1_13e N2_DEn_P2_13e N2_DEn_P3_13e ...
    N2_DEn_Q1_13e N2_DEn_Q2_13e N2_DEn_Q3_13e ...
    N2_DEn_R1_13f N2_DEn_R2_13f N2_DEn_R3_13f N2_DEn_P1_13f N2_DEn_P2_13f N2_DEn_P3_13f ...
    N2_DEn_Q1_13f N2_DEn_Q2_13f N2_DEn_Q3_13f;
global  N2_DEn_R1_24e N2_DEn_R2_24e N2_DEn_R3_24e N2_DEn_P1_24e N2_DEn_P2_24e N2_DEn_P3_24e ...
    N2_DEn_Q1_24e N2_DEn_Q2_24e N2_DEn_Q3_24e ...
    N2_DEn_R1_24f N2_DEn_R2_24f N2_DEn_R3_24f N2_DEn_P1_24f N2_DEn_P2_24f N2_DEn_P3_24f ...
    N2_DEn_Q1_24f N2_DEn_Q2_24f N2_DEn_Q3_24f;
global V_C_3 V_C_2 V_C_4 V_C_1 V_C_0;
global A_00 A_11;
global A_CN_00 A_CN_11 A_CN_22 A_CN_33 A_CN_44 A_CN_55 A_CN_66 ...
    A_CN_77 A_CN_88 A_CN_99 A_CN_1010;
global A_N2_36 A_N2_25 A_N2_47 A_N2_14 A_N2_03 A_N2_02 A_N2_13 A_N2_24; 

lambda = lambda(:);
Trot = beta(1);
Tvib = beta(1);
TrotCN = beta(1);
TrotN2 = beta(1);
TvibCN = beta(1);
TvibN2 = beta(1);
CN_ratio = 1;
N2p_ratio = beta(2);
N2_ratio = beta(4);
gamma = beta(5);
delta = beta(6);
eta = 0.0;
sum = zeros(length(lambda),1);
delta_lambda_inst = beta(3);
f_DT1 = -c*h_p*CN_DT1*1E2*(1-gamma + eta*CN_DT1)/kb/TvibCN+delta;
f_DT2 = -c*h_p*CN_DT2*1E2*(1-gamma + eta*CN_DT2)/kb/TvibCN+delta;
f_DT3 = -c*h_p*CN_DT3*1E2*(1-gamma + eta*CN_DT3)/kb/TvibCN+delta;
f_DT4 = -c*h_p*CN_DT4*1E2*(1-gamma + eta*CN_DT4)/kb/TvibCN+delta;
f_DT5 = -c*h_p*CN_DT5*1E2*(1-gamma + eta*CN_DT5)/kb/TvibCN+delta;
f_DT6 = -c*h_p*CN_DT6*1E2*(1-gamma + eta*CN_DT6)/kb/TvibCN+delta;
f_DT7 = -c*h_p*CN_DT7*1E2*(1-gamma + eta*CN_DT7)/kb/TvibCN+delta;
f_DT8 = -c*h_p*CN_DT8*1E2*(1-gamma + eta*CN_DT8)/kb/TvibCN+delta;
f_DT9 = -c*h_p*CN_DT9*1E2*(1-gamma + eta*CN_DT9)/kb/TvibCN+delta;
    for J=1:N_J
        if(mod(J,2)==1)
            gP = 2*(2/3)*(2*J+1);
            gR = 2*(1/3)*(2*(J-1)+1);
            gP_N2 = 3*(2/3)*(2*J+1);
            gR_N2 = 3*(1/3)*(2*(J-1)+1);
            gQ_N2 = 3*(2/3)*(2*J+1);
        else
            gP = 2*(1/3)*(2*J+1);
            gR = 2*(2/3)*(2*(J-1)+1);
            gP_N2 = 3*(1/3)*(2*J+1);
            gR_N2 = 3*(2/3)*(2*(J-1)+1);
            gQ_N2 = 3*(1/3)*(2*J+1);
        end
        gCN = 2*(2*J+1);
        
        % N2(+) frequencies and wavelengths
        nu_B1 = (DEn_R1(J))*1E2*c;
        nu_B2 = (DEn_R2(J))*1E2*c;
        nu_RQ21 = (DEn_RQ21(J))*1E2*c;
        nu_B1_1 = (DEn_R1_1(J))*1E2*c;
        nu_B2_1 = (DEn_R2_1(J))*1E2*c;
        nu_RQ21_1 = (DEn_RQ21_1(J))*1E2*c;
        lambda_B1 = 1E7/DEn_R1(J);
        lambda_B2 = 1E7/DEn_R2(J);
        lambda_RQ21 = 1E7/DEn_RQ21(J);
        lambda_B1_1 = 1E7/DEn_R1_1(J);
        lambda_B2_1 = 1E7/DEn_R2_1(J);
        lambda_RQ21_1 = 1E7/DEn_RQ21_1(J);
        
        % CN frequencies and wavelengths
        CN_nu_B1 = (CN_DEn_R1(J))*1E2*c;
        CN_nu_B2 = (CN_DEn_R2(J))*1E2*c;
        CN_nu_RQ21 = (CN_DEn_RQ21(J))*1E2*c;
        CN_nu_B1_1 = (CN_DEn_R1_1(J))*1E2*c;
        CN_nu_B1_2 = (CN_DEn_R1_2(J))*1E2*c;
        CN_nu_B1_3 = (CN_DEn_R1_3(J))*1E2*c;
        CN_nu_B1_4 = (CN_DEn_R1_4(J))*1E2*c;
        CN_nu_B1_5 = (CN_DEn_R1_5(J))*1E2*c;
        CN_nu_B1_6 = (CN_DEn_R1_6(J))*1E2*c;
        CN_nu_B1_7 = (CN_DEn_R1_7(J))*1E2*c;
        CN_nu_B1_8 = (CN_DEn_R1_8(J))*1E2*c;
        CN_nu_B1_9 = (CN_DEn_R1_9(J))*1E2*c;
        CN_nu_B1_10 = (CN_DEn_R1_10(J))*1E2*c;
        CN_nu_B2_1 = (CN_DEn_R2_1(J))*1E2*c;
        CN_nu_B2_2 = (CN_DEn_R2_2(J))*1E2*c;
        CN_nu_B2_3 = (CN_DEn_R2_3(J))*1E2*c;
        CN_nu_B2_4 = (CN_DEn_R2_4(J))*1E2*c;
        CN_nu_B2_5 = (CN_DEn_R2_5(J))*1E2*c;
        CN_nu_B2_6 = (CN_DEn_R2_6(J))*1E2*c;
        CN_nu_B2_7 = (CN_DEn_R2_7(J))*1E2*c;
        CN_nu_B2_8 = (CN_DEn_R2_8(J))*1E2*c;
        CN_nu_B2_9 = (CN_DEn_R2_9(J))*1E2*c;
        CN_nu_B2_10 = (CN_DEn_R2_10(J))*1E2*c;
        CN_nu_RQ21_1 = (CN_DEn_RQ21_1(J))*1E2*c;
        CN_nu_RQ21_2 = (CN_DEn_RQ21_2(J))*1E2*c;
        CN_nu_RQ21_3 = (CN_DEn_RQ21_3(J))*1E2*c;
        CN_nu_RQ21_4 = (CN_DEn_RQ21_4(J))*1E2*c;
        CN_nu_RQ21_5 = (CN_DEn_RQ21_5(J))*1E2*c;
        CN_nu_RQ21_6 = (CN_DEn_RQ21_6(J))*1E2*c;
        CN_nu_RQ21_7 = (CN_DEn_RQ21_7(J))*1E2*c;
        CN_nu_RQ21_8 = (CN_DEn_RQ21_8(J))*1E2*c;
        CN_nu_RQ21_9 = (CN_DEn_RQ21_9(J))*1E2*c;
        CN_nu_RQ21_10 = (CN_DEn_RQ21_10(J))*1E2*c;
        CN_lambda_B1 = 1E7/CN_DEn_R1(J);
        CN_lambda_B2 = 1E7/CN_DEn_R2(J);
        CN_lambda_RQ21 = 1E7/CN_DEn_RQ21(J);
        CN_lambda_B1_1 = 1E7/CN_DEn_R1_1(J);
        CN_lambda_B1_2 = 1E7/CN_DEn_R1_2(J);
        CN_lambda_B1_3 = 1E7/CN_DEn_R1_3(J);
        CN_lambda_B1_4 = 1E7/CN_DEn_R1_4(J);
        CN_lambda_B1_5 = 1E7/CN_DEn_R1_5(J);
        CN_lambda_B1_6 = 1E7/CN_DEn_R1_6(J);
        CN_lambda_B1_7 = 1E7/CN_DEn_R1_7(J);
        CN_lambda_B1_8 = 1E7/CN_DEn_R1_8(J);
        CN_lambda_B1_9 = 1E7/CN_DEn_R1_9(J);
        CN_lambda_B1_10 = 1E7/CN_DEn_R1_10(J);
        CN_lambda_B2_1 = 1E7/CN_DEn_R2_1(J);
        CN_lambda_B2_2 = 1E7/CN_DEn_R2_2(J);
        CN_lambda_B2_3 = 1E7/CN_DEn_R2_3(J);
        CN_lambda_B2_4 = 1E7/CN_DEn_R2_4(J);
        CN_lambda_B2_5 = 1E7/CN_DEn_R2_5(J);
        CN_lambda_B2_6 = 1E7/CN_DEn_R2_6(J);
        CN_lambda_B2_7 = 1E7/CN_DEn_R2_7(J);
        CN_lambda_B2_8 = 1E7/CN_DEn_R2_8(J);
        CN_lambda_B2_9 = 1E7/CN_DEn_R2_9(J);
        CN_lambda_B2_10 = 1E7/CN_DEn_R2_10(J);
        CN_lambda_RQ21_1 = 1E7/CN_DEn_RQ21_1(J);
        CN_lambda_RQ21_2 = 1E7/CN_DEn_RQ21_2(J);
        CN_lambda_RQ21_3 = 1E7/CN_DEn_RQ21_3(J);
        CN_lambda_RQ21_4 = 1E7/CN_DEn_RQ21_4(J);
        CN_lambda_RQ21_5 = 1E7/CN_DEn_RQ21_5(J);
        CN_lambda_RQ21_6 = 1E7/CN_DEn_RQ21_6(J);
        CN_lambda_RQ21_7 = 1E7/CN_DEn_RQ21_7(J);
        CN_lambda_RQ21_8 = 1E7/CN_DEn_RQ21_8(J);
        CN_lambda_RQ21_9 = 1E7/CN_DEn_RQ21_9(J);
        CN_lambda_RQ21_10 = 1E7/CN_DEn_RQ21_10(J);
        
            
       % N2(+) frequencies and wavelengths
       nu_P_B1 = (DEn_P1(J))*1E2*c;
       nu_P_B2 = (DEn_P2(J))*1E2*c;
       nu_PQ12 = (DEn_PQ12(J))*1E2*c;
       nu_P_B1_1 = (DEn_P1_1(J))*1E2*c;
       nu_P_B2_1 = (DEn_P2_1(J))*1E2*c;
       nu_PQ12_1 = (DEn_PQ12_1(J))*1E2*c;
       lambda_P_B1 = 1E7/DEn_P1(J);
       lambda_P_B2 = 1E7/DEn_P2(J);
       lambda_PQ12 = 1E7/DEn_PQ12(J);
       lambda_P_B1_1 = 1E7/DEn_P1_1(J);
       lambda_P_B2_1 = 1E7/DEn_P2_1(J);
       lambda_PQ12_1 = 1E7/DEn_PQ12_1(J);
            
        % CN frequencies and wavelengths
        CN_nu_P_B1 = (CN_DEn_P1(J))*1E2*c;
        CN_nu_P_B2 = (CN_DEn_P2(J))*1E2*c;
        CN_nu_PQ12 = (CN_DEn_PQ12(J))*1E2*c;
        CN_nu_P_B1_1 = (CN_DEn_P1_1(J))*1E2*c;
        CN_nu_P_B1_2 = (CN_DEn_P1_2(J))*1E2*c;
        CN_nu_P_B1_3 = (CN_DEn_P1_3(J))*1E2*c;
        CN_nu_P_B1_4 = (CN_DEn_P1_4(J))*1E2*c;
        CN_nu_P_B1_5 = (CN_DEn_P1_5(J))*1E2*c;
        CN_nu_P_B1_6 = (CN_DEn_P1_6(J))*1E2*c;
        CN_nu_P_B1_7 = (CN_DEn_P1_7(J))*1E2*c;
        CN_nu_P_B1_8 = (CN_DEn_P1_8(J))*1E2*c;
        CN_nu_P_B1_9 = (CN_DEn_P1_9(J))*1E2*c;
        CN_nu_P_B1_10 = (CN_DEn_P1_10(J))*1E2*c;
        CN_nu_P_B2_1 = (CN_DEn_P2_1(J))*1E2*c;
        CN_nu_P_B2_2 = (CN_DEn_P2_2(J))*1E2*c;
        CN_nu_P_B2_3 = (CN_DEn_P2_3(J))*1E2*c;
        CN_nu_P_B2_4 = (CN_DEn_P2_4(J))*1E2*c;
        CN_nu_P_B2_5 = (CN_DEn_P2_5(J))*1E2*c;
        CN_nu_P_B2_6 = (CN_DEn_P2_6(J))*1E2*c;
        CN_nu_P_B2_7 = (CN_DEn_P2_7(J))*1E2*c;
        CN_nu_P_B2_8 = (CN_DEn_P2_8(J))*1E2*c;
        CN_nu_P_B2_9 = (CN_DEn_P2_9(J))*1E2*c;
        CN_nu_P_B2_10 = (CN_DEn_P2_10(J))*1E2*c;
        CN_nu_PQ12_1 = (CN_DEn_PQ12_1(J))*1E2*c;
        CN_nu_PQ12_2 = (CN_DEn_PQ12_2(J))*1E2*c;
        CN_nu_PQ12_3 = (CN_DEn_PQ12_3(J))*1E2*c;
        CN_nu_PQ12_4 = (CN_DEn_PQ12_4(J))*1E2*c;
        CN_nu_PQ12_5 = (CN_DEn_PQ12_5(J))*1E2*c;
        CN_nu_PQ12_6 = (CN_DEn_PQ12_6(J))*1E2*c;
        CN_nu_PQ12_7 = (CN_DEn_PQ12_7(J))*1E2*c;
        CN_nu_PQ12_8 = (CN_DEn_PQ12_9(J))*1E2*c;
        CN_nu_PQ12_9 = (CN_DEn_PQ12_10(J))*1E2*c;
        CN_nu_PQ12_10 = (CN_DEn_PQ12_8(J))*1E2*c;
        CN_lambda_P_B1 = 1E7/CN_DEn_P1(J);
        CN_lambda_P_B2 = 1E7/CN_DEn_P2(J);
        CN_lambda_PQ12 = 1E7/CN_DEn_PQ12(J);
        CN_lambda_P_B1_1 = 1E7/CN_DEn_P1_1(J);
        CN_lambda_P_B1_2 = 1E7/CN_DEn_P1_2(J);
        CN_lambda_P_B1_3 = 1E7/CN_DEn_P1_3(J);
        CN_lambda_P_B1_4 = 1E7/CN_DEn_P1_4(J);
        CN_lambda_P_B1_5 = 1E7/CN_DEn_P1_5(J);
        CN_lambda_P_B1_6 = 1E7/CN_DEn_P1_6(J);
        CN_lambda_P_B1_7 = 1E7/CN_DEn_P1_7(J);
        CN_lambda_P_B1_8 = 1E7/CN_DEn_P1_8(J);
        CN_lambda_P_B1_9 = 1E7/CN_DEn_P1_9(J);
        CN_lambda_P_B1_10 = 1E7/CN_DEn_P1_10(J);
        CN_lambda_P_B2_1 = 1E7/CN_DEn_P2_1(J);
        CN_lambda_P_B2_2 = 1E7/CN_DEn_P2_2(J);
        CN_lambda_P_B2_3 = 1E7/CN_DEn_P2_3(J);
        CN_lambda_P_B2_4 = 1E7/CN_DEn_P2_4(J);
        CN_lambda_P_B2_5 = 1E7/CN_DEn_P2_5(J);
        CN_lambda_P_B2_6 = 1E7/CN_DEn_P2_6(J);
        CN_lambda_P_B2_7 = 1E7/CN_DEn_P2_7(J);
        CN_lambda_P_B2_8 = 1E7/CN_DEn_P2_8(J);
        CN_lambda_P_B2_9 = 1E7/CN_DEn_P2_9(J);
        CN_lambda_P_B2_10 = 1E7/CN_DEn_P2_10(J);
        CN_lambda_PQ12_1 = 1E7/CN_DEn_PQ12_1(J);
        CN_lambda_PQ12_2 = 1E7/CN_DEn_PQ12_2(J);  
        CN_lambda_PQ12_3 = 1E7/CN_DEn_PQ12_3(J);
        CN_lambda_PQ12_4 = 1E7/CN_DEn_PQ12_4(J);
        CN_lambda_PQ12_5 = 1E7/CN_DEn_PQ12_5(J);
        CN_lambda_PQ12_6 = 1E7/CN_DEn_PQ12_6(J);
        CN_lambda_PQ12_7 = 1E7/CN_DEn_PQ12_7(J); 
        CN_lambda_PQ12_8 = 1E7/CN_DEn_PQ12_8(J);
        CN_lambda_PQ12_9 = 1E7/CN_DEn_PQ12_9(J); 
        CN_lambda_PQ12_10 = 1E7/CN_DEn_PQ12_10(J);  
        
        % N2 frequencies and wavelengths
        % (3-6)
        N2_nu_R1e = (N2_DEn_R1_36e(J))*1E2*c;
        N2_nu_R2e = (N2_DEn_R2_36e(J))*1E2*c;
        N2_nu_R3e = (N2_DEn_R3_36e(J))*1E2*c;
        N2_nu_P1e = (N2_DEn_P1_36e(J))*1E2*c;
        N2_nu_P2e = (N2_DEn_P2_36e(J))*1E2*c;
        N2_nu_P3e = (N2_DEn_P3_36e(J))*1E2*c;
        N2_nu_Q1e = (N2_DEn_Q1_36e(J))*1E2*c;
        N2_nu_Q2e = (N2_DEn_Q2_36e(J))*1E2*c;
        N2_nu_Q3e = (N2_DEn_Q3_36e(J))*1E2*c;
        N2_nu_R1f = (N2_DEn_R1_36f(J))*1E2*c;
        N2_nu_R2f = (N2_DEn_R2_36f(J))*1E2*c;
        N2_nu_R3f = (N2_DEn_R3_36f(J))*1E2*c;
        N2_nu_P1f = (N2_DEn_P1_36f(J))*1E2*c;
        N2_nu_P2f = (N2_DEn_P2_36f(J))*1E2*c;
        N2_nu_P3f = (N2_DEn_P3_36f(J))*1E2*c;
        N2_nu_Q1f = (N2_DEn_Q1_36f(J))*1E2*c;
        N2_nu_Q2f = (N2_DEn_Q2_36f(J))*1E2*c;
        N2_nu_Q3f = (N2_DEn_Q3_36f(J))*1E2*c;     
        if(J>2)
            N2_lambda_Q1e = 1E7/N2_DEn_Q1_36e(J);
            N2_lambda_P1e = 1E7/N2_DEn_P1_36e(J);
            N2_lambda_R1e = 1E7/N2_DEn_R1_36e(J);
            N2_lambda_Q1f = 1E7/N2_DEn_Q1_36f(J);
            N2_lambda_P1f = 1E7/N2_DEn_P1_36f(J);
            N2_lambda_R1f = 1E7/N2_DEn_R1_36f(J);         
        end
        if(J>1)
            N2_lambda_Q2e = 1E7/N2_DEn_Q2_36e(J);
            N2_lambda_P2e = 1E7/N2_DEn_P2_36e(J);
            N2_lambda_R2e = 1E7/N2_DEn_R2_36e(J);
            N2_lambda_Q2f = 1E7/N2_DEn_Q2_36f(J);
            N2_lambda_P2f = 1E7/N2_DEn_P2_36f(J);
            N2_lambda_R2f = 1E7/N2_DEn_R2_36f(J);          
        end
        N2_lambda_Q3e = 1E7/N2_DEn_Q3_36e(J);            
        N2_lambda_P3e = 1E7/N2_DEn_P3_36e(J);
        N2_lambda_R3e = 1E7/N2_DEn_R3_36e(J);
        N2_lambda_Q3f = 1E7/N2_DEn_Q3_36f(J);            
        N2_lambda_P3f = 1E7/N2_DEn_P3_36f(J);
        N2_lambda_R3f = 1E7/N2_DEn_R3_36f(J);     
        
        % (2-5)
        N2_nu_R1_25e = (N2_DEn_R1_25e(J))*1E2*c;
        N2_nu_R2_25e = (N2_DEn_R2_25e(J))*1E2*c;
        N2_nu_R3_25e = (N2_DEn_R3_25e(J))*1E2*c;
        N2_nu_P1_25e = (N2_DEn_P1_25e(J))*1E2*c;
        N2_nu_P2_25e = (N2_DEn_P2_25e(J))*1E2*c;
        N2_nu_P3_25e = (N2_DEn_P3_25e(J))*1E2*c;
        N2_nu_Q1_25e = (N2_DEn_Q1_25e(J))*1E2*c;
        N2_nu_Q2_25e = (N2_DEn_Q2_25e(J))*1E2*c;
        N2_nu_Q3_25e = (N2_DEn_Q3_25e(J))*1E2*c;
        N2_nu_R1_25f = (N2_DEn_R1_25f(J))*1E2*c;
        N2_nu_R2_25f = (N2_DEn_R2_25f(J))*1E2*c;
        N2_nu_R3_25f = (N2_DEn_R3_25f(J))*1E2*c;
        N2_nu_P1_25f = (N2_DEn_P1_25f(J))*1E2*c;
        N2_nu_P2_25f = (N2_DEn_P2_25f(J))*1E2*c;
        N2_nu_P3_25f = (N2_DEn_P3_25f(J))*1E2*c;
        N2_nu_Q1_25f = (N2_DEn_Q1_25f(J))*1E2*c;
        N2_nu_Q2_25f = (N2_DEn_Q2_25f(J))*1E2*c;
        N2_nu_Q3_25f = (N2_DEn_Q3_25f(J))*1E2*c;  
        if(J>2)
            N2_lambda_Q1_25e = 1E7/N2_DEn_Q1_25e(J);
            N2_lambda_P1_25e = 1E7/N2_DEn_P1_25e(J);
            N2_lambda_R1_25e = 1E7/N2_DEn_R1_25e(J);
            N2_lambda_Q1_25f = 1E7/N2_DEn_Q1_25f(J);
            N2_lambda_P1_25f = 1E7/N2_DEn_P1_25f(J);
            N2_lambda_R1_25f = 1E7/N2_DEn_R1_25f(J);           
        end
        if(J>1)
            N2_lambda_Q2_25e = 1E7/N2_DEn_Q2_25e(J);
            N2_lambda_P2_25e = 1E7/N2_DEn_P2_25e(J);
            N2_lambda_R2_25e = 1E7/N2_DEn_R2_25e(J);
            N2_lambda_Q2_25f = 1E7/N2_DEn_Q2_25f(J);
            N2_lambda_P2_25f = 1E7/N2_DEn_P2_25f(J);
            N2_lambda_R2_25f = 1E7/N2_DEn_R2_25f(J);          
        end
        N2_lambda_Q3_25e = 1E7/N2_DEn_Q3_25e(J);            
        N2_lambda_P3_25e = 1E7/N2_DEn_P3_25e(J);
        N2_lambda_R3_25e = 1E7/N2_DEn_R3_25e(J);
        N2_lambda_Q3_25f = 1E7/N2_DEn_Q3_25f(J);            
        N2_lambda_P3_25f = 1E7/N2_DEn_P3_25f(J);
        N2_lambda_R3_25f = 1E7/N2_DEn_R3_25f(J);      
        
        % (4-7)
        N2_nu_R1_47e = (N2_DEn_R1_47e(J))*1E2*c;
        N2_nu_R2_47e = (N2_DEn_R2_47e(J))*1E2*c;
        N2_nu_R3_47e = (N2_DEn_R3_47e(J))*1E2*c;
        N2_nu_P1_47e = (N2_DEn_P1_47e(J))*1E2*c;
        N2_nu_P2_47e = (N2_DEn_P2_47e(J))*1E2*c;
        N2_nu_P3_47e = (N2_DEn_P3_47e(J))*1E2*c;
        N2_nu_Q1_47e = (N2_DEn_Q1_47e(J))*1E2*c;
        N2_nu_Q2_47e = (N2_DEn_Q2_47e(J))*1E2*c;
        N2_nu_Q3_47e = (N2_DEn_Q3_47e(J))*1E2*c;
        N2_nu_R1_47f = (N2_DEn_R1_47f(J))*1E2*c;
        N2_nu_R2_47f = (N2_DEn_R2_47f(J))*1E2*c;
        N2_nu_R3_47f = (N2_DEn_R3_47f(J))*1E2*c;
        N2_nu_P1_47f = (N2_DEn_P1_47f(J))*1E2*c;
        N2_nu_P2_47f = (N2_DEn_P2_47f(J))*1E2*c;
        N2_nu_P3_47f = (N2_DEn_P3_47f(J))*1E2*c;
        N2_nu_Q1_47f = (N2_DEn_Q1_47f(J))*1E2*c;
        N2_nu_Q2_47f = (N2_DEn_Q2_47f(J))*1E2*c;
        N2_nu_Q3_47f = (N2_DEn_Q3_47f(J))*1E2*c;    
        if(J>2)
            N2_lambda_Q1_47e = 1E7/N2_DEn_Q1_47e(J);
            N2_lambda_P1_47e = 1E7/N2_DEn_P1_47e(J);
            N2_lambda_R1_47e = 1E7/N2_DEn_R1_47e(J);
            N2_lambda_Q1_47f = 1E7/N2_DEn_Q1_47f(J);
            N2_lambda_P1_47f = 1E7/N2_DEn_P1_47f(J);
            N2_lambda_R1_47f = 1E7/N2_DEn_R1_47f(J);           
        end
        if(J>1)
            N2_lambda_Q2_47e = 1E7/N2_DEn_Q2_47e(J);
            N2_lambda_P2_47e = 1E7/N2_DEn_P2_47e(J);
            N2_lambda_R2_47e = 1E7/N2_DEn_R2_47e(J);
            N2_lambda_Q2_47f = 1E7/N2_DEn_Q2_47f(J);
            N2_lambda_P2_47f = 1E7/N2_DEn_P2_47f(J);
            N2_lambda_R2_47f = 1E7/N2_DEn_R2_47f(J);         
        end
        N2_lambda_Q3_47e = 1E7/N2_DEn_Q3_47e(J);            
        N2_lambda_P3_47e = 1E7/N2_DEn_P3_47e(J);
        N2_lambda_R3_47e = 1E7/N2_DEn_R3_47e(J);
        N2_lambda_Q3_47f = 1E7/N2_DEn_Q3_47f(J);            
        N2_lambda_P3_47f = 1E7/N2_DEn_P3_47f(J);
        N2_lambda_R3_47f = 1E7/N2_DEn_R3_47f(J);     

        % (1-4)
        N2_nu_R1_14e = (N2_DEn_R1_14e(J))*1E2*c;
        N2_nu_R2_14e = (N2_DEn_R2_14e(J))*1E2*c;
        N2_nu_R3_14e = (N2_DEn_R3_14e(J))*1E2*c;
        N2_nu_P1_14e = (N2_DEn_P1_14e(J))*1E2*c;
        N2_nu_P2_14e = (N2_DEn_P2_14e(J))*1E2*c;
        N2_nu_P3_14e = (N2_DEn_P3_14e(J))*1E2*c;
        N2_nu_Q1_14e = (N2_DEn_Q1_14e(J))*1E2*c;
        N2_nu_Q2_14e = (N2_DEn_Q2_14e(J))*1E2*c;
        N2_nu_Q3_14e = (N2_DEn_Q3_14e(J))*1E2*c;
        N2_nu_R1_14f = (N2_DEn_R1_14f(J))*1E2*c;
        N2_nu_R2_14f = (N2_DEn_R2_14f(J))*1E2*c;
        N2_nu_R3_14f = (N2_DEn_R3_14f(J))*1E2*c;
        N2_nu_P1_14f = (N2_DEn_P1_14f(J))*1E2*c;
        N2_nu_P2_14f = (N2_DEn_P2_14f(J))*1E2*c;
        N2_nu_P3_14f = (N2_DEn_P3_14f(J))*1E2*c;
        N2_nu_Q1_14f = (N2_DEn_Q1_14f(J))*1E2*c;
        N2_nu_Q2_14f = (N2_DEn_Q2_14f(J))*1E2*c;
        N2_nu_Q3_14f = (N2_DEn_Q3_14f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_14e = 1E7/N2_DEn_Q1_14e(J);
            N2_lambda_P1_14e = 1E7/N2_DEn_P1_14e(J);
            N2_lambda_R1_14e = 1E7/N2_DEn_R1_14e(J);
            N2_lambda_Q1_14f = 1E7/N2_DEn_Q1_14f(J);
            N2_lambda_P1_14f = 1E7/N2_DEn_P1_14f(J);
            N2_lambda_R1_14f = 1E7/N2_DEn_R1_14f(J);           
        end
        if(J>1)
            N2_lambda_Q2_14e = 1E7/N2_DEn_Q2_14e(J);
            N2_lambda_P2_14e = 1E7/N2_DEn_P2_14e(J);
            N2_lambda_R2_14e = 1E7/N2_DEn_R2_14e(J);
            N2_lambda_Q2_14f = 1E7/N2_DEn_Q2_14f(J);
            N2_lambda_P2_14f = 1E7/N2_DEn_P2_14f(J);
            N2_lambda_R2_14f = 1E7/N2_DEn_R2_14f(J);         
        end
        N2_lambda_Q3_14e = 1E7/N2_DEn_Q3_14e(J);            
        N2_lambda_P3_14e = 1E7/N2_DEn_P3_14e(J);
        N2_lambda_R3_14e = 1E7/N2_DEn_R3_14e(J);
        N2_lambda_Q3_14f = 1E7/N2_DEn_Q3_14f(J);            
        N2_lambda_P3_14f = 1E7/N2_DEn_P3_14f(J);
        N2_lambda_R3_14f = 1E7/N2_DEn_R3_14f(J);      
        
        % (0-3)
        N2_nu_R1_03e = (N2_DEn_R1_03e(J))*1E2*c;
        N2_nu_R2_03e = (N2_DEn_R2_03e(J))*1E2*c;
        N2_nu_R3_03e = (N2_DEn_R3_03e(J))*1E2*c;
        N2_nu_P1_03e = (N2_DEn_P1_03e(J))*1E2*c;
        N2_nu_P2_03e = (N2_DEn_P2_03e(J))*1E2*c;
        N2_nu_P3_03e = (N2_DEn_P3_03e(J))*1E2*c;
        N2_nu_Q1_03e = (N2_DEn_Q1_03e(J))*1E2*c;
        N2_nu_Q2_03e = (N2_DEn_Q2_03e(J))*1E2*c;
        N2_nu_Q3_03e = (N2_DEn_Q3_03e(J))*1E2*c;
        N2_nu_R1_03f = (N2_DEn_R1_03f(J))*1E2*c;
        N2_nu_R2_03f = (N2_DEn_R2_03f(J))*1E2*c;
        N2_nu_R3_03f = (N2_DEn_R3_03f(J))*1E2*c;
        N2_nu_P1_03f = (N2_DEn_P1_03f(J))*1E2*c;
        N2_nu_P2_03f = (N2_DEn_P2_03f(J))*1E2*c;
        N2_nu_P3_03f = (N2_DEn_P3_03f(J))*1E2*c;
        N2_nu_Q1_03f = (N2_DEn_Q1_03f(J))*1E2*c;
        N2_nu_Q2_03f = (N2_DEn_Q2_03f(J))*1E2*c;
        N2_nu_Q3_03f = (N2_DEn_Q3_03f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_03e = 1E7/N2_DEn_Q1_03e(J);
            N2_lambda_P1_03e = 1E7/N2_DEn_P1_03e(J);
            N2_lambda_R1_03e = 1E7/N2_DEn_R1_03e(J);
            N2_lambda_Q1_03f = 1E7/N2_DEn_Q1_03f(J);
            N2_lambda_P1_03f = 1E7/N2_DEn_P1_03f(J);
            N2_lambda_R1_03f = 1E7/N2_DEn_R1_03f(J);           
        end
        if(J>1)
            N2_lambda_Q2_03e = 1E7/N2_DEn_Q2_03e(J);
            N2_lambda_P2_03e = 1E7/N2_DEn_P2_03e(J);
            N2_lambda_R2_03e = 1E7/N2_DEn_R2_03e(J);
            N2_lambda_Q2_03f = 1E7/N2_DEn_Q2_03f(J);
            N2_lambda_P2_03f = 1E7/N2_DEn_P2_03f(J);
            N2_lambda_R2_03f = 1E7/N2_DEn_R2_03f(J);         
        end
        N2_lambda_Q3_03e = 1E7/N2_DEn_Q3_03e(J);            
        N2_lambda_P3_03e = 1E7/N2_DEn_P3_03e(J);
        N2_lambda_R3_03e = 1E7/N2_DEn_R3_03e(J);
        N2_lambda_Q3_03f = 1E7/N2_DEn_Q3_03f(J);            
        N2_lambda_P3_03f = 1E7/N2_DEn_P3_03f(J);
        N2_lambda_R3_03f = 1E7/N2_DEn_R3_03f(J); 
        
        % (0-2)
        N2_nu_R1_02e = (N2_DEn_R1_02e(J))*1E2*c;
        N2_nu_R2_02e = (N2_DEn_R2_02e(J))*1E2*c;
        N2_nu_R3_02e = (N2_DEn_R3_02e(J))*1E2*c;
        N2_nu_P1_02e = (N2_DEn_P1_02e(J))*1E2*c;
        N2_nu_P2_02e = (N2_DEn_P2_02e(J))*1E2*c;
        N2_nu_P3_02e = (N2_DEn_P3_02e(J))*1E2*c;
        N2_nu_Q1_02e = (N2_DEn_Q1_02e(J))*1E2*c;
        N2_nu_Q2_02e = (N2_DEn_Q2_02e(J))*1E2*c;
        N2_nu_Q3_02e = (N2_DEn_Q3_02e(J))*1E2*c;
        N2_nu_R1_02f = (N2_DEn_R1_02f(J))*1E2*c;
        N2_nu_R2_02f = (N2_DEn_R2_02f(J))*1E2*c;
        N2_nu_R3_02f = (N2_DEn_R3_02f(J))*1E2*c;
        N2_nu_P1_02f = (N2_DEn_P1_02f(J))*1E2*c;
        N2_nu_P2_02f = (N2_DEn_P2_02f(J))*1E2*c;
        N2_nu_P3_02f = (N2_DEn_P3_02f(J))*1E2*c;
        N2_nu_Q1_02f = (N2_DEn_Q1_02f(J))*1E2*c;
        N2_nu_Q2_02f = (N2_DEn_Q2_02f(J))*1E2*c;
        N2_nu_Q3_02f = (N2_DEn_Q3_02f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_02e = 1E7/N2_DEn_Q1_02e(J);
            N2_lambda_P1_02e = 1E7/N2_DEn_P1_02e(J);
            N2_lambda_R1_02e = 1E7/N2_DEn_R1_02e(J);
            N2_lambda_Q1_02f = 1E7/N2_DEn_Q1_02f(J);
            N2_lambda_P1_02f = 1E7/N2_DEn_P1_02f(J);
            N2_lambda_R1_02f = 1E7/N2_DEn_R1_02f(J);           
        end
        if(J>1)
            N2_lambda_Q2_02e = 1E7/N2_DEn_Q2_02e(J);
            N2_lambda_P2_02e = 1E7/N2_DEn_P2_02e(J);
            N2_lambda_R2_02e = 1E7/N2_DEn_R2_02e(J);
            N2_lambda_Q2_02f = 1E7/N2_DEn_Q2_02f(J);
            N2_lambda_P2_02f = 1E7/N2_DEn_P2_02f(J);
            N2_lambda_R2_02f = 1E7/N2_DEn_R2_02f(J);         
        end
        N2_lambda_Q3_02e = 1E7/N2_DEn_Q3_02e(J);            
        N2_lambda_P3_02e = 1E7/N2_DEn_P3_02e(J);
        N2_lambda_R3_02e = 1E7/N2_DEn_R3_02e(J);
        N2_lambda_Q3_02f = 1E7/N2_DEn_Q3_02f(J);            
        N2_lambda_P3_02f = 1E7/N2_DEn_P3_02f(J);
        N2_lambda_R3_02f = 1E7/N2_DEn_R3_02f(J);         
        
        % (1-3)
        N2_nu_R1_13e = (N2_DEn_R1_13e(J))*1E2*c;
        N2_nu_R2_13e = (N2_DEn_R2_13e(J))*1E2*c;
        N2_nu_R3_13e = (N2_DEn_R3_13e(J))*1E2*c;
        N2_nu_P1_13e = (N2_DEn_P1_13e(J))*1E2*c;
        N2_nu_P2_13e = (N2_DEn_P2_13e(J))*1E2*c;
        N2_nu_P3_13e = (N2_DEn_P3_13e(J))*1E2*c;
        N2_nu_Q1_13e = (N2_DEn_Q1_13e(J))*1E2*c;
        N2_nu_Q2_13e = (N2_DEn_Q2_13e(J))*1E2*c;
        N2_nu_Q3_13e = (N2_DEn_Q3_13e(J))*1E2*c;
        N2_nu_R1_13f = (N2_DEn_R1_13f(J))*1E2*c;
        N2_nu_R2_13f = (N2_DEn_R2_13f(J))*1E2*c;
        N2_nu_R3_13f = (N2_DEn_R3_13f(J))*1E2*c;
        N2_nu_P1_13f = (N2_DEn_P1_13f(J))*1E2*c;
        N2_nu_P2_13f = (N2_DEn_P2_13f(J))*1E2*c;
        N2_nu_P3_13f = (N2_DEn_P3_13f(J))*1E2*c;
        N2_nu_Q1_13f = (N2_DEn_Q1_13f(J))*1E2*c;
        N2_nu_Q2_13f = (N2_DEn_Q2_13f(J))*1E2*c;
        N2_nu_Q3_13f = (N2_DEn_Q3_13f(J))*1E2*c;       
        if(J>2)
            N2_lambda_Q1_13e = 1E7/N2_DEn_Q1_13e(J);
            N2_lambda_P1_13e = 1E7/N2_DEn_P1_13e(J);
            N2_lambda_R1_13e = 1E7/N2_DEn_R1_13e(J);
            N2_lambda_Q1_13f = 1E7/N2_DEn_Q1_13f(J);
            N2_lambda_P1_13f = 1E7/N2_DEn_P1_13f(J);
            N2_lambda_R1_13f = 1E7/N2_DEn_R1_13f(J);           
        end
        if(J>1)
            N2_lambda_Q2_13e = 1E7/N2_DEn_Q2_13e(J);
            N2_lambda_P2_13e = 1E7/N2_DEn_P2_13e(J);
            N2_lambda_R2_13e = 1E7/N2_DEn_R2_13e(J);
            N2_lambda_Q2_13f = 1E7/N2_DEn_Q2_13f(J);
            N2_lambda_P2_13f = 1E7/N2_DEn_P2_13f(J);
            N2_lambda_R2_13f = 1E7/N2_DEn_R2_13f(J);         
        end
        N2_lambda_Q3_13e = 1E7/N2_DEn_Q3_13e(J);            
        N2_lambda_P3_13e = 1E7/N2_DEn_P3_13e(J);
        N2_lambda_R3_13e = 1E7/N2_DEn_R3_13e(J);
        N2_lambda_Q3_13f = 1E7/N2_DEn_Q3_13f(J);            
        N2_lambda_P3_13f = 1E7/N2_DEn_P3_13f(J);
        N2_lambda_R3_13f = 1E7/N2_DEn_R3_13f(J);        
        
        % (2-4)
        N2_nu_R1_24e = (N2_DEn_R1_24e(J))*1E2*c;
        N2_nu_R2_24e = (N2_DEn_R2_24e(J))*1E2*c;
        N2_nu_R3_24e = (N2_DEn_R3_24e(J))*1E2*c;
        N2_nu_P1_24e = (N2_DEn_P1_24e(J))*1E2*c;
        N2_nu_P2_24e = (N2_DEn_P2_24e(J))*1E2*c;
        N2_nu_P3_24e = (N2_DEn_P3_24e(J))*1E2*c;
        N2_nu_Q1_24e = (N2_DEn_Q1_24e(J))*1E2*c;
        N2_nu_Q2_24e = (N2_DEn_Q2_24e(J))*1E2*c;
        N2_nu_Q3_24e = (N2_DEn_Q3_24e(J))*1E2*c;
        N2_nu_R1_24f = (N2_DEn_R1_24f(J))*1E2*c;
        N2_nu_R2_24f = (N2_DEn_R2_24f(J))*1E2*c;
        N2_nu_R3_24f = (N2_DEn_R3_24f(J))*1E2*c;
        N2_nu_P1_24f = (N2_DEn_P1_24f(J))*1E2*c;
        N2_nu_P2_24f = (N2_DEn_P2_24f(J))*1E2*c;
        N2_nu_P3_24f = (N2_DEn_P3_24f(J))*1E2*c;
        N2_nu_Q1_24f = (N2_DEn_Q1_24f(J))*1E2*c;
        N2_nu_Q2_24f = (N2_DEn_Q2_24f(J))*1E2*c;
        N2_nu_Q3_24f = (N2_DEn_Q3_24f(J))*1E2*c;  
        if(J>2)
            N2_lambda_Q1_24e = 1E7/N2_DEn_Q1_24e(J);
            N2_lambda_P1_24e = 1E7/N2_DEn_P1_24e(J);
            N2_lambda_R1_24e = 1E7/N2_DEn_R1_24e(J);
            N2_lambda_Q1_24f = 1E7/N2_DEn_Q1_24f(J);
            N2_lambda_P1_24f = 1E7/N2_DEn_P1_24f(J);
            N2_lambda_R1_24f = 1E7/N2_DEn_R1_24f(J);           
        end
        if(J>1)
            N2_lambda_Q2_24e = 1E7/N2_DEn_Q2_24e(J);
            N2_lambda_P2_24e = 1E7/N2_DEn_P2_24e(J);
            N2_lambda_R2_24e = 1E7/N2_DEn_R2_24e(J);
            N2_lambda_Q2_24f = 1E7/N2_DEn_Q2_24f(J);
            N2_lambda_P2_24f = 1E7/N2_DEn_P2_24f(J);
            N2_lambda_R2_24f = 1E7/N2_DEn_R2_24f(J);          
        end
        N2_lambda_Q3_24e = 1E7/N2_DEn_Q3_24e(J);            
        N2_lambda_P3_24e = 1E7/N2_DEn_P3_24e(J);
        N2_lambda_R3_24e = 1E7/N2_DEn_R3_24e(J);
        N2_lambda_Q3_24f = 1E7/N2_DEn_Q3_24f(J);            
        N2_lambda_P3_24f = 1E7/N2_DEn_P3_24f(J);
        N2_lambda_R3_24f = 1E7/N2_DEn_R3_24f(J);
        
        
            % caso geral N2+ 1 negativo e CN violeta
            SJR1 = SRDoublet1(J,0,0,0)/(2*(2*J+1));
            SJR2 = SRDoublet2(J,0,0,0)/(2*(2*J+1));
            %SJQ1 = SQDoublet1(J,0,0,0)/(2*(2*J+1));
            SJRQ21 = 2*((J+0.5)+1)/(4*(J-0.5)*(J+0.5))/(2*(2*J+1));
            
            % caso geral N2 2 positivo
            SJR1T = SRTriplet1(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJR2T = SRTriplet2(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJR3T = SRTriplet3(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ1T = SQTriplet1(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ2T = SQTriplet2(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            SJQ3T = SQTriplet3(J,1,N2_Y_B,N2_Y_C)/(3*(2*J+1));
            
            SJR1T_25 = SRTriplet1(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJR2T_25 = SRTriplet2(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJR3T_25 = SRTriplet3(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ1T_25 = SQTriplet1(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ2T_25 = SQTriplet2(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));
            SJQ3T_25 = SQTriplet3(J,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*J+1));         
            
            SJR1T_47 = SRTriplet1(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJR2T_47 = SRTriplet2(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJR3T_47 = SRTriplet3(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ1T_47 = SQTriplet1(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ2T_47 = SQTriplet2(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1));
            SJQ3T_47 = SQTriplet3(J,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*J+1)); 
            
            SJR1T_14 = SRTriplet1(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJR2T_14 = SRTriplet2(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJR3T_14 = SRTriplet3(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ1T_14 = SQTriplet1(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ2T_14 = SQTriplet2(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1));
            SJQ3T_14 = SQTriplet3(J,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*J+1)); 
            
            SJR1T_03 = SRTriplet1(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJR2T_03 = SRTriplet2(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJR3T_03 = SRTriplet3(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ1T_03 = SQTriplet1(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ2T_03 = SQTriplet2(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1));
            SJQ3T_03 = SQTriplet3(J,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*J+1)); 
            
            SJR1T_02 = SRTriplet1(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJR2T_02 = SRTriplet2(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJR3T_02 = SRTriplet3(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ1T_02 = SQTriplet1(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ2T_02 = SQTriplet2(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1));
            SJQ3T_02 = SQTriplet3(J,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*J+1)); 
            
            SJR1T_13 = SRTriplet1(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJR2T_13 = SRTriplet2(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJR3T_13 = SRTriplet3(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ1T_13 = SQTriplet1(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ2T_13 = SQTriplet2(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1));
            SJQ3T_13 = SQTriplet3(J,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*J+1)); 
            
            SJR1T_24 = SRTriplet1(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJR2T_24 = SRTriplet2(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJR3T_24 = SRTriplet3(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ1T_24 = SQTriplet1(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ2T_24 = SQTriplet2(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            SJQ3T_24 = SQTriplet3(J,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*J+1));
            
            % (nu 0-0) transition in N2(+) R1, R2 and RQ21 branches
            
            sum = sum + N2p_ratio*A_00*gR*SJR1*nu_B1*exp(-c*h_p*En_B1(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B1).^2*log(2)/delta_lambda_inst^2);
            sum = sum + N2p_ratio*A_00*gR*SJR2*nu_B2*exp(-c*h_p*En_B2(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B2).^2*log(2)/delta_lambda_inst^2);
            sum = sum + N2p_ratio*A_00*gR*SJRQ21*nu_RQ21*exp(-c*h_p*En_B2(J)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_RQ21).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in N2(+) R1, R2 and RQ21 branches
            
            sum = sum + N2p_ratio*A_11*gP*SJR1*nu_B1_1*exp(-c*h_p*((En_B1_1(J)-DT_1)*1E2/kb/Trot)).*...
                exp(-(lambda-lambda_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            sum = sum + N2p_ratio*A_11*gP*SJR2*nu_B2_1*exp(-c*h_p*(En_B2_1(J)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            sum = sum + N2p_ratio*A_11*gP*SJRQ21*nu_RQ21_1*exp(-c*h_p*(En_B2_1(J)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_RQ21_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*(DT_1*1E2/kb/Tvib));
            
            % (nu 3-6) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR1T*N2_nu_R1e*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR1T*N2_nu_R1f*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ1T*N2_nu_Q1e*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ1T*N2_nu_Q1f*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));     
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR2T*N2_nu_R2e*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR2T*N2_nu_R2f*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));        
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T*N2_nu_Q2e*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T*N2_nu_Q2f*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            end
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR3T*N2_nu_R3e*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gR_N2*SJR3T*N2_nu_R3f*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ3T*N2_nu_Q3e*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ3T*N2_nu_Q3f*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            % (nu 2-5) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR1T_25*N2_nu_R1_25e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR1T_25*N2_nu_R1_25f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_25*gQ_N2*SJQ1T_25*N2_nu_Q1_25e*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_25*gQ_N2*SJQ1T_25*N2_nu_Q1_25f*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));        
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR2T_25*N2_nu_R2_25e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR2T_25*N2_nu_R2_25f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));            
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T_25*N2_nu_Q2_25e*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T_25*N2_nu_Q2_25f*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR3T_25*N2_nu_R3_25e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_25*gR_N2*SJR3T_25*N2_nu_R3_25f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));         
            sum = sum + N2_ratio*A_N2_25*gQ_N2*SJQ3T_25*N2_nu_Q3_25e*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_25*gQ_N2*SJQ3T_25*N2_nu_Q3_25f*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            
            % (nu 4-7) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_47*gR_N2*SJR1T_47*N2_nu_R1_47e*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_47*gR_N2*SJR1T_47*N2_nu_R1_47f*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ1T_47*N2_nu_Q1_47e*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ1T_47*N2_nu_Q1_47f*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_47*gR_N2*SJR2T_47*N2_nu_R2_47e*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gR_N2*SJR2T_47*N2_nu_R2_47f*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ2T_47*N2_nu_Q2_47e*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ2T_47*N2_nu_Q2_47f*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            end
            sum = sum + N2_ratio*A_N2_47*gR_N2*SJR3T_47*N2_nu_R3_47e*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gR_N2*SJR3T_47*N2_nu_R3_47f*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ3T_47*N2_nu_Q3_47e*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_47*gQ_N2*SJQ3T_47*N2_nu_Q3_47f*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));             
            % (nu 1-4) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR1T_14*N2_nu_R1_14e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR1T_14*N2_nu_R1_14f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ1T_14*N2_nu_Q1_14f*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ1T_14*N2_nu_Q1_14e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR2T_14*N2_nu_R2_14e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR2T_14*N2_nu_R2_14f*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ2T_14*N2_nu_Q2_14e*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ2T_14*N2_nu_Q2_14f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));         
            end
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR3T_14*N2_nu_R3_14e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gR_N2*SJR3T_14*N2_nu_R3_14f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ3T_14*N2_nu_Q3_14f*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gQ_N2*SJQ3T_14*N2_nu_Q3_14e*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));         
            
            % (nu 0-3) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR1T_03*N2_nu_R1_03e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR1T_03*N2_nu_R1_03f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ1T_03*N2_nu_Q1_03f*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ1T_03*N2_nu_Q1_03e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR2T_03*N2_nu_R2_03e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR2T_03*N2_nu_R2_03f*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ2T_03*N2_nu_Q2_03e*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ2T_03*N2_nu_Q2_03f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));         
            end
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR3T_03*N2_nu_R3_03e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gR_N2*SJR3T_03*N2_nu_R3_03f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ3T_03*N2_nu_Q3_03f*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gQ_N2*SJQ3T_03*N2_nu_Q3_03e*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));    
            
            % (nu 0-2) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR1T_02*N2_nu_R1_02e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR1T_02*N2_nu_R1_02f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ1T_02*N2_nu_Q1_02f*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ1T_02*N2_nu_Q1_02e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR2T_02*N2_nu_R2_02e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR2T_02*N2_nu_R2_02f*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ2T_02*N2_nu_Q2_02e*exp(-c*h_p*(N2_En_C2_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ2T_02*N2_nu_Q2_02f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));         
            end
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR3T_02*N2_nu_R3_02e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gR_N2*SJR3T_02*N2_nu_R3_02f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ3T_02*N2_nu_Q3_02f*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gQ_N2*SJQ3T_02*N2_nu_Q3_02e*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));   
            
            % (nu 1-3) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR1T_13*N2_nu_R1_13e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR1T_13*N2_nu_R1_13f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ1T_13*N2_nu_Q1_13f*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ1T_13*N2_nu_Q1_13e*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR2T_13*N2_nu_R2_13e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR2T_13*N2_nu_R2_13f*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ2T_13*N2_nu_Q2_13e*exp(-c*h_p*(N2_En_C2_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ2T_13*N2_nu_Q2_13f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));         
            end
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR3T_13*N2_nu_R3_13e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gR_N2*SJR3T_13*N2_nu_R3_13f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ3T_13*N2_nu_Q3_13f*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gQ_N2*SJQ3T_13*N2_nu_Q3_13e*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));  
            
            % (nu 2-4) transition in N2 R1, R2, R3, Q1, Q2 and Q3 branches
            % R
            if(J>2)
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR1T_24*N2_nu_R1_24e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR1T_24*N2_nu_R1_24f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_24*gQ_N2*SJQ1T_24*N2_nu_Q1_24e*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gQ_N2*SJQ1T_24*N2_nu_Q1_24f*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));        
            end
            if(J>1)
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR2T_24*N2_nu_R2_24e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR2T_24*N2_nu_R2_24f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));            
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T_24*N2_nu_Q2_24e*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gQ_N2*SJQ2T_24*N2_nu_Q2_24f*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR3T_24*N2_nu_R3_24e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gR_N2*SJR3T_24*N2_nu_R3_24f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_R3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));         
            sum = sum + N2_ratio*A_N2_24*gQ_N2*SJQ3T_24*N2_nu_Q3_24e*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gQ_N2*SJQ3T_24*N2_nu_Q3_24f*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_Q3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2)); 
            
            % (nu 0-0) transition in CN R1, R2 and RQ21 branches
            
            sum = sum + gCN*CN_ratio*A_CN_00*SJR1*CN_nu_B1*exp(-c*h_p*CN_En_B1(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1).^2*log(2)/delta_lambda_inst^2);
            sum = sum + gCN*CN_ratio*A_CN_00*SJR2*CN_nu_B2*exp(-c*h_p*CN_En_B2(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2).^2*log(2)/delta_lambda_inst^2);
            sum = sum + gCN*CN_ratio*A_CN_00*SJRQ21*CN_nu_RQ21*exp(-c*h_p*CN_En_B2(J)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in CN R1, R2 and RQ21 branches
        
            sum = sum + gCN*CN_ratio*A_CN_11*SJR1*CN_nu_B1_1*exp(-c*h_p*(CN_En_B1_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
            sum = sum + gCN*CN_ratio*A_CN_11*SJR2*CN_nu_B2_1*exp(-c*h_p*(CN_En_B2_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
            sum = sum + gCN*CN_ratio*A_CN_11*SJRQ21*CN_nu_RQ21_1*exp(-c*h_p*(CN_En_B2_1(J)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);
        % (nu 2-2) transition in CN R1, R2 and RQ21 branches
        
            sum = sum + gCN*CN_ratio*A_CN_22*SJR1*CN_nu_B1_2*exp(-c*h_p*(CN_En_B1_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
            sum = sum + gCN*CN_ratio*A_CN_22*SJR2*CN_nu_B2_2*exp(-c*h_p*(CN_En_B2_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
            sum = sum + gCN*CN_ratio*A_CN_22*SJRQ21*CN_nu_RQ21_2*exp(-c*h_p*(CN_En_B2_2(J)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);
        % (nu 3-3) transition in CN R1, R2 and RQ21 branches
        
            sum = sum + gCN*CN_ratio*A_CN_33*SJR1*CN_nu_B1_3*exp(-c*h_p*(CN_En_B1_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);
            sum = sum + gCN*CN_ratio*A_CN_33*SJR2*CN_nu_B2_3*exp(-c*h_p*(CN_En_B2_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);
            sum = sum + gCN*CN_ratio*A_CN_33*SJRQ21*CN_nu_RQ21_3*exp(-c*h_p*(CN_En_B2_3(J)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);    
          % (nu 4-4) transition in CN R1, R2 and RQ21 branches
         
            sum = sum + gCN*CN_ratio*A_CN_44*SJR1*CN_nu_B1_4*exp(-c*h_p*(CN_En_B1_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);
            sum = sum + gCN*CN_ratio*A_CN_44*SJR2*CN_nu_B2_4*exp(-c*h_p*(CN_En_B2_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);
            sum = sum + gCN*CN_ratio*A_CN_44*SJRQ21*CN_nu_RQ21_4*exp(-c*h_p*(CN_En_B2_4(J)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);            
        % (nu 5-5) transition in CN R1, R2 and RQ21 branches
        
            sum = sum + gCN*CN_ratio*A_CN_55*SJR1*CN_nu_B1_5*exp(-c*h_p*(CN_En_B1_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
            sum = sum + gCN*CN_ratio*A_CN_55*SJR2*CN_nu_B2_5*exp(-c*h_p*(CN_En_B2_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
            sum = sum + gCN*CN_ratio*A_CN_55*SJRQ21*CN_nu_RQ21_5*exp(-c*h_p*(CN_En_B2_5(J)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);
        % (nu 6-6) transition in CN R1, R2 and RQ21 branches
        
            sum = sum + gCN*CN_ratio*A_CN_66*SJR1*CN_nu_B1_6*exp(-c*h_p*(CN_En_B1_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);
            sum = sum + gCN*CN_ratio*A_CN_66*SJR2*CN_nu_B2_6*exp(-c*h_p*(CN_En_B2_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);
            sum = sum + gCN*CN_ratio*A_CN_66*SJRQ21*CN_nu_RQ21_6*exp(-c*h_p*(CN_En_B2_6(J)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);
            
        % (nu 7-7) transition in CN R1, R2 and RQ21 branches
            
            sum = sum + gCN*CN_ratio*A_CN_77*SJR1*CN_nu_B1_7*exp(-c*h_p*(CN_En_B1_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7); 
            sum = sum + gCN*CN_ratio*A_CN_77*SJR2*CN_nu_B2_7*exp(-c*h_p*(CN_En_B2_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);
            sum = sum + gCN*CN_ratio*A_CN_77*SJRQ21*CN_nu_RQ21_7*exp(-c*h_p*(CN_En_B2_7(J)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);
            
        % (nu 8-8) transition in CN R1, R2 and RQ21 branches
            
            sum = sum + gCN*CN_ratio*A_CN_88*SJR1*CN_nu_B1_8*exp(-c*h_p*(CN_En_B1_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            sum = sum + gCN*CN_ratio*A_CN_88*SJR2*CN_nu_B2_8*exp(-c*h_p*(CN_En_B2_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            sum = sum + gCN*CN_ratio*A_CN_88*SJRQ21*CN_nu_RQ21_8*exp(-c*h_p*(CN_En_B2_8(J)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);
            
        % (nu 9-9) transition in CN R1, R2 and RQ21 branches
            
            sum = sum + gCN*CN_ratio*A_CN_99*SJR1*CN_nu_B1_9*exp(-c*h_p*(CN_En_B1_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B1_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            sum = sum + gCN*CN_ratio*A_CN_99*SJR2*CN_nu_B2_9*exp(-c*h_p*(CN_En_B2_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_B2_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            sum = sum + gCN*CN_ratio*A_CN_99*SJRQ21*CN_nu_RQ21_9*exp(-c*h_p*(CN_En_B2_9(J)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_RQ21_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);
            
        % (nu 10-10) transition in CN R1, R2 and RQ21 branches
            
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJR1*CN_nu_B1_10*exp(-c*h_p*(CN_En_B1_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_B1_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJR2*CN_nu_B2_10*exp(-c*h_p*(CN_En_B2_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_B2_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJRQ21*CN_nu_RQ21_10*exp(-c*h_p*(CN_En_B2_10(J)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_RQ21_10).^2*log(2)/delta_lambda_inst^2).*...
%                 exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);

        if(J>1)
            SJP1 = SPDoublet1(J-1,0,0,0)/(2*(2*(J-1)+1));
            SJP2 = SPDoublet2(J-1,0,0,0)/(2*(2*(J-1)+1));
            %SJQ2 = SQDoublet2(J-1,0,0,0)/(2*(2*(J-1)+1));
            SJPQ12 = (2*(J+0.5)+1)/(4*(J-0.5)*(J+0.5))/(2*(2*(J-1)+1));
            
            SJP1T = SPTriplet1(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            SJP2T = SPTriplet2(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            SJP3T = SPTriplet3(J-1,1,N2_Y_B,N2_Y_C)/(3*(2*(J-1)+1));
            
            SJP1T_25 = SPTriplet1(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP2T_25 = SPTriplet2(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP3T_25 = SPTriplet3(J-1,1,N2_Y_B_5,N2_Y_C_2)/(3*(2*(J-1)+1));  
            
            SJP1T_47 = SPTriplet1(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));
            SJP2T_47 = SPTriplet2(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));
            SJP3T_47 = SPTriplet3(J-1,1,N2_Y_B_7,N2_Y_C_4)/(3*(2*(J-1)+1));  
            
            SJP1T_14 = SPTriplet1(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP2T_14 = SPTriplet2(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP3T_14 = SPTriplet3(J-1,1,N2_Y_B_4,N2_Y_C_1)/(3*(2*(J-1)+1));
            
            SJP1T_03 = SPTriplet1(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP2T_03 = SPTriplet2(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP3T_03 = SPTriplet3(J-1,1,N2_Y_B_3,N2_Y_C_0)/(3*(2*(J-1)+1)); 
            
            SJP1T_02 = SPTriplet1(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP2T_02 = SPTriplet2(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            SJP3T_02 = SPTriplet3(J-1,1,N2_Y_B_2,N2_Y_C_0)/(3*(2*(J-1)+1));
            
            SJP1T_13 = SPTriplet1(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP2T_13 = SPTriplet2(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));
            SJP3T_13 = SPTriplet3(J-1,1,N2_Y_B_3,N2_Y_C_1)/(3*(2*(J-1)+1));

            SJP1T_24 = SPTriplet1(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP2T_24 = SPTriplet2(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            SJP3T_24 = SPTriplet3(J-1,1,N2_Y_B_4,N2_Y_C_2)/(3*(2*(J-1)+1));
            
            % (nu 0-0) transition in N2(+) P1, P2 and PQ12 branches
            sum = sum + N2p_ratio*A_00*gP*SJP1*nu_P_B1*exp(-c*h_p*En_B1(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B1).^2*log(2)/delta_lambda_inst^2); 
            sum = sum + N2p_ratio*A_00*gP*SJP2*nu_P_B2*exp(-c*h_p*En_B2(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B2).^2*log(2)/delta_lambda_inst^2);
             sum = sum + N2p_ratio*A_00*gP*SJPQ12*nu_PQ12*exp(-c*h_p*En_B1(J-1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_PQ12).^2*log(2)/delta_lambda_inst^2); 
            
            % (nu 1-1) transition in N2(+) P1, P2 and PQ12 branches
            sum = sum + N2p_ratio*A_11*gR*SJP1*nu_P_B1_1*exp(-c*h_p*(En_B1_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib); 
            sum = sum + N2p_ratio*A_11*gR*SJP2*nu_P_B2_1*exp(-c*h_p*(En_B2_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_P_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib);
            sum = sum + N2p_ratio*A_11*gR*SJPQ12*nu_PQ12_1*exp(-c*h_p*(En_B1_1(J-1)-DT_1)*1E2/kb/Trot).*...
                exp(-(lambda-lambda_PQ12_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(-c*h_p*DT_1*1E2/kb/Tvib);  
            
            % (nu 3-6) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP1T*N2_nu_P1e*exp(-c*h_p*(N2_En_C1_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP1T*N2_nu_P1f*exp(-c*h_p*(N2_En_C1_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));           
            end
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP2T*N2_nu_P2e*exp(-c*h_p*(N2_En_C2_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP2T*N2_nu_P2f*exp(-c*h_p*(N2_En_C2_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));            
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP3T*N2_nu_P3e*exp(-c*h_p*(N2_En_C3_3e(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2)); 
            sum = sum + N2_ratio*A_N2_36*gP_N2*SJP3T*N2_nu_P3f*exp(-c*h_p*(N2_En_C3_3f(J)-V_C_3)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_3*1E2/kb/TvibN2));          
            
            % (nu 2-5) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_25*gP_N2*SJP1T_25*N2_nu_P1_25e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_25*gP_N2*SJP1T_25*N2_nu_P1_25f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            end
            sum = sum + N2_ratio*A_N2_25*gP_N2*SJP2T_25*N2_nu_P2_25e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_25*gP_N2*SJP2T_25*N2_nu_P2_25f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_25*gP_N2*SJP3T_25*N2_nu_P3_25e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_25e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));   
            sum = sum + N2_ratio*A_N2_25*gP_N2*SJP3T_25*N2_nu_P3_25f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_25f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));             
            
            % (nu 4-7) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP1T_47*N2_nu_P1_47e*exp(-c*h_p*(N2_En_C1_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP1T_47*N2_nu_P1_47f*exp(-c*h_p*(N2_En_C1_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP2T_47*N2_nu_P2_47e*exp(-c*h_p*(N2_En_C2_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP2T_47*N2_nu_P2_47f*exp(-c*h_p*(N2_En_C2_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));           
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP3T_47*N2_nu_P3_47e*exp(-c*h_p*(N2_En_C3_4e(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_47e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2)); 
            sum = sum + N2_ratio*A_N2_47*gP_N2*SJP3T_47*N2_nu_P3_47f*exp(-c*h_p*(N2_En_C3_4f(J)-V_C_4)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_47f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_4*1E2/kb/TvibN2));            
            
            % (nu 1-4) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP1T_14*N2_nu_P1_14e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP1T_14*N2_nu_P1_14f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP2T_14*N2_nu_P2_14e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP2T_14*N2_nu_P2_14f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP3T_14*N2_nu_P3_14e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_14e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));       
            sum = sum + N2_ratio*A_N2_14*gP_N2*SJP3T_14*N2_nu_P3_14f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_14f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            
            % (nu 0-3) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP1T_03*N2_nu_P1_03e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP1T_03*N2_nu_P1_03f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP2T_03*N2_nu_P2_03e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP2T_03*N2_nu_P2_03f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP3T_03*N2_nu_P3_03e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_03e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));       
            sum = sum + N2_ratio*A_N2_03*gP_N2*SJP3T_03*N2_nu_P3_03f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_03f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2)); 
            
            % (nu 0-2) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP1T_02*N2_nu_P1_02e*exp(-c*h_p*(N2_En_C1_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP1T_02*N2_nu_P1_02f*exp(-c*h_p*(N2_En_C1_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP2T_02*N2_nu_P2_02e*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP2T_02*N2_nu_P2_02f*exp(-c*h_p*(N2_En_C2_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP3T_02*N2_nu_P3_02e*exp(-c*h_p*(N2_En_C3_0e(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_02e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));       
            sum = sum + N2_ratio*A_N2_02*gP_N2*SJP3T_02*N2_nu_P3_02f*exp(-c*h_p*(N2_En_C3_0f(J)-V_C_0)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_02f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_0*1E2/kb/TvibN2));   
            
            % (nu 1-3) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP1T_13*N2_nu_P1_13e*exp(-c*h_p*(N2_En_C1_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP1T_13*N2_nu_P1_13f*exp(-c*h_p*(N2_En_C1_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            end
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP2T_13*N2_nu_P2_13e*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP2T_13*N2_nu_P2_13f*exp(-c*h_p*(N2_En_C2_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP3T_13*N2_nu_P3_13e*exp(-c*h_p*(N2_En_C3_1e(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_13e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));       
            sum = sum + N2_ratio*A_N2_13*gP_N2*SJP3T_13*N2_nu_P3_13f*exp(-c*h_p*(N2_En_C3_1f(J)-V_C_1)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_13f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_1*1E2/kb/TvibN2));
            
            % (nu 2-4) transition in N2 P1, P2, P3 branches
            % P
            if(J>2)
            sum = sum + N2_ratio*A_N2_24*gP_N2*SJP1T_24*N2_nu_P1_24e*exp(-c*h_p*(N2_En_C1_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
             sum = sum + N2_ratio*A_N2_24*gP_N2*SJP1T_24*N2_nu_P1_24f*exp(-c*h_p*(N2_En_C1_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P1_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));           
            end
            sum = sum + N2_ratio*A_N2_24*gP_N2*SJP2T_24*N2_nu_P2_24e*exp(-c*h_p*(N2_En_C2_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));
            sum = sum + N2_ratio*A_N2_24*gP_N2*SJP2T_24*N2_nu_P2_24f*exp(-c*h_p*(N2_En_C2_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P2_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));          
            sum = sum + N2_ratio*A_N2_24*gP_N2*SJP3T_24*N2_nu_P3_24e*exp(-c*h_p*(N2_En_C3_2e(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_24e).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2));   
            sum = sum + N2_ratio*A_N2_24*gP_N2*SJP3T_24*N2_nu_P3_24f*exp(-c*h_p*(N2_En_C3_2f(J)-V_C_2)*1E2/kb/TrotN2).*...
                exp(-(lambda-N2_lambda_P3_24f).^2*log(2)/delta_lambda_inst^2)*...
                exp(-c*h_p*(V_C_2*1E2/kb/TvibN2)); 
            
            % (nu 0-0) transition in CN, P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_00*SJP1*CN_nu_P_B1*exp(-c*h_p*CN_En_B1(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1).^2*log(2)/delta_lambda_inst^2);
            sum = sum + gCN*CN_ratio*A_CN_00*SJP2*CN_nu_P_B2*exp(-c*h_p*CN_En_B2(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2).^2*log(2)/delta_lambda_inst^2);
            sum = sum + gCN*CN_ratio*A_CN_00*SJPQ12*CN_nu_PQ12*exp(-c*h_p*CN_En_B1(J-1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12).^2*log(2)/delta_lambda_inst^2);
            
            % (nu 1-1) transition in CN P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_11*SJP1*CN_nu_P_B1_1*exp(-c*h_p*(CN_En_B1_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1); 
            sum = sum + gCN*CN_ratio*A_CN_11*SJP2*CN_nu_P_B2_1*exp(-c*h_p*(CN_En_B2_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1); 
            sum = sum + gCN*CN_ratio*A_CN_11*SJPQ12*CN_nu_PQ12_1*exp(-c*h_p*(CN_En_B1_1(J-1)-CN_DT1)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_1).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT1);        
            
            % (nu 2-2) transition in CN P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_22*SJP1*CN_nu_P_B1_2*exp(-c*h_p*(CN_En_B1_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2);  
            sum = sum + gCN*CN_ratio*A_CN_22*SJP2*CN_nu_P_B2_2*exp(-c*h_p*(CN_En_B2_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2); 
            sum = sum + gCN*CN_ratio*A_CN_22*SJPQ12*CN_nu_PQ12_2*exp(-c*h_p*(CN_En_B1_2(J-1)-CN_DT2)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_2).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT2); 
            % (nu 3-3) transition in CN P1, P2 and PQ12 branches
            
            sum = sum + gCN*CN_ratio*A_CN_33*SJP1*CN_nu_P_B1_3*exp(-c*h_p*(CN_En_B1_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3);  
            sum = sum + gCN*CN_ratio*A_CN_33*SJP2*CN_nu_P_B2_3*exp(-c*h_p*(CN_En_B2_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3); 
            sum = sum + gCN*CN_ratio*A_CN_33*SJPQ12*CN_nu_PQ12_3*exp(-c*h_p*(CN_En_B1_3(J-1)-CN_DT3)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_3).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT3); 
        
            % (nu 4-4) transition in CN P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_44*SJP1*CN_nu_P_B1_4*exp(-c*h_p*(CN_En_B1_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);  
            sum = sum + gCN*CN_ratio*A_CN_44*SJP2*CN_nu_P_B2_4*exp(-c*h_p*(CN_En_B2_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4); 
            sum = sum + gCN*CN_ratio*A_CN_44*SJPQ12*CN_nu_PQ12_4*exp(-c*h_p*(CN_En_B1_4(J-1)-CN_DT4)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_4).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT4);    
        
            % (nu 5-5) transition in CN P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_55*SJP1*CN_nu_P_B1_5*exp(-c*h_p*(CN_En_B1_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5);  
            sum = sum + gCN*CN_ratio*A_CN_55*SJP2*CN_nu_P_B2_5*exp(-c*h_p*(CN_En_B2_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5); 
            sum = sum + gCN*CN_ratio*A_CN_55*SJPQ12*CN_nu_PQ12_5*exp(-c*h_p*(CN_En_B1_5(J-1)-CN_DT5)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_5).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT5); 
        
            % (nu 6-6) transition in CN P1, P2 and PQ12 branches
            sum = sum + gCN*CN_ratio*A_CN_66*SJP1*CN_nu_P_B1_6*exp(-c*h_p*(CN_En_B1_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);  
            sum = sum + gCN*CN_ratio*A_CN_66*SJP2*CN_nu_P_B2_6*exp(-c*h_p*(CN_En_B2_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6);  
            sum = sum + gCN*CN_ratio*A_CN_66*SJPQ12*CN_nu_PQ12_6*exp(-c*h_p*(CN_En_B1_6(J-1)-CN_DT6)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_6).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT6); 

            % (nu 7-7) transition in CN P1, P2 and PQ12 branches
            
            sum = sum + gCN*CN_ratio*A_CN_77*SJP1*CN_nu_P_B1_7*exp(-c*h_p*(CN_En_B1_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);   
            sum = sum + gCN*CN_ratio*A_CN_77*SJP2*CN_nu_P_B2_7*exp(-c*h_p*(CN_En_B2_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);  
            sum = sum + gCN*CN_ratio*A_CN_77*SJPQ12*CN_nu_PQ12_7*exp(-c*h_p*(CN_En_B1_7(J-1)-CN_DT7)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_7).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT7);    
        
            % (nu 8-8) transition in CN P1, P2 and PQ12 branches
            
            sum = sum + gCN*CN_ratio*A_CN_88*SJP1*CN_nu_P_B1_8*exp(-c*h_p*(CN_En_B1_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);   
            sum = sum + gCN*CN_ratio*A_CN_88*SJP2*CN_nu_P_B2_8*exp(-c*h_p*(CN_En_B2_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8);  
            sum = sum + gCN*CN_ratio*A_CN_88*SJPQ12*CN_nu_PQ12_8*exp(-c*h_p*(CN_En_B1_8(J-1)-CN_DT8)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_8).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT8); 
        
            % (nu 9-9) transition in CN P1, P2 and PQ12 branches
            
            sum = sum + gCN*CN_ratio*A_CN_99*SJP1*CN_nu_P_B1_9*exp(-c*h_p*(CN_En_B1_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B1_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9);   
            sum = sum + gCN*CN_ratio*A_CN_99*SJP2*CN_nu_P_B2_9*exp(-c*h_p*(CN_En_B2_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_P_B2_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9); 
            sum = sum + gCN*CN_ratio*A_CN_99*SJPQ12*CN_nu_PQ12_9*exp(-c*h_p*(CN_En_B1_9(J-1)-CN_DT9)*1E2/kb/TrotCN).*...
                exp(-(lambda-CN_lambda_PQ12_9).^2*log(2)/delta_lambda_inst^2).*...
                exp(f_DT9); 
        
            % (nu 10-10) transition in CN P1, P2 and PQ12 branches
            
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJP1*CN_nu_P_B1_10*exp(-c*h_p*(CN_En_B1_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_P_B1_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);  
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJP2*CN_nu_P_B2_10*exp(-c*h_p*(CN_En_B2_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_P_B2_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN); 
%             sum = sum + gCN*CN_ratio*A_CN_1010*SJPQ12*CN_nu_PQ12_10*exp(-c*h_p*(CN_En_B1_10(J-1)-CN_DT10)*1E2/kb/TrotCN).*...
%                 exp(-(lambda-CN_lambda_PQ12_10).^2*log(2)/delta_lambda_inst^2).*...
%             exp(-c*h_p*CN_DT10*1E2/kb/TvibCN);     
        end
    end
    I = sum/max(sum);
    % add artificial gaussian shape profile to model the unknown line 
%    I = I + beta(9)*exp(-(lambda-395.028).^2*log(2)/delta_lambda_inst^2);
%    I = I + beta(10)*exp(-(lambda-394.89).^2*log(2)/delta_lambda_inst^2);
    %    beta(6)*exp(-(lambda-383.3).^2*log(2)/beta(7)^2);
end

function [ beta, Res ] = fitGM( X, y, s, model, beta0, maxiter)
% [ beta, Res ] = mmqGM( X, y, s, model, beta0 )
% Argumentos de entrada:  
%   beta0: parâmetros iniciais
%
%
%
%  inspirado na rotina do nlinfit do toolbox do matlab
%  Primeira versao: Otaviano (14-9-01)
%  Ultima modificacao: Marco (05-05-09)

n = length(y);
%Transforma os vetores de dados em vetores coluna
y = y(:);
X = X(:);

p = length(beta0);
%Transforma os vetores de parâmetros iniciais em vetores coluna
beta0 = beta0(:);
if length(s)==1
    S = s*ones(n,1);
else
    S = s;
end

J = zeros(n,p);
beta = beta0;
betanew = beta + 1;
%maxiter = 40; % ==============
iter = 0;
betatol = 1.0E-2;
rtol = 1.0E-4;
qui2 = 1;
qui2old = qui2;
S = S(:);
v = eye(p);
convergiu = 0;
fatorGM = 1;
% disp( 'new fit' )
usandoGM = 0;
while and( convergiu<2, iter < maxiter )
    
    if iter > 0
        beta = betanew;
        yfit = yfitnew;
        r = rnew;
        qui2old = qui2;
    else
        yfit = feval(model,beta,X);
        r = ( y(:) - yfit(:) ) ./ S;
        qui2old = r'*r;
    end
   iter = iter + 1;
   J = zeros(size(y,1), p);
   for k = 1:p,
       delta = zeros( size( beta) );
       delta(k) = sqrt(eps)*beta(k);
       if abs(delta(k)) < 100*eps
          delta(k) = 100*eps;
       end
       if (abs(delta(k)) < 0.001*abs(beta(k)))
           delta(k) = 0.001*beta(k);
       end
       betaplus = beta + delta;
       betaminus = beta - delta;
       yplus = feval(model, betaplus , X);
       yminus = feval(model, betaminus, X);
       J(:,k) = ( ( yplus(:) - yminus(:) ) ./ S) / (delta(k));
       %J(:,k) = ( ( yplus(:) - yminus(:) ) ./ s ) / (betaplus(k)-betaminus(k));
   end
   JJ = J'*J;

   JJplus = JJ + (diag(diag(JJ)))*fatorGM;
   step = inverta(JJplus)*(J'*r);
   

   betanew = beta + step;
   yfitnew = feval(model,betanew,X);
   rnew = (y(:) - yfitnew(:))./S;
   
   qui2 = rnew'*rnew;
   if ( qui2 > qui2old*(1+1E-4))
       convergiu = 0;
       while and( qui2>qui2old*(1+1E-4), fatorGM<=1E2 )
           fatorGM = fatorGM * sqrt(10);
           JJplus = JJ + diag(diag(JJ))*fatorGM;
           step = inverta(JJplus)*(J'*r);
           betanew = beta + step;
           yfitnew = feval(model,betanew,X);
           rnew = (y(:) - yfitnew(:))./S;
           qui2 = rnew'*rnew;
       end
       fatorGM = 1;
       %       disp( sprintf( '%3d', iter) );
   else
       fatorGM = fatorGM / 10;
       if abs( qui2old - qui2 ) < 1E-5*(qui2+qui2old) %originalmente 1E-5*
           v = inverta( JJ )*eye(p);
           if max( abs( step )./sqrt( diag(v)*(qui2/(n-p)) ) )<1E-1 
               convergiu = convergiu +1;
               fatorGM = 1E-6;
           end
       end
   end
end
if convergiu==0, v = inverta( JJ )*eye(p); end
if iter == maxiter
   disp('mmqGM nao convergiu. Retornando os resultados da ultima iteracao.');
end
Res.A = beta;
Res.SA = sqrt( diag(v) );
Res.VA = v;
Res.qui2 = qui2;
Res.pqui2 = 100*pqui2(qui2, n-p);
Res.ngl = n-p;
Res.F = yfit;
Res.iter = iter;
Res.si = S;
Res.x = X;
Res.y = y;
Res.Dif = y(:) - yfit(:);
Res.Res = r;
end

function [ invM ] = inverta( M )
% [ invM ] = inverta( M )
%
% usar no lugar da função inv

% (c) Zwinglio Guimaraes-Filho (2001-2002)

% para verificar se há elementos na diagonal principal menores do que a
% precisão numérica (eps):
%
if min( sqrt( diag(M) ) ) < eps
M
diag(M)
error('ops')
end

dM = 1./( sqrt( diag(M) ) );
MdM = dM * dM';

invMuni = inv( MdM .* M );
invM = MdM .* invMuni;

end

function resultado = pqui2(chi2,N)

if (chi2<=0 || ~isscalar(chi2))
    disp('Chi2 deve ser estritamente positivo e escalar.');
    return;
end

if (N<1 || ~isscalar(N))
    disp('Número de graus de liberdade N deve ser inteiro e escalar.');
    return;
end

    try
        resultado = 1 - chi2cdf(chi2 , N );
    catch
        F = @(x) ((0.5*x^2)^(0.5*N-1))*exp(-0.5*x^2)/gamma(0.5*N);
        P = 1 - quadl(F, 0, 10*N);
    end
end

function p = chi2cdf(x,v)
%CHI2CDF Chi-square cumulative distribution function.
%   P = CHI2CDF(X,V) returns the chi-square cumulative distribution
%   function with V degrees of freedom at the values in X.
%   The chi-square density function with V degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of P is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.

if   nargin < 2, 
    error('Requires two input arguments.');
end

[errorcode, x, v] = distchck(2,x,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end
    
% Call the gamma distribution function. 
p = gamcdf(x,v/2,2);

% Return NaN if the degrees of freedom is not a positive integer.
k = find(v < 0  |  round(v) ~= v);
if any(k)
   tmp = NaN;
    p(k) = tmp(ones(size(k)));
end

end

function u = u1P(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*J^2)+L*(Y-2);
end

function u = u1M(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*J^2)-L*(Y-2);
end

function u = u3P(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+1)^2)+L*(Y-2);
end

function u = u3M(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+1)^2)-L*(Y-2);
end

function C = C1(J,Y,L)
    C = L^2*Y*(Y-4)*(J-L+1)*(J+L)+2*(2*J+1)*(J-L)*J*(J+L);
end

function C = C2(J,Y,L)
    C = L^2*Y*(Y-4)+4*J*(J+1);
end

function C = C3(J,Y,L)
    C = L^2*Y*(Y-4)*(J+L+1)*(J-L)+2*(2*J+1)*(J-L+1)*(J+1)*(J+L+1);
end


function S = SPTriplet1(J, L, Y1, Y2)

S = (J-L)*(J+L)*((J-L+1)*(J+L-1)*u1P(J-1,Y1,L)*u1P(J,Y2,L)...
    +(J-L-1)*(J+L+1)*u1M(J-1,Y1,L)*u1M(J,Y2,L)...
    +8*(J-L-1)*(J+L-1)*(J-L)*(J+L))^2/...
    (16*J*C1(J-1,Y1,L)*C1(J,Y2,L));

end

function S = SQTriplet1(J, L, Y1, Y2)

S = (2*J+1)*((J-L+1)*(L-1)*(J+L)*u1P(J,Y1,L)*u1P(J,Y2,L)...
    +(J+L+1)*(J-L)*(L+1)*u1M(J,Y1,L)*u1M(J,Y2,L)...
    +8*L*(J-L)^2*(J+L)^2)^2/...
    (16*J*(J+1)*C1(J,Y1,L)*C1(J,Y2,L));

end

function S = SRTriplet1(J, L, Y1, Y2)

S = (J-L+1)*(J+L+1)*((J-L+2)*(J+L)*u1P(J+1,Y1,L)*u1P(J,Y2,L)...
    +(J-L)*(J+L+2)*u1M(J+1,Y1,L)*u1M(J,Y2,L)...
    +8*(J-L+1)*(J+L+1)*(J-L)*(J+L))^2/...
    (16*(J+1)*C1(J+1,Y1,L)*C1(J,Y2,L));

end

function S = SPTriplet2(J, L, Y1, Y2)

S = 4*(J-L)*(J+L)*(0.5*L^2*(Y1-2)*(Y2-2)+(J-L-1)*(J+L+1)...
    +(J-L+1)*(J+L-1))^2/(J*C2(J-1,Y1,L)*C2(J,Y2,L));

end

function S = SQTriplet2(J, L, Y1, Y2)

S = 4*(2*J+1)*(0.5*L^3*(Y1-2)*(Y2-2)+(L+1)*(J-L)*(J+L+1)...
    +(L-1)*(J-L+1)*(J+L))^2/(J*(J+1)*C2(J,Y1,L)*C2(J,Y2,L));

end

function S = SRTriplet2(J, L, Y1, Y2)

S = 4*(J-L+1)*(J+L+1)*(0.5*L^2*(Y1-2)*(Y2-2)+(J-L)*(J+L+2)...
    +(J-L+2)*(J+L))^2/((J+1)*C2(J+1,Y1,L)*C2(J,Y2,L));

end

function S = SPTriplet3(J, L, Y1, Y2)

S = (J-L)*(J+L)*((J-L+1)*(J+L-1)*u3M(J-1,Y1,L)*u3M(J,Y2,L)...
    +(J-L-1)*(J+L+1)*u3P(J-1,Y1,L)*u3P(J,Y2,L)...
    +8*(J-L+1)*(J+L+1)*(J-L)*(J+L))^2/...
    (16*J*C3(J-1,Y1,L)*C3(J,Y2,L));

end

function S = SQTriplet3(J, L, Y1, Y2)

S = (2*J+L)*((L-1)*(J-L+1)*(J+L)*u3M(J,Y1,L)*u3M(J,Y2,L)...
    +(L+1)*(J-L)*(J+L+1)*u3P(J,Y1,L)*u3P(J,Y2,L)...
    +8*L*(J-L+1)^2*(J+L+1))^2/...
    (16*J*(J+1)*C3(J,Y1,L)*C3(J,Y2,L));

end

function S = SRTriplet3(J, L, Y1, Y2)

S = (J-L+1)*(J+L+1)*((J-L+2)*(J+L)*u3M(J+1,Y1,L)*u3M(J,Y2,L)...
    +(J-L)*(J+L+2)*u3P(J+1,Y1,L)*u3P(J,Y2,L)...
    +8*(J-L+1)*(J-L+2)*(J+L+1)*(J+L+2))^2/...
    (16*(J+1)*C3(J+1,Y1,L)*C3(J,Y2,L));

end

function u = uP(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+0.5)^2)+L*(Y-2);
end

function u = uM(J,Y,L)
    u = sqrt(L^2*Y*(Y-4)+4*(J+0.5)^2)-L*(Y-2);
end

function C = CP(J,Y,L)
    C = 0.5*(uP(J,Y,L)^2+4*((J+0.5)^2-L^2));
end

function C = CM(J,Y,L)
    C = 0.5*(uM(J,Y,L)^2+4*((J+0.5)^2-L^2));
end

function S = SPDoublet1(J, L, Y1, Y2)

S = (J-L-0.5)*(J+L+0.5)*(uM(J-1,Y1,L)*uM(J-1,Y2,L)+...
    4*(J-L+0.5)*(J+L-0.5))^2/...
    (4*J*CM(J-1,Y1,L)*CM(J,Y2,L));

end

function S = SPDoublet2(J, L, Y1, Y2)

S = (J-L-0.5)*(J+L+0.5)*(uP(J-1,Y1,L)*uP(J,Y2,L)+...
    4*(J-L+0.5)*(J+L-0.5))^2/(4*J*CM(J-1,Y1,L)*CM(J,Y2,L));

end

function S = SRDoublet1(J, L, Y1, Y2)

S = (J-L+0.5)*(J+L+1.5)*(uM(J+1,Y1,L)*uM(J,Y2,L)+...
    4*(J-L+1.5)*(J+L+0.5))^2/(4*(J+1)*CM(J+1,Y1,L)*CM(J,Y2,L));

end

function S = SRDoublet2(J, L, Y1, Y2)

S = (J-L+0.5)*(J+L+1.5)*(uP(J+1,Y1,L)*uP(J,Y2,L)+...
    4*(J-L+1.5)*(J+L+0.5))^2/(4*(J+1)*CP(J+1,Y1,L)*CP(J,Y2,L));

end

function S = SQDoublet1(J, L, Y1, Y2)

S = (J+0.5)*((L+0.5)*uM(J,Y1,L)*uM(J,Y2,L)...
    +4*(L-0.5)*(J-L+0.5)*(J+L+0.5))^2/...
    (2*J*(J+1)*CM(J,Y1,L)*CM(J,Y2,L));

end

function S = SQDoublet2(J, L, Y1, Y2)

S = (J+0.5)*((L+0.5)*uP(J,Y1,L)*uP(J,Y2,L)...
    +4*(L-0.5)*(J-L+0.5)*(J+L+0.5))^2/...
    (2*J*(J+1)*CP(J,Y1,L)*CP(J,Y2,L));

end

% computation of the secular equation for N2 C3Pi and N2 B3Pi

function [E1, E2, E3] = N2_Energy_e(T,A,B,D,l,g,o,p,q,J)

x = J*(J+1);
H = zeros(3,3);  
H(1,1) = T - A + B*(x+2) + 2*l/3 - 2*g - D*(x^2 + 6*x + 4) - (o+p+q);
H(2,2) = T - 4*l/3 + B*(x+2) - 2*g - 0.5*q*x - D*(x^2+8*x);
H(3,3) = T + A + B*(x-2) + 2*l/3 - D*(x^2-2*x);
H(1,2) = -sqrt(2*x)*(B-0.5*g-0.5*(p+2*q)-2*D*(x+2));
H(2,1) = H(1,2);
H(1,3) = -sqrt(x*(x-2))*(2*D+0.5*q);
H(1,3) = H(3,1);
H(2,3) = -sqrt(2*(x-2))*(B-0.5*g-2*D*x);
H(3,2) = H(2,3);
E = eig(H);
E1 = E(1);
E2 = E(2);
E3 = E(3);
end

% Solve the secular equation for N2 C3Pi and N2 B3Pi states
function [E1, E2, E3] = N2_Energy_f(T,A,B,D,l,g,o,p,q,J)

x = J*(J+1);
H = zeros(3,3);  
H(1,1) = T - A + B*(x+2) + 2*l/3 - 2*g + D*(x^2 + 6*x + 4) - (o+p+q);
H(2,2) = T - 4*l/3 + B*(x+2) - 2*g + 0.5*q*x - D*(x^2+8*x);
H(3,3) = T + A + B*(x-2) + 2*l/3 - D*(x^2-2*x);
H(1,2) = -sqrt(2*x)*(B-0.5*g+0.5*(p+2*q)-2*D*(x+2));
H(2,1) = H(1,2);
H(1,3) = -sqrt(x*(x-2))*(2*D-0.5*q);
H(1,3) = H(3,1);
H(2,3) = -sqrt(2*(x-2))*(B-0.5*g-2*D*x);
H(3,2) = H(2,3);
E = eig(H);
E1 = E(1);
E2 = E(2);
E3 = E(3);
end

