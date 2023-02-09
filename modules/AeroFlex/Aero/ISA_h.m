%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 - 2011 
% 
% Sergio Ricci (sergio.ricci@polimi.it)
%
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale
% Via La Masa 34, 20156 Milano - ITALY
% 
% This file is part of NeoCASS Software (www.neocass.org)
%
% NeoCASS is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2, or (at your option) any later version.
%
% NeoCASS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public
% License along with NeoCASS; see the file GNU GENERAL 
% PUBLIC LICENSE.TXT.  If not, write to the Free Software 
% Foundation, 59 Temple Place -Suite 330, Boston, MA
% 02111-1307, USA.
%
%
%
%*******************************************************************************
%  SimSAC Project
%
%  NeoCASS
%  Next generation Conceptual Aero Structural Sizing  
%
%                      Sergio Ricci             <ricci@aero.polimi.it>
%                      Luca Cavagna             <cavagna@aero.polimi.it>
%                      Alessandro De Gaspari    <degaspari@aero.polimi.it>
%                      Luca Riccobene           <riccobene@aero.polimi.it>
%
%  Department of Aerospace Engineering - Politecnico di Milano (DIAPM)
%  Warning: This code is released only to be used by SimSAC partners.
%  Any usage without an explicit authorization may be persecuted.
%
%*******************************************************************************
% FUNCTION           [rho, p, T, a, mu] = ISA_h(altitude)
%
%*******************************************************************************
function [rho, p, T, a, mu] = ISA_h(altitude)

% Sea Level Conditions
rho_ssl = 1.225;         % Kg/m^3
p_ssl   = 101300;        % Pa
T_ssl   = 15;            % C
a_ssl   = 340.3;         % m/s
mu_ssl  = 1.789*10^-5;   % Kg/m/s
%
nu_ssl  = 1.460*10^-5;   % m^2/s

% Check input
if altitude < 0
  rho = 1.225;
  p   = P_ssl;
  T   = T_ssl + 273.15; 
  a   = a_ssl;
  mu  = mu_ssl;
  return
end

if altitude >= 21336
  error('Altitude required greater than 21336 m. No data available from ISA.')
end

% Atmosperic Table
%----------------------------------------------------------------------
%Altitude  Temperature                                Kinematic   Speed
%  m      ft     oC    Pressure   Density   Viscosity Viscosity     of
%                      Ratio      Ratio      Ratio    Ratio       Sound
%-----------------------------------------------------------------------
A=[ 0     0    15.2   1.0000     1.0000     1.0000    1.0000    340.3 
  152   500    14.2   0.9821     0.9855     0.9973    1.0121    339.7 
  304  1000    13.2   0.9644     0.9711     0.9947    1.0243    339.1 
  457  1500    12.2   0.9470     0.9568     0.9920    1.0367    338.5 
  609  2000    11.2   0.9298     0.9428     0.9893    1.0493    338.0 
  762  2500    10.2   0.9129     0.9289     0.9866    1.0622    337.4 
  914  3000     9.3   0.8962     0.9151     0.9839    1.0752    336.8 
 1066  3500     8.3   0.8798     0.9015     0.9812    1.0884    336.2 
 1219  4000     7.3   0.8637     0.8881     0.9785    1.1018    335.6 
 1371  4500     6.3   0.8477     0.8748     0.9758    1.1155    335.0 
 1524  5000     5.3   0.8320     0.8617     0.9731    1.1293    334.4      
 1676  5500     4.3   0.8166     0.8487     0.9704    1.1434    333.8 
 1828  6000     3.3   0.8014     0.8359     0.9677    1.1577    333.2 
 1981  6500     2.3   0.7864     0.8232     0.9649    1.1722    332.6 
 2133  7000     1.3   0.7716     0.8106     0.9622    1.1870    332.0 
 2286  7500     0.3   0.7571     0.7983     0.9595    1.2020    331.4 
 2438  8000    -0.6   0.7428     0.7860     0.9567    1.2172    330.8 
 2590  8500    -1.6   0.7287     0.7739     0.9540    1.2327    330.2 
 2743  9000    -2.6   0.7148     0.7620     0.9512    1.2484    329.6 
 2895  9500    -3.6   0.7012     0.7501     0.9485    1.2644    329.0 
 3048 10000    -4.6   0.6877     0.7385     0.9457    1.2807    328.4       
 3200 10500    -5.6   0.6745     0.7269     0.9430    1.2972    327.8 
 3352 11000    -6.6   0.6614     0.7155     0.9402    1.3140    327.2 
 3505 11500    -7.6   0.6486     0.7043     0.9374    1.3310    326.6 
 3657 12000    -8.6   0.6360     0.6932     0.9347    1.3484    326.0 
 3810 12500    -9.6   0.6236     0.6822     0.9319    1.3660    325.4 
 3962 13000   -10.6   0.6113     0.6713     0.9291    1.3840    324.7  
 4114 13500   -11.5   0.5993     0.6606     0.9263    1.4022    324.1 
 4267 14000   -12.5   0.5875     0.6500     0.9235    1.4207    323.5 
 4419 14500   -13.5   0.5758     0.6396     0.9207    1.4396    322.9 
 4572 15000   -14.5   0.5643     0.6292     0.9179    1.4588    322.3       
 4724 15500   -15.5   0.5531     0.6190     0.9151    1.4783    321.7 
 4876 16000   -16.5   0.5420     0.6090     0.9123    1.4981    321.0 
 5029 16500   -17.5   0.5311     0.5990     0.9094    1.5183    320.4 
 5181 17000   -18.5   0.5203     0.5892     0.9066    1.5388    319.8  
 5334 17500   -19.5   0.5098     0.5795     0.9038    1.5596    319.2  
 5486 18000   -20.5   0.4994     0.5699     0.9009    1.5809    318.5 
 5638 18500   -21.5   0.4892     0.5604     0.8981    1.6025    317.9 
 5791 19000   -22.4   0.4791     0.5511     0.8953    1.6244    317.3 
 5943 19500   -23.4   0.4693     0.5419     0.8924    1.6468    316.7 
 6096 20000   -24.4   0.4595     0.5328     0.8895    1.6696    316.0       
 6248 20500   -25.4   0.4500     0.5238     0.8867    1.6927    315.4 
 6400 21000   -26.4   0.4406     0.5150     0.8838    1.7163    314.8 
 6553 21500   -27.4   0.4314     0.5062     0.8809    1.7403    314.1 
 6705 22000   -28.4   0.4223     0.4976     0.8781    1.7647    313.5 
 6858 22500   -29.4   0.4134     0.4891     0.8752    1.7895    312.9 
 7010 23000   -30.4   0.4046     0.4806     0.8723    1.8148    312.2 
 7162 23500   -31.4   0.3960     0.4723     0.8694    1.8406    311.6  
 7315 24000   -32.3   0.3876     0.4642     0.8665    1.8668    311.0 
 7467 24500   -33.3   0.3793     0.4561     0.8636    1.8935    310.3 
 7620 25000   -34.3   0.3711     0.4481     0.8607    1.9207    309.7       
 7772 25500   -35.3   0.3631     0.4402     0.8578    1.9484    309.0 
 7924 26000   -36.3   0.3552     0.4325     0.8548    1.9766    308.4 
 8077 26500   -37.3   0.3474     0.4248     0.8519    2.0053    307.7 
 8229 27000   -38.3   0.3398     0.4173     0.8490    2.0345    307.1 
 8382 27500   -39.3   0.3324     0.4098     0.8460    2.0643    306.4 
 8534 28000   -40.3   0.3250     0.4025     0.8431    2.0947    305.8 
 8686 28500   -41.3   0.3178     0.3953     0.8402    2.1256    305.1 
 8839 29000   -42.3   0.3107     0.3881     0.8372    2.1571    304.5 
 8991 29500   -43.2   0.3038     0.3811     0.8342    2.1892    303.8 
 9144 30000   -44.2   0.2970     0.3741     0.8313    2.2219    303.2       
 9296 30500   -45.2   0.2903     0.3673     0.8283    2.2553    302.5 
 9448 31000   -46.2   0.2837     0.3605     0.8253    2.2892    301.9
 9601 31500   -47.2   0.2772     0.3539     0.8223    2.3239    301.2
 9753 32000   -48.2   0.2709     0.3473     0.8194    2.3592    300.5
 9906 32500   -49.2   0.2647     0.3408     0.8164    2.3952    299.9
10058 33000   -50.2   0.2586     0.3345     0.8134    2.4318    299.2 
10210 33500   -51.2   0.2526     0.3282     0.8104    2.4692    298.6 
10363 34000   -52.2   0.2467     0.3220     0.8073    2.5074    297.9 
10515 34500   -53.2   0.2410     0.3159     0.8043    2.5463    297.2 
10668 35000   -54.1   0.2353     0.3099     0.8013    2.5859    296.5       
10820 35500   -55.1   0.2298     0.3039     0.7983    2.6264    295.9 
10972 36000   -56.1   0.2243     0.2981     0.7952    2.6677    295.2 
10999 36089   -56.3   0.2234     0.2971     0.7947    2.6751    295.1 
11277 37000   -56.3   0.2138     0.2843     0.7947    2.7948    295.1 
11582 38000   -56.3   0.2038     0.2710     0.7947    2.9324    295.1 
11887 39000   -56.3   0.1942     0.2583     0.7947    3.0768    295.1 
12192 40000   -56.3   0.1851     0.2462     0.7947    3.2283    295.1        
12496 41000   -56.3   0.1764     0.2346     0.7947    3.3872    295.1 
12801 42000   -56.3   0.1681     0.2236     0.7947    3.5540    295.1 
13106 43000   -56.3   0.1602     0.2131     0.7947    3.7290    295.1 
13411 44000   -56.3   0.1527     0.2031     0.7947    3.9126    295.1 
13716 45000   -56.3   0.1456     0.1936     0.7947    4.1052    295.1       
14020 46000   -56.3   0.1387     0.1845     0.7947    4.3073    295.1 
14325 47000   -56.3   0.1322     0.1758     0.7947    4.5194    295.1 
14630 48000   -56.3   0.1260     0.1676     0.7947    4.7419    295.1 
14935 49000   -56.3   0.1201     0.1597     0.7947    4.9754    295.1 
15240 50000   -56.3   0.1145     0.1522     0.7947    5.2203    295.1      
15544 51000   -56.3   0.1091     0.1451     0.7947    5.4773    295.1 
15849 52000   -56.3   0.1040     0.1383     0.7947    5.7470    295.1 
16154 53000   -56.3   0.09909    0.1318     0.7947    6.0300    295.1 
16459 54000   -56.3   0.09444    0.1256     0.7947    6.3268    295.1 
16764 55000   -56.3   0.09001    0.1197     0.7947    6.6383    295.1     
17068 56000   -56.3   0.08579    0.1141     0.7947    6.9652    295.1 
17373 57000   -56.3   0.08176    0.1087     0.7947    7.3081    295.1 
17678 58000   -56.3   0.07793    0.1036     0.7947    7.6679    295.1 
17983 59000   -56.3   0.07427    0.09878    0.7947    8.0454    295.1 
18288 60000   -56.3   0.07079    0.09414    0.7947    8.4416    295.1       
18592 61000   -56.3   0.06746    0.08972    0.7947    8.8572    295.1 
18897 62000   -56.3   0.06430    0.08551    0.7947    9.2932    295.1 
19202 63000   -56.3   0.06128    0.08150    0.7947    9.7508    295.1 
19507 64000   -56.3   0.05841    0.07768    0.7947    10.231    295.1 
19812 65000   -56.3   0.05566    0.07403    0.7947    10.735    295.1 
20116 66000   -56.3   0.05305    0.07056    0.7947    11.263    295.1 
20421 67000   -56.3   0.05056    0.06725    0.7947    11.818    295.1 
20726 68000   -56.3   0.04819    0.06409    0.7947    12.399    295.1 
21031 69000   -56.3   0.04593    0.06108    0.7947    13.010    295.1 
21336 70000   -56.3   0.04377    0.05822    0.7947    13.650    295.1 ];

% density
rho = interp1(A(:,1), A(:,5), altitude);
rho = rho .* rho_ssl;
% pressure
p = interp1(A(:,1), A(:,4), altitude);
p = p .* p_ssl;
% Temperature
T = interp1(A(:,1), A(:,3), altitude);
T = T + 273.15;
% a
a = interp1(A(:,1), A(:,8), altitude);
% mu
mu = interp1(A(:,1), A(:,6), altitude);
mu = mu .* mu_ssl;
