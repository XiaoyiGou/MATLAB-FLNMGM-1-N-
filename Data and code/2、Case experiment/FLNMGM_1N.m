%Program name:Fractional Lag Nonlinear Mixed-frequency Grey Prediction Model£¬abbreviation FLNMGM(1,N)
%r(1¡ªRow2) Lag parameters¡¢r(Row2+1¡ª2*Row2) Power indexs¡¢r(2*Row2+1¡ª2*Row2+2) Nakagami parameters¡¢r(2*Row2+3) Nonlinear trend of adjacent time points


function [result] = FLNMGM_1N(r)

format long g;  

%% Step 1: Data preparation
% % % RawData of Jiangsu EC in 2014-2024
Y0=[5012.54	5114.70	5458.95	5807.89	6128.27	6264.36	6373.71	7101.16	7399.55	7832.96	8487.00];
% % % The raw data of Jiangsu's population, GDP and import and export trade in 2014-2024 and the data needed to predict EC's trend in 2025-2030.
% % % To enable trend forecasting for the period 2025¨C2030, the required future
% % % data of the independent variables have been predicted according to the
% % % method described in the text and stored in the array X0
% % % Due to the differing statistical frequencies of the independent variables, their data lengths vary. 
% % % To ensure consistency, zero-padding is applied at the end of each variable to align their lengths.
X0=[8281.09	8315.11	8381.47	8423.5	8446.19	8469.09	8477.26	8505.4	8515	8526	8526	8577.2941	8599.1127	8620.9868	8642.9165	8664.902	8686.9434	8709.0409	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
12892.85	17222.15	15527.77	19445.55	14620.67	19306.23	17275.47	18914.01	16509.04	20022.69	18749.73	20804.71	18822.6	21998.62	21783.24	23296.48	20440.34	22702.98	23410.33	26653.9	21711.72	24152.19	25111.5	27681.41	21003.74	25634.9	26959.44	29209.6	26036	29281	29911	32164	27824	29004	31717	33544.3	29330	30982	32645	35265.2	31422	33002	34811	37773	33088.6	34327.92566	36921.04504	39676.17379	34907.35346	35774.05629	38912.71528	41747.45423	36762.6552	37309.98326	40879.24097	43800.30397	38649.42418	38881.2304	42848.65678	45855.30793	40550.27651	40467.59381	44827.45516	47919.79584	42462.85977	42064.72328	46817.75934	49996.2439	44385.2969	43670.14018	48818.19914	52083.23478	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
29312339	20100504	28302675	30146352	30899466	28981694	32627712	30466209	29938313	28238359	28046249	29583118	27931663	22596683	25381604	27615801	28896669	28614918	29642631	29142043	31044999	28416065	28622612	31060909	26308573.4	21202307.31	28692452.04	28122974.9	28114087.42	27408522.18	28226231.05	30891916.67	29128722.14	28247736.84	31337836.62	33220889	31074078.19	25637217.13	32138765.88	31082361.64	33019259.44	34672755.96	34199924.1	35727841.47	35819520.36	32637801.89	37064183.19	37259320.76	34739145.27	29815261.94	32725423.92	34293131	36119498	36738439	36155022	38629118	42112414	38808083	40607657	38045767	37065443	27169181	35377427	33159219	36856688	37642361	38957704	37141744	37641076	38014585	38015766	38520013	35613453.93	20357301.47	34026820	36627657	36080845	37849861	40471548	40433561.24	40982248	39490944	40871551	43260629	39802004	33020863	40162688	41787676	41328782.61	43639291.21	42356379	46480317	45978678.25	46295949.5	50264142	51240049.7	47932350.94	36044711.66	44871048.94	37038552	46883698.41	51112543	51694412.29	48944763	47713809.37	43784669.3	45387667.76	45050556.13	40748504	34709110.02	44425822.08	42147888.63	43020341.7	44112760.46	44266489.42	44562034.94	46813163.64	45878344.78	48253484.04	48143323.67	47647415.36	37396353.27	44662307.91	45877305.61	46798252.51	46247484.03	46862887.18	49446441.95	47818251.69	48984919.26	48327155.8	52571468.68	49442177.44	36985523.22	49215756.18	49072132.8	48496989.63	50012389.54	50221874.44	50254672.88	49980508.22	49687024.64	51303850.06	52272996.28	49259180.74	37820880.33	47704431.1	45091665.7	48460153.18	50527197.93	50909768.32	50893959.1	50729124.46	50469478.22	51580298.38	52509172.62	49545407.62	38083048.19	48396306.39	45639249.99	49109469.93	51260581.07	51631711.28	51558430.88	51365810.08	51178201.23	52352139.92	53235770.41	50293273.32	38660064.3	49080922.42	46186540.14	49726314.57	51934224.18	52300518.89	52176229.61	51952174.63	51828134.42	52948260.05	53831138.51	50867062.35	39074272.87	49680874.54	46674853.29	50293216.12	52553359.61	52915542.6	52744555.42	52491646.17	52425597.97	53520694.06	54396419.17	51423123.39	39486648.66	50231591.55	47121493.27	50813324.34	53121143.29	53479738.7	53266059.67	52986941.98	52974475.73	54041091.54	54911913.8	51925929.67	39856739.48	50731930.24	47527359.34	51287383.2	53638548.24	53993980.89	53741379.92	53438503.22	53475494.81	54516907.61	55383018.45];

% % % RawData of Guizhou EC in 2014-2024
% Y0=[1173.74	1174.21	1241.78	1384.89	1482.12	1540.68	1586	1743	1743	1783	1902];
% % % The raw data of Guizhou's population, GDP and import and export trade in 2014-2024 and the data needed to predict EC's trend in 2025-2030.
% % % To enable trend forecasting for the period 2025¨C2030, the required future
% % % data of the independent variables have been predicted according to the
% % % method described in the text and stored in the array X0
% % % Due to the differing statistical frequencies of the independent variables, their data lengths vary. 
% % % To ensure consistency, zero-padding is applied at the end of each variable to align their lengths.
% X0=[3677	3708	3758	3803	3822	3848	3858	3852	3856	3865	3860	3905.9348	3921.2506	3936.6264	3952.0625	3967.5591	3983.1166	3998.735	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0
% 1421.15	2359.16	2389.22	3081.48	1659.39	2691.68	2844.51	3306.98	2102.97	2833.64	3198.77	3599.05	2504.83	3227.52	3767.18	4041.3	3330.22	3717.73	3850.96	4454.3	3690.48	4105.31	4268.36	4705.19	3708.04	4271.33	4647.57	5233.47	4297	4659	4870	5633	4795	4970	4908	5337.4	4922	5258	5112	5621.25	5360	5716	5595	5996.12	5598.45	5950.453126	5821.895585	6344.617445	6015.164729	6221.993713	6043.054799	6594.717625	6358.444227	6502.318653	6265.500184	6828.044517	6678.206701	6777.656312	6486.330816	7056.962734	6983.770827	7044.640781	6702.765623	7281.366953	7279.658737	7303.437884	6913.695303	7500.330342	7567.014916	7554.378899	7118.582649	7713.148392	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0
% 229367	217272	391338	414485	348148	731584	472970	460643	1290782	1260867	585463	249135	266517	183920	261900	249084	1065442	1695851	287721	361753	601114	1100891	853166	991873	655704.58	471200.73	339947.1	514090.91	380715.87	341728.71	353222.06	345031.79	354865.3	320534.1	376449.05	376820	317152.33	210315.39	298676.06	345282.94	433379.58	515887.61	473368.15	573481	621933.11	469731.19	573711.67	704012.53	493679.4	228697.14	430676.43	409405	428974	433438	421774	435399	401722	470957	467278	418638	340430	181653	355536	360097	380891	404273	345129	375601	403305	413643	430362	531703	306501.44	131154.85	311273	381326	384593	405754	607383	642889.98	679144	478691	536157	656733	503170	332278	552559	591910	687693.84	586494.15	441681	445630	509523.5	467163.23	577531	845420.56	422812.94	242206.54	519058.41	478559	546826.79	954512	784575.23	726313	583661.8	638672.71	868624.19	1261434.16	604985	498464.51	792183.73	701496.07	679572.76	603877.22	502619.91	576819.96	437707.42	530339.07	896437.67	1369736.22	667531.71	259535.41	965520.9	604794.49	523042.66	668059.42	598387.05	603270.85	512492.2	560355.32	910561.17	1792261.18	630676.23	499923.46	769090.75	645152.22	552576.44	1662550.375	1647783	1634234.25	1621742.5	1610175	1599420.5	1589385.625	1579991.75	1571170.875	1562865.125	1555024.25	1547604.5	1540567.625	1533880.125	1527512.25	1521437.875	1515633.625	1510078.375	1504753.5	1499642	1494728.875	1490000.25	1485443.75	1481047.875	1476802.25	1472697.375	1468724.625	1464876	1461144.125	1457522.375	1454004.5	1450584.625	1447257.375	1444017.75	1440861.375	1437783.75	1434780.875	1431848.875	1428984.375	1426184.25	1423445.125	1420764.375	1418139.125	1415566.875	1413045.125	1410571.75	1408144.625	1405761.5	1403420.625	1401120.25	1398858.5	1396633.875	1394444.75	1392289.625	1390167.125	1388076	1386014.875	1383982.625	1381978.125	1380000.25	1378047.75	1376119.875	1374215.625	1372334	1370474.125	1368635.125	1366816.25	1365016.625	1363235.625	1361472.5	1359726.5	1357996.875	1356283.125	1354584.625	1352900.625	1351230.625	1349574	1347930.375	1346299.125];


[Row1,Col1]=size(Y0);
[Row2,~]=size(X0);  

%Initial value 
Ys = Y0 ./ Y0(:,1);
Xs = X0 ./ X0(:,1);

x=[1 1 1 1];  %Dependent and independent variables order
bjz=0.5;       %Background value

%Order discrimination
for i=1:Row1+Row2
     if (x(i)==0)||(x(i)==fix(x(i))) 
         x(i) =x(i)+1*10^-10; 
     end
end

Sim_len=Col1-4;  %Simulation interval length
Pre_len=Col1-Sim_len;
% If Steps>0, the initial array X0 must include the corresponding independent variable data for the future interval of the dependent variable to be predicted
Steps=6;     %Future forecast lengt. 
S=[1 4 12];  %Frequency multiple difference between dependent variable and independent variable

%% Step 2: Parameter Estimation
disp(' ')
disp(['¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡¾Parameter Estimation¡¿¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª']);

Ysr_AGO=cumsum(Ys, 2);
Xsr_AGO=cumsum(Xs, 2);

% Construction matrix Y
Y = Ys(2:Sim_len)';


%Construction matrix B
Zr=[];
for i=2:Sim_len;
    tmp=bjz*Ysr_AGO(i)+(1-bjz)*Ysr_AGO(i-1);
    Zr=[Zr,-tmp]; 
end
B1=Zr';


f=[r(2*Row2+1),r(2*Row2+2)];  %Nakagami function parameters
M=[];
N=[];
const = 2 * f(1)^f(1) / (gamma(f(1)) * f(2)^f(1));  % Unified coefficient constant
for i = 1:Row2
    ri = abs(r(i));
    ti = floor(ri);
    t(i) = ti;

    j = 0:ti;
    % Calculate the weight corresponding to item j.
    tmp1 = const * (j.^(2*f(1)-1)) .* exp(-f(1) * j.^2 / f(2));
    tmp2 = sum(tmp1);
    M(i,1:ti+1) = tmp1;  

    if ri == ti  % The lag parameter is an integer
        N(i) = tmp2;
    else  % Lag parameter is non-integer
        tmp3 = const * (ri.^(2*f(1)-1)) * exp(-f(1) * ri.^2 / f(2));
        M(i,ti+2) = tmp3;
        N(i) = tmp2 + tmp3;
    end
end


Beta=[];
for i = 1:Row2
    t(i) = floor(abs(r(i)));
    
    for j = 0:t(i)
        Beta(i, j+1) = M(i, j+1) / N(i);
    end

    % Lag parameter is non-integer
    if r(i) ~= floor(r(i))
        Beta(i, t(i)+2) = M(i, t(i)+2) / N(i);
    end
end


T2=[];
for i=1:Row2
   t(i)=floor(abs(r(i)));
   if r(i)>=0     %The lag parameter is positive.
        for k=2:Sim_len
          tmp6=0;
            for j=0:t(i)
                tmp6=tmp6+Beta(i,j+1)*Xsr_AGO(i,k*S(i)-j);  
            end
            if abs(r(i))==t(i)     %The lag parameter is an integer
                 tmp7=tmp6;  
            else            % Considering the nonlinear development trend between adjacent time points, the corresponding value of non-integer lag period is solved
                 tmp7=tmp6+Beta(i,t(i)+2)*((Xsr_AGO(i,k*S(i)-floor(r(i)))-Xsr_AGO(i,k*S(i)-(floor(r(i))+1)))*power(floor(r(i))+1-r(i),r(2*Row2+3))+Xsr_AGO(i,k*S(i)-(floor(r(i))+1)));   
            end
            T2(i,k-1)=tmp7^r(Row2+i);
       end
   else  %r(i)<0  Lag parameter is negative
        for k=2:Sim_len
          tmp6=0;
            for j=0:t(i)
                tmp6=tmp6+Beta(i,j+1)*Xsr_AGO(i,k*S(i)+j);
            end
            if abs(r(i))==t(i)
                 tmp7=tmp6;
            else
                tmp7=tmp6+Beta(i,t(i)+2)*( (Xsr_AGO(i,k*S(i)-floor(r(i)))-Xsr_AGO(i,k*S(i)-(floor(r(i))+1)))*power(floor(r(i))+1-r(i),r(2*Row2+3))+Xsr_AGO(i,k*S(i)-(floor(r(i))+1)));   
            end
            T2(i,k-1)=tmp7^r(Row2+i); 
        end    
   end
end
B2=T2'; 


T3=[];
T4=[];
for k=2:Sim_len
    T3=[T3,k];
    T4=[T4,1];
end
B3=T3';
B4=T4';
B=[B1,B2,B3,B4]   

% Calculating model parameters:a,b1,b2,...,bn,c,d
fprintf(['Model parameters a,b1,b2,...,bn,c,d:']);
if Sim_len>Row2+4
   Ps=(B'*B)^(-1)*B'*Y
elseif Sim_len==Row2+4
  Ps=B^(-1)*Y
else 
   Ps=B'*(B*B')^(-1)*Y
end    

% Calculate the values of four alternative parameters in the time response formula
len=length(Ps);
miu_1=(1-(1-0.5)*Ps(1))/(1+0.5*Ps(1));  
miu_2N=[];
for i=1:Row2
    miu_2N(i)=Ps(1+i)/(1+0.5*Ps(1));
end
miu_2N;
miu_3=Ps(len-1)/(1+0.5*Ps(1));
miu_4=Ps(len)/(1+0.5*Ps(1));


Simu_Ys_r=[];
Simu_Ys_r(1,1)=Y0(1,1);
 for k=2:Sim_len+Pre_len+Steps
     tmp1=power(miu_1,k-1)*Ys(1,1);
     tmp2=0;
    for w=2:k
      tmp2=tmp2+power(miu_1,k-w)*(miu_3*w+miu_4);
    end  
      tmp6=0;
    for v=2:k
        tmp5=0; 
        for i=1:Row2
            t(i)=floor(abs(r(i)));
            tmp3=0;
            if r(i)>=0
                for j=0:t(i)
                    tmp3=tmp3+Beta(i,j+1)*Xsr_AGO(i,v*S(i)-j);     % Integer lag part
                end
                if abs(r(i))==t(i)
                    tmp4=power(tmp3,r(Row2+i));
                else
                    tmp4_1=tmp3+Beta(i,t(i)+2)*((Xsr_AGO(i,k*S(i)-floor(r(i)))-Xsr_AGO(i,k*S(i)-(floor(r(i))+1)))*power(floor(r(i))+1-r(i),r(2*Row2+3))+Xsr_AGO(i,k*S(i)-(floor(r(i))+1)));
                    tmp4=power(tmp4_1,r(Row2+i));
                end
                tmp5=tmp5+tmp4*miu_2N(i);
            else %r(i)<0
                for j=0:t(i)
                    tmp3=tmp3+Beta(i,j+1)*Xsr_AGO(i,v*S(i)+j);   % Integer lag part
                end
                if abs(r(i))==t(i)
                     tmp4=power(tmp3,r(Row2+i));
                else   
                     tmp4_1=tmp3+Beta(i,t(i)+2)*((Xsr_AGO(i,k*S(i)-floor(r(i)))-Xsr_AGO(i,k*S(i)-(floor(r(i))+1)))*power(floor(r(i))+1-r(i),r(2*Row2+3))+Xsr_AGO(i,k*S(i)-(floor(r(i))+1)));   
                     tmp4=power(tmp4_1,r(Row2+i));
                end   
                tmp5=tmp5+tmp4*miu_2N(i);
            end
       end
       tmp6=tmp6+power(miu_1,k-v)*tmp5;
    end
        tmp=tmp1+tmp2+tmp6;
        Simu_Ys_r=[Simu_Ys_r,tmp*Y0(1,1)];  %Initial value reduction
end

% All_cal_data: Used to store data calculated by the model.
% Include simulation data, prediction test data and future prediction data.
All_cal_data = [Simu_Ys_r(1), diff(Simu_Ys_r)];

disp(' ')
fprintf(['¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡¾Calculate simulation value, residual error and relative error¡¿¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª']);
SNO=[];
Raw_data=[];
Sim_data=[];
Residual_sim_error=[];
Relative_sim_error=[];
Sum_sim_error=0;
 for k=2:Sim_len
    canCha=All_cal_data(1,k)-Y0(1,k);
    xdError=canCha/Y0(1,k);
    Sum_sim_error=Sum_sim_error+abs(xdError*100);
    SNO=[SNO,k];
    Raw_data=[Raw_data,Y0(1,k)];
    Sim_data=[Sim_data,All_cal_data(1,k)];
    Residual_sim_error=[Residual_sim_error,canCha];
    Relative_sim_error=[Relative_sim_error,xdError];
 end

% Display the prediction results in a table
Col_name={'No','Raw_data','Simulated_data','Residual_error','Percentage_error'};
Simulation=table(SNO',Raw_data',Sim_data',Residual_sim_error',Relative_sim_error','VariableNames',Col_name)
Mean_sim_error=Sum_sim_error/(Sim_len-1);
str=['MAPE_Train =' num2str(Mean_sim_error) '%'];
disp(str);


Sum_pre_error=0;
if Pre_len>0
    disp(' ')
    fprintf(['¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡¾Calculate prediction value, residual error and relative error¡¿¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª']);
    SNO_pre=[];
    Pre_data=[];
    Raw_data=[];
    Residual_pre_error=[];
    Relative_pre_error=[];
  for k=Sim_len+1:Sim_len+Pre_len
    canCha=All_cal_data(1,k)-Y0(1,k);
    xdError=canCha/Y0(1,k);
    Sum_pre_error=Sum_pre_error+abs(xdError*100);
    SNO_pre=[SNO_pre,k];
    Raw_data=[Raw_data,Y0(1,k)];
    Pre_data=[Pre_data,All_cal_data(1,k)];
    Residual_pre_error=[Residual_pre_error,canCha];
    Relative_pre_error=[Relative_pre_error,xdError];
  end

  %Display the prediction results in a table
  Col_name={'No','Raw_data','Predicted_data','Residual_error','Percentage_error'};
  Prediction=table(SNO_pre',Raw_data',Pre_data',Residual_pre_error',Relative_pre_error','VariableNames',Col_name)
  Mean_pre_error=Sum_pre_error/Pre_len;
  str=['MAPE_Test =' num2str(Mean_pre_error) '%'];
  disp(str);    
end

disp(['¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª']);
MAPE=(Sum_sim_error+Sum_pre_error)/(Sim_len+Pre_len-1);
str2=['MAPE_All =' num2str(MAPE) '%'];
disp(str2);    
Relative_all = [Relative_sim_error(:); Relative_pre_error(:)] ;
STD=sqrt(mean((Relative_all(:) - MAPE/100).^2));
str3=['STD_All =' num2str(STD) ];
disp(str3); 

if Steps>0
   disp(' ')
   fprintf(['¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡¾Applying FLNMGM(1,N) model to forecast the future trend¡¿¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª']);
   SNO=[];
   Pre_steps_data=[];
  for k=Sim_len+Pre_len+1:Sim_len+Pre_len+Steps
     pre=All_cal_data(1,k);
     SNO=[SNO,k];
     Pre_steps_data=[Pre_steps_data,pre];
  end

 %Display the prediction results in a table
Col_name={'No','Predicted_data'};
Prediction_steps=table(SNO',Pre_steps_data','VariableNames',Col_name)    
end
    result = MAPE;
end