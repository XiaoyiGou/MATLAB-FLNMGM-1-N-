
function [result] = FLNMGM_1N(r,Y0,X1,X2)

% clc,clear all;
format long g;  

%% Step 1: Data input
X0 = [X1; X2];
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

Sim_len=Col1-3;  % 90% 3;85% 5; 75% 7
Pre_len=Col1-Sim_len;
Steps=0;
S=[1 4];

%% Step 2: Parameter Estimation

Ysr_AGO=cumsum(Ys, 2);
Xsr_AGO=cumsum(Xs, 2);

% Construction matrix Y
Y = Ys(2:Sim_len)';

%Construction matrix B
Zr=[];
for i=2:Sim_len;
    tmp=bjz*Ysr_AGO(i)+(1-bjz)*Ysr_AGO(i-1);
    Zr=[Zr,-tmp]; 
end;
Zr;
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
B=[B1,B2,B3,B4];

% Calculating model parameters:a,b1,b2,...,bn,c,d
fprintf(['Model parameters a,b1,b2,...,bn,c,d:']);
if Sim_len>Row2+4
   Ps=(B'*B)^(-1)*B'*Y;
elseif Sim_len==Row2+4
  Ps=B^(-1)*Y;
else 
   Ps=B'*(B*B')^(-1)*Y;
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

%% Step 3: Calculate simulated value and predicted value
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
                    tmp3=tmp3+Beta(i,j+1)*Xsr_AGO(i,v*S(i)-j);
                end
                if abs(r(i))==t(i)
                    tmp4=tmp3;
                else   
                    tmp4_1=tmp3+Beta(i,t(i)+2)*((Xsr_AGO(i,k*S(i)-floor(r(i)))-Xsr_AGO(i,k*S(i)-(floor(r(i))+1)))*power(floor(r(i))+1-r(i),r(2*Row2+3))+Xsr_AGO(i,k*S(i)-(floor(r(i))+1)));   
                    tmp4=power(tmp4_1,r(Row2+i));
                end
                    tmp5=tmp5+tmp4*miu_2N(i);
            else %r(i)<0
                for j=0:t(i)
                    tmp3=tmp3+Beta(i,j+1)*Xsr_AGO(i,v*S(i)+j);
                end
                if abs(r(i))==t(i)
                     tmp4=tmp3;
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
        Simu_Ys_r=[Simu_Ys_r,tmp*Y0(1,1)];
end

%% Step 4:Decremental reduction
% All_cal_data: Used to store data calculated by the model.
% Include simulation data, prediction test data and future prediction data.
All_cal_data = [Simu_Ys_r(1,1), diff(Simu_Ys_r(1,:))];


%% Step 5:Calculate residuals and relative errors, and display the modeling results.
disp(' ')
fprintf(['！！！！！！！！！！！！！！！！‐Calculate simulation value, residual error and relative error／！！！！！！！！！！！！！！！！！']);
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
    fprintf(['！！！！！！！！！！！！！！！！‐Calculate prediction value, residual error and relative error／！！！！！！！！！！！！！！！！！']);
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

disp(['！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！']);
MAPE=(Sum_sim_error+Sum_pre_error)/(Sim_len+Pre_len-1);
str2=['MAPE_All =' num2str(MAPE) '%'];
disp(str2);    
% Relative_all = [Relative_sim_error(:); Relative_pre_error(:)] ;
% STD=sqrt(mean((Relative_all(:) - MAPE/100).^2));
% str3=['STD_All =' num2str(STD) ];
% disp(str3);    

result = MAPE;
end