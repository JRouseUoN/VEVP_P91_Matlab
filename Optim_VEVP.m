function [x]=Optim_VEVP(x0,x_scales,lb,ub,ideal_load_cycle_data_1,ideal_load_cycle_data_2,ideal_load_cycle_data_3,...
    sigm_check_1,sigm_check_2,sigm_check_3,Rel_Tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trapezoid Function Optimisation (Fitting)
%
%James Rouse - 6/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Define virtual function for lsqnonlin and set a handle for the function
format long
fun=@(x)objective_fun_JPR(x,x_scales,ideal_load_cycle_data_1,ideal_load_cycle_data_2,ideal_load_cycle_data_3,...
    sigm_check_1,sigm_check_2,sigm_check_3,Rel_Tol);

Tol_Fun=10^-5;%Difference limit on function values between evaluations  
Tol_X=10^-5;%Difference limit on optimisation variable values between evaluations 
Max_Fun_Evals=10^6;%Maximum number of evaluations  
Max_Iter=400;

%Define optimisation options in least square method 
options=optimset('display','on','TolFun',Tol_Fun,'TolX',Tol_X,'MaxFunEvals',Max_Fun_Evals,'MaxIter',Max_Iter);

%Call optimisation program (least squares method) 
[x,resnorm,residual]=lsqnonlin(fun,x0,lb,ub,options);
%[x,resnorm,residual]=fminunc(fun,x0,options);
 
res=residual;%Residual of objective functions 

squr=resnorm;%Square of residual of objective functions 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=objective_fun_JPR(x,x_scales,ideal_load_cycle_data_1,ideal_load_cycle_data_2,ideal_load_cycle_data_3,...
    sigm_check_1,sigm_check_2,sigm_check_3,Rel_Tol);


% ideal_load_cycle_data_1=ideal_strain_profiles{1,1};
% ideal_load_cycle_data_2=ideal_strain_profiles{2,1};
% ideal_load_cycle_data_3=ideal_strain_profiles{3,1};
% sigm_check_1=stress_check{1,1};
% sigm_check_2=stress_check{2,1};
% sigm_check_3=stress_check{3,1};
%
xf=x.*x_scales;
xf(8)=0;
xf(9)=0;
xf(11)=0;
%
%
[state_var_1]=Chab_ODE_Solve(xf,ideal_load_cycle_data_1,Rel_Tol); 
[state_var_2]=Chab_ODE_Solve(xf,ideal_load_cycle_data_2,Rel_Tol); 
[state_var_3]=Chab_ODE_Solve(xf,ideal_load_cycle_data_3,Rel_Tol); 

for ii=1:1:size(sigm_check_1,1)
    state_var_temp=cell2mat(state_var_1(ii,1));
    sigm_check_temp=cell2mat(sigm_check_1(ii,1));
    F_temp=sigm_check_temp(:,1)-state_var_temp(1:length(state_var_temp),4);
    if ii==1
        F_1=F_temp;
    else
        F_1=[F_1;F_temp];
    end    
    clear state_var_temp sigm_check_temp F_temp
end
clear F_temp
for ii=1:1:size(sigm_check_2,1)
    state_var_temp=cell2mat(state_var_2(ii,1));
    sigm_check_temp=cell2mat(sigm_check_2(ii,1));
    F_temp=sigm_check_temp(:,1)-state_var_temp(1:length(state_var_temp),4);
    if ii==1
        F_2=F_temp;
    else
        F_2=[F_2;F_temp];
    end    
    clear state_var_temp sigm_check_temp F_temp
end
clear F_temp
for ii=1:1:size(sigm_check_3,1)
    state_var_temp=cell2mat(state_var_3(ii,1));
    sigm_check_temp=cell2mat(sigm_check_3(ii,1));
    F_temp=sigm_check_temp(:,1)-state_var_temp(1:length(state_var_temp),4);
    if ii==1
        F_3=F_temp;
    else
        F_3=[F_3;F_temp];
    end    
    clear state_var_temp sigm_check_temp F_temp
end
clear F_temp

F=[F_1;F_2;F_3];

for ii=1:1:size(F,1)
    if isnan(F(ii,1))==1
        F(ii,1)=1000;
    end
    if isinf(F(ii,1))==1
        F(ii,1)=1000;
    end
end

end
