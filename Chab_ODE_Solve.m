function [state_var]=Chab_ODE_Solve(x,ideal_load_cycle_data,Rel_Tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to solve the Chaboche ODEs
%
%James Rouse - 20/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note:
%The following convention is used for the state_var array
%state_var(:,1,:)=plastic (inelastic) strain component (%)
%state_var(:,2,:)=accumulated plastic (inelastic) strain component (%)
%state_var(:,3,:)=Drag Stress (MPa)
%state_var(:,4,:)=Stress (MPa)
%state_var(:,5,:)=Back Stress Component 1 (MPa)
%state_var(:,6,:)=Back Stress Component 2 (MPa)
%state_var(:,7,:)=ev1 strain
%state_var(:,8,:)=ev2 strain
%state_var(:,9,:)=ev3 strain

%Define options set for ODE45
options=odeset('RelTol',Rel_Tol);

state_var=cell(size(ideal_load_cycle_data,1),1);

for jj=1:1:size(ideal_load_cycle_data,1)
    ideal_load_cycle_data_temp=cell2mat(ideal_load_cycle_data(jj,1));
    if jj==1
        state_var_1=zeros(size(ideal_load_cycle_data_temp,1),9);
        state_var_1(1,5)=0;
    else
        state_var_1=zeros(size(ideal_load_cycle_data_temp,1),9);
        state_var_temp=cell2mat(state_var(jj-1,1));
        state_var_1(1,:)=state_var_temp(end,:);
        clear state_var_temp
    end
    
    for ii=2:1:size(ideal_load_cycle_data_temp,1)
        t0=ideal_load_cycle_data_temp(ii-1,1);
        t_final=ideal_load_cycle_data_temp(ii,1);
        strain_rate=(ideal_load_cycle_data_temp(ii,2)-ideal_load_cycle_data_temp(ii-1,2))/(100*(t_final-t0));
        y0=[(state_var_1(ii-1,1)/100),(state_var_1(ii-1,2)/100),state_var_1(ii-1,4),state_var_1(ii-1,5),state_var_1(ii-1,6),(state_var_1(ii-1,7)/100),(state_var_1(ii-1,8)/100),(state_var_1(ii-1,9)/100)];
        [t,y]=ode45('Chab_ODEs',[t0,t_final],y0,options,x,strain_rate); %Call ODE45 to solve Chaboche ODEs
        state_var_1(ii,1)=y(end,1)*100;%Plastic Strain
        state_var_1(ii,2)=y(end,2)*100;%Accumulated Plastic Strain
        state_var_1(ii,3)=((x(9)*(1-exp((-1)*x(8)*y(end,2))))+(x(11)*y(end,2)));%Drag Stress
        state_var_1(ii,4)=y(end,3);%Stress
        state_var_1(ii,5)=y(end,4);%Back Stress Component 1
        state_var_1(ii,6)=y(end,5);%Back Stress Component 2
        state_var_1(ii,7)=y(end,6)*100;%ev1
        state_var_1(ii,8)=y(end,7)*100;%ev1
        state_var_1(ii,9)=y(end,8)*100;%ev1
        clear t y t0 t_final strain_rate y0
    end

    state_var{jj,1}=state_var_1;
    clear ideal_load_cycle_data_temp state_var_1
    jj
end
