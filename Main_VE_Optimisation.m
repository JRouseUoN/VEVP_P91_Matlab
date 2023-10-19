clear all
close all
clc
plotting_cols_master=[1,0,0;
                      0,0,1;
                      0,1,0;
                      0,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify isothermal data to be used in viscoelastic parameter estiamtion
%and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VE_perlim_400.mat')
%
Select_data_sets=cell(3,1);
Select_data_sets{1,1}='int_dwell_0_1_400.mat';
Select_data_sets{2,1}='int_dwell_0_01_400.mat';
Select_data_sets{3,1}='int_dwell_0_001_400.mat';
%
Selected_strain_rate=zeros(3,1);
Selected_strain_rate(1,1)=0.1;
Selected_strain_rate(2,1)=0.01;
Selected_strain_rate(3,1)=0.001;
%
NAMES=cell(6,1);
NAMES{1,1}='Exp. 0.1%/s';
NAMES{2,1}='Exp. 0.01%/s';
NAMES{3,1}='Exp. 0.001%/s';
NAMES{4,1}='Pred. 0.1%/s';
NAMES{5,1}='Pred. 0.01%/s';
NAMES{6,1}='Pred. 0.001%/s';
%
cut_off_time=18000;
cut_off_strain=0.19;
%
%E=1.2254e11; %600C
%E=1.4619e11; %500C
E=1.680e11; %400C
%
E=E*1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create ideal strain profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ideal_strain_profiles=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    [ideal_load_cycle_data_1]=ideal_load_cycle_gen_1(Selected_strain_rate(ii,1)); %"Special", interupted dwell, cycle
    %[ideal_load_cycle_data_2]=ideal_load_cycle_gen_2(strain_rate); %Saw tooth
    %
    ideal_load_cycle_data_1=ideal_load_cycle_data_1(find(ideal_load_cycle_data_1(:,1)<=18000),:);
    %
    ideal_load_temp=cell(1,1);
    ideal_load_temp{1,1}=ideal_load_cycle_data_1;
    ideal_strain_profiles{ii,1}=ideal_load_temp;
    clear ideal_load_cycle_data_1 ideal_load_temp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create stress check profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stress_check=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{1,1};
    data_temp(1,1)=0;
    ideal_load_temp=cell2mat(ideal_strain_profiles{ii,1});
    %
    stress_check_temp=cell(1,1);
    stress_check_temp{1,1}=interp1(data_temp(:,1),data_temp(:,3),ideal_load_temp(:,1));
    %
    stress_check{ii,1}=stress_check_temp;
    clear data data_proc data_temp ideal_load_temp stress_check_temp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve ideal strain profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x(1)=C1;
%x(2)=gamma1;
%x(3)=C2;
%x(4)=gamma2;
%x(5)=A;
%x(6)=B;
%x(7)=D;
%x(8)=b;
%x(9)=Q;
%x(10)=k;
%x(11)=H;
% x(12)=E0;
% x(13)=E1;
% x(14)=eta1;
% x(15)=E2;
% x(16)=eta2;
% x(17)=E3;
% x(18)=eta3;
%
%Initial Conditions
x0=zeros(7,1);
x0(1)=E/1e6;
x0(2)=E_VE(1);
x0(3)=eta_VE(1);
x0(4)=E_VE(2);
x0(5)=eta_VE(2);
x0(6)=E_VE(3);
x0(7)=eta_VE(3);
%
x_scales=[1e5;1e5;1e10;1e5;1e9;1e5;1e9];
x0=x0./x_scales;
%
lb=x0.*1e-3;
ub=x0.*1e3;
%
%Define differential equation solver parameters
Rel_Tol=1e-4;
% [x]=Optim_VE(x0,x_scales,lb,ub,ideal_strain_profiles{1,1},ideal_strain_profiles{2,1},ideal_strain_profiles{3,1},...
%     stress_check{1,1},stress_check{2,1},stress_check{3,1},Rel_Tol);
[x]=Optim_VE(x0,x_scales,lb,ub,ideal_strain_profiles{3,1},ideal_strain_profiles{2,1},...
    stress_check{3,1},stress_check{2,1},Rel_Tol);
%
x=x.*x_scales;
x_op=zeros(18,1);
x_op(7)=1;
x_op(10)=1000000;
x_op(12)=x(1);
x_op(13)=x(2);
x_op(14)=x(3);
x_op(15)=x(4);
x_op(16)=x(5);
x_op(17)=x(6);
x_op(18)=x(7);
%
Results=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    %Solve Chaboche equations for each point in ideal_load_cycle_data
    Results{ii,1}=Chab_ODE_Solve(x_op,ideal_strain_profiles{ii,1},Rel_Tol);      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
Select_cycle=1;
%
figure = figure('Color',[1 1 1]);
title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{Select_cycle,1};
    plot(data_temp(:,2),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
    load_temp=cell2mat(ideal_strain_profiles{ii,1});
    results_temp=cell2mat(Results{ii,1});
    if ii==1
        plot(load_temp(:,2),results_temp(:,4), 'k-','MarkerSize', 2, 'LineWidth', 2);
    elseif ii==2
        plot(load_temp(:,2),results_temp(:,4), 'k--','MarkerSize', 2, 'LineWidth', 2);
    elseif ii==3
        plot(load_temp(:,2),results_temp(:,4), 'k:','MarkerSize', 2, 'LineWidth', 2);
    end
    clear load_temp results_temp
%     if ii==1
%         max_strain=max(data_temp(:,2));
%         max_stress=max(data_temp(:,3));
%     else
%         max_strain=max([max_strain;max(data_temp(:,2))]);
%         max_stress=max([max_stress;max(data_temp(:,3))]);
%     end
    clear data_temp
end
%
xlabel('\epsilon (%)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
xlim([-0.6 0.6])
ylim([-500 500])
hold off
clear figure  plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
Select_cycle=1;
%
figure = figure('Color',[1 1 1]);
title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{Select_cycle,1};
    plot(data_temp(:,1),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
    load_temp=cell2mat(ideal_strain_profiles{ii,1});
    results_temp=cell2mat(Results{ii,1});
    if ii==1
        plot(load_temp(:,1),results_temp(:,4), 'k-','MarkerSize', 2, 'LineWidth', 2);
    elseif ii==2
        plot(load_temp(:,1),results_temp(:,4), 'k--','MarkerSize', 2, 'LineWidth', 2);
    elseif ii==3
        plot(load_temp(:,1),results_temp(:,4), 'k:','MarkerSize', 2, 'LineWidth', 2);
    end
    clear load_temp results_temp
%     if ii==1
%         max_strain=max(data_temp(:,2));
%         max_stress=max(data_temp(:,3));
%     else
%         max_strain=max([max_strain;max(data_temp(:,2))]);
%         max_stress=max([max_stress;max(data_temp(:,3))]);
%     end
    clear data_temp
end
%
xlabel('t (s)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
ylim([-500 500])
hold off
clear figure plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%