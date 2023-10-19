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
load('VEVP_perlim_500.mat')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_in(12)=1.70e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create ideal strain profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ideal_strain_profiles=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    [ideal_load_cycle_data_1]=ideal_load_cycle_gen_1(Selected_strain_rate(ii,1)); %"Special", interupted dwell, cycle
    %[ideal_load_cycle_data_2]=ideal_load_cycle_gen_2(strain_rate); %Saw tooth
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

%Initial Conditions
x0=X_in;
x_scales=10.^floor(log10(x0));
x0=x0./x_scales;
lb=x0.*1e-3;
ub=x0.*1e3;
%Define differential equation solver parameters
Rel_Tol=1e-4;
[x]=Optim_VEVP(x0,x_scales,lb,ub,ideal_strain_profiles{1,1},ideal_strain_profiles{2,1},ideal_strain_profiles{3,1},...
    stress_check{1,1},stress_check{2,1},stress_check{3,1},Rel_Tol);
x_op=x.*x_scales;
x_op(8)=0;
x_op(9)=0;
x_op(11)=0;
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