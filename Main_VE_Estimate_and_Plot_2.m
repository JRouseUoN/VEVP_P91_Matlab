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
load('VE_perlim_600.mat')
%
Select_data_sets=cell(3,1);
Select_data_sets{1,1}='int_dwell_0_1_600.mat';
Select_data_sets{2,1}='int_dwell_0_01_600.mat';
Select_data_sets{3,1}='int_dwell_0_001_600.mat';
%
Selected_strain_rate=zeros(3,1);
Selected_strain_rate(1,1)=0.1;
Selected_strain_rate(2,1)=0.01;
Selected_strain_rate(3,1)=0.001;
%
NAMES=cell(3,1);
NAMES{1,1}='Exp. 0.1%/s';
NAMES{2,1}='Exp. 0.01%/s';
NAMES{3,1}='Exp. 0.001%/s';
%
cut_off_time=18000;
cut_off_strain=0.19;
%
E=1.2254e11; %600C
%E=1.4619e11; %500C
%E=1.780e11; %400C
%
E=E*1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Open figure and section in to 3 regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % openfig('strain_VE_400.fig')
% % [region_x,region_y]=ginput(6);
% % 
% % eta_VE=zeros(3,1);
% % E_VE=zeros(3,1);
% % 
% % for ii=1:1:3
% %     x_temp=region_x((2*(ii-1))+1:(2*ii));
% %     y_temp=region_y((2*(ii-1))+1:(2*ii));
% %     %
% %     hold on
% %     plot(x_temp,y_temp,'k-','MarkerSize', 2, 'LineWidth', 2)
% %     hold off
% %     %
% %     pp=fit(x_temp,y_temp,'poly1');
% %     eta_VE(ii)=1/(pp.p1);
% %     E_VE(ii)=((-1)*pp.p2)*eta_VE(ii);
% %     clear x_temp y_temp pp
% % end
% % tau_VE=eta_VE./E_VE;
% % save('VE_perlim_400.mat','eta_VE','E_VE','tau_VE')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x=zeros(18,1);
x(7)=1;
x(10)=1000000;
x(12)=E/1e6;
x(13)=E_VE(1);
x(14)=eta_VE(1);
x(15)=E_VE(2);
x(16)=eta_VE(2);
x(17)=E_VE(3);
x(18)=eta_VE(3);
%Define differential equation solver parameters
Rel_Tol=1e-4;
Results=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    %Solve Chaboche equations for each point in ideal_load_cycle_data
    Results{ii,1}=Chab_ODE_Solve(x,ideal_strain_profiles{ii,1},Rel_Tol);      
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
    plot(load_temp(:,2),results_temp(:,4), 'k-','MarkerSize', 2, 'LineWidth', 2);
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
%legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
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
    plot(load_temp(:,1),results_temp(:,4), 'k-','MarkerSize', 2, 'LineWidth', 2);
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
%legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
ylim([-500 500])
hold off
clear figure Select_data_sets NAMES plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%