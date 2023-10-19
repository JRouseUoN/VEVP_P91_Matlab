clear all
close all
clc
plotting_cols_master=[1,0,0;
                      0,0,1;
                      0,1,0;
                      0,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify isothermal data to be used in viscoplastic parameter estiamtion
%and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VE_optim_600.mat','x_op')
VE_x=x_op;
clear x_op
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
cut_off_time=45000;
cut_off_strain=0.49;
%
%Approx. Yield
sig_y=225*0.7; %600
%sig_y=275*0.7; %500
%sig_y=325*0.7; %400
%
%Assumed K_VP
K_VP=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find first stress relaxation period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stress_relax_data=cell(size(Select_data_sets,1),1);
Loading_data=cell(size(Select_data_sets,1),1);
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{1,1};
    clear data data_proc
    %Limit time (cut out data after first hold period)
    cut_off_index=find(data_temp(:,1)<=cut_off_time);
    data_temp=data_temp(1:max(cut_off_index),1:3);
    clear cut_off_index
    %Limit strain (cut loading data)
    cut_off_index=find(data_temp(:,2)<=cut_off_strain);
    if sum(diff(cut_off_index)==1)==size(cut_off_index,1)-1
        cut_off_index=max(cut_off_index);
    else
        cut_off_index=min(find(diff(cut_off_index)~=1));
    end
    %
    Loading_data{ii,1}=data_temp(1:cut_off_index(1),:);
    %
    data_temp=data_temp(cut_off_index(1)+1:end,:);
    clear cut_off_index
    data_temp(:,1)=data_temp(:,1)-data_temp(1,1);
    Stress_relax_data{ii,1}=data_temp;
    clear data_temp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot stress relaxation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    data_temp=Stress_relax_data{ii,1};
    plot(data_temp(:,1),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
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
xlim([0 9000])
ylim([100 200])
hold off
clear figure plotting_cols xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find e_VP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_VP_data_raw=cell(size(Stress_relax_data));
%1 Time
%2 VP Strain
%3 Overstress
for ii=1:1:size(E_VP_data_raw,1)
    data_temp=Stress_relax_data{ii,1};
    E_VP=zeros(size(data_temp,1),3);
    E_VP(:,1)=data_temp(:,1);
    %
    E_E=data_temp(:,3)./VE_x(12);
    E_VE_1=(data_temp(:,3)./VE_x(13)).*...
        (1-exp((-data_temp(:,1).*VE_x(13))./VE_x(14)));
    E_VE_2=(data_temp(:,3)./VE_x(15)).*...
        (1-exp((-data_temp(:,1).*VE_x(15))./VE_x(16)));
    E_VE_3=(data_temp(:,3)./VE_x(17)).*...
        (1-exp((-data_temp(:,1).*VE_x(17))./VE_x(18)));
    E_E_E_VE=E_E;%+E_VE_1+E_VE_2+E_VE_3;
    %E_E_E_VE=E_E_E_VE.*100;
    %
    E_VP(:,2)=0.005-E_E_E_VE;
    %
    E_VP(:,3)=data_temp(:,3)-sig_y;
    %
    E_VP_data_raw{ii,1}=E_VP;
    clear data_temp E_VP E_E E_VE_1 E_VE_2 E_VE_3 E_E_E_VE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot viscoplastic data, smooth, and estimate rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    data_temp=E_VP_data_raw{ii,1};
    %
    data_temp(find(data_temp(:,3)<=0),:)=[];
    %
    plot(data_temp(:,1),data_temp(:,2), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
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
ylabel('\epsilon_{VP} (abs)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
%xlim([0 9000])
%ylim([50 150])
hold off
clear figure plotting_cols xlim ylim
[VP_fit_x,VP_fit_y]=ginput;
VP_fit_y(find(VP_fit_x<=0))=[];
VP_fit_x(find(VP_fit_x<=0))=[];
VP_smooth=fit(VP_fit_x,VP_fit_y,'exp2');
%
E_VP_data_proc=cell(size(Stress_relax_data));
%1 Time
%2 VP Strain Rate
%3 Overstress
for ii=1:1:size(E_VP_data_proc,1)
    data_temp=E_VP_data_raw{ii,1};
    data_temp(find(data_temp(:,3)<=0),:)=[];
    %
    data_temp(:,2)=((VP_smooth.a*VP_smooth.b).*exp(data_temp(:,1).*VP_smooth.b))+...
        ((VP_smooth.c*VP_smooth.d).*exp(data_temp(:,1).*VP_smooth.d));
    %
    E_VP_data_proc{ii,1}=data_temp;
    clear data_temp 
end
%
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    data_temp=E_VP_data_raw{ii,1};
    %
    data_temp(find(data_temp(:,3)<=0),:)=[];
    %
    plot(data_temp(:,1),data_temp(:,2), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
%     if ii==1
%         max_strain=max(data_temp(:,2));
%         max_stress=max(data_temp(:,3));
%     else
%         max_strain=max([max_strain;max(data_temp(:,2))]);
%         max_stress=max([max_stress;max(data_temp(:,3))]);
%     end
    clear data_temp
end
data_temp=E_VP_data_raw{1,1};
data_temp(find(data_temp(:,3)<=0),:)=[];
plot(data_temp(:,1),feval(VP_smooth,data_temp(:,1)), 'k--','MarkerSize', 2, 'LineWidth', 2);
clear data_temp
%
xlabel('t (s)','fontSize',14,'fontWeight','bold');
ylabel('\epsilon_{VP} (abs)','fontSize',14,'fontWeight','bold');
legend([NAMES;'Fit'],'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
%xlim([0 9000])
%ylim([50 150])
hold off
clear figure plotting_cols xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot viscoplastic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    data_temp=E_VP_data_proc{ii,1};
    plot(log(sinh(data_temp(:,3)./K_VP)),log(data_temp(:,2)), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    clear data_temp
end
%
xlabel('ln(sinh(\sigma_{VP}/K_{VP}))','fontSize',14,'fontWeight','bold');
ylabel('ln(\epsilon_{VP}^.)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
%xlim([0 9000])
%ylim([50 150])
hold off
clear figure plotting_cols xlim ylim
%
[VP_proc_fit_x,VP_proc_fit_y]=ginput;
FIT_MAT=ones(size(VP_proc_fit_x,1),2);
%FIT_MAT(:,2)=VP_proc_fit_x-log(2);
FIT_MAT(:,2)=VP_proc_fit_x;
res=FIT_MAT\VP_proc_fit_y;
ep_VP_0=exp(res(1));
m_VP=res(2);
clear res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot flow rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    data_temp=E_VP_data_proc{ii,1};
    plot(log(data_temp(:,3)),log(data_temp(:,2)), 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    clear data_temp
end
data_temp=transpose(linspace(0,150,100));
data_temp(:,2)=ep_VP_0.*(sinh(data_temp(:,1)/K_VP).^m_VP);
plot(log(data_temp(:,1)),log(data_temp(:,2)), 'k--','MarkerSize', 2, 'LineWidth', 2);
clear data_temp
%
xlabel('ln(\sigma_{VP})','fontSize',14,'fontWeight','bold');
ylabel('ln(\epsilon_{VP}^.)','fontSize',14,'fontWeight','bold');
legend([NAMES;'sinh Fit'],'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
xlim([0 4])
%ylim([50 150])
hold off
clear figure plotting_cols xlim ylim
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
% % % % % % X_in=VE_x;
% % % % % % 
% % % % % % X_in(1)=7.54e3;
% % % % % % X_in(2)=68.48;
% % % % % % X_in(3)=26.2e3;
% % % % % % X_in(4)=1157.8;
% % % % % % 
% % % % % % X_in(5)=ep_VP_0;
% % % % % % X_in(6)=1/m_VP;
% % % % % % X_in(7)=K_VP;
% % % % % % 
% % % % % % X_in(10)=sig_y;
% % % % % % 
% % % % % % clearvars -except X_in

