clear all
close all
clc
plotting_cols_master=[1,0,0;
                      0,0,1;
                      0,1,0;
                      0,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Determine .mat files in folder
% listing = dir;
% listing=struct2cell(listing);
% listing=listing(:,3:end);
% %Determine if any other file formats (apart from .csv) are present and
% %remove from list
% jj=1;
% for ii=1:1:size(listing,2)
%     name_temp=listing{1,ii};
%     if sum(name_temp(1,size(name_temp,2)-3:end)=='.mat')~=4
%         name_check(jj,1)=ii;
%         jj=jj+1;
%     end
% end
% if exist('name_check')==1
%     listing(:,name_check)=[];
% end
% clearvars -except listing
% %
% for ii=1:1:size(listing,2)
%     name_temp=listing{1,ii};
%     load(name_temp);
%     data_proc=cell(max(data(:,4)),1);
%     data(:,2)=data(:,2)-data(1,2);
%     data(:,3)=data(:,3)-data(1,3);
%     for jj=1:1:size(data_proc,1)
%         data_proc{jj,1}=data(find(data(:,4)==jj),:);
%     end
%     save(name_temp,'data','data_proc')
%     clear name_temp data data_proc
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot stress/strain hyst. loop for n data sets (selected cycle)
Select_data_sets=cell(3,1);
Select_data_sets{1,1}='int_dwell_0_1_500.mat';
Select_data_sets{2,1}='int_dwell_0_01_500.mat';
Select_data_sets{3,1}='int_dwell_0_001_500.mat';
%
NAMES=cell(3,1);
NAMES{1,1}='Exp. 0.1%/s';
NAMES{2,1}='Exp. 0.01%/s';
NAMES{3,1}='Exp. 0.001%/s';
%
%plotting_cols=varycolor(size(Select_data_sets,1));
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
clear figure Select_data_sets NAMES plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot stress/time data for n data sets (selected cycle)
Select_data_sets=cell(3,1);
Select_data_sets{1,1}='int_dwell_0_1_500.mat';
Select_data_sets{2,1}='int_dwell_0_01_500.mat';
Select_data_sets{3,1}='int_dwell_0_001_500.mat';
%
NAMES=cell(3,1);
NAMES{1,1}='Exp. 0.1%/s';
NAMES{2,1}='Exp. 0.01%/s';
NAMES{3,1}='Exp. 0.001%/s';
%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
Select_cycle=2;
%
figure = figure('Color',[1 1 1]);
title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    load(Select_data_sets{ii,1});
    data_temp=data_proc{Select_cycle,1};
    data_temp(:,1)=data_temp(:,1)-data_temp(1,1);
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
ylim([-500 500])
hold off
clear figure Select_data_sets NAMES plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot stress/strain hyst. loop for n data sets (selected cycle)
Select_data_sets=cell(4,1);
Select_data_sets{1,1}='IP12.mat';
Select_data_sets{2,1}='IP34.mat';
Select_data_sets{3,1}='OP12.mat';
Select_data_sets{4,1}='OP12_low_strain.mat';
%
NAMES=cell(4,1);
NAMES{1,1}='IP 1';
NAMES{2,1}='IP 3';
NAMES{3,1}='OP 1';
NAMES{4,1}='OP 1 low strain';
%
%plotting_cols=varycolor(size(Select_data_sets,1));
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
    plot(data_temp(:,2),data_temp(:,3), 'bo','MarkerSize', 2, 'LineWidth', 2,...
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
xlabel('\epsilon (%)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
% xlim([-0.8 0.8])
% ylim([-600 400])
hold off
clear figure Select_data_sets NAMES plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
