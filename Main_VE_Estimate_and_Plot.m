clear all
close all
clc



plotting_cols_master=[1,0,0;
                      0,0,1;
                      0,1,0;
                      0,1,1];

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %Specify isothermal data to be used in viscoelastic parameter estiamtion
% % % %and parameters
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Select_data_sets=cell(3,1);
% % % Select_data_sets{1,1}='int_dwell_0_1_600.mat';
% % % Select_data_sets{2,1}='int_dwell_0_01_600.mat';
% % % Select_data_sets{3,1}='int_dwell_0_001_600.mat';
% % % %
% % % Selected_strain_rate=zeros(3,1);
% % % Selected_strain_rate(1,1)=0.1;
% % % Selected_strain_rate(2,1)=0.01;
% % % Selected_strain_rate(3,1)=0.001;
% % % %
% % % NAMES=cell(3,1);
% % % NAMES{1,1}='Exp. 0.1%/s';
% % % NAMES{2,1}='Exp. 0.01%/s';
% % % NAMES{3,1}='Exp. 0.001%/s';
% % % %
% % % cut_off_time=9000;
% % % cut_off_strain=0.0;
% % % strain_tol=0.009;
% % % %
% % % %Test Young's moduli - no longer used
% % % %E=1.2254e11; %600C
% % % %E=1.4619e11; %500C
% % % %E=1.780e11; %400C
% % % %E=E*1.3;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %Find target stress relaxation period
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Stress_relax_data=cell(size(Select_data_sets,1),1);
% % % Loading_data=cell(size(Select_data_sets,1),1);
% % % Total_data=cell(size(Select_data_sets,1),1);
% % % for ii=1:1:size(Select_data_sets,1)
% % %     load(Select_data_sets{ii,1});
% % %     data_temp=data_proc{1,1};
% % %     clear data data_proc
% % %     %Limit time (cut out data after target hold period)
% % %     cut_off_index=find(data_temp(:,1)<=cut_off_time);
% % %     data_temp=data_temp(1:max(cut_off_index),1:3);
% % %     clear cut_off_index
% % % 
% % % 
% % %     %Limit strain (cut loading data)
% % %     cut_off_index=abs(data_temp(:,2)-cut_off_strain);
% % %     cut_off_index(find(cut_off_index<=strain_tol))=0;
% % %     cut_off_index=max(find(cut_off_index==0));
% % % %     cut_off_index=find(data_temp(:,2)<=cut_off_strain);
% % % %     if sum(diff(cut_off_index)==1)==size(cut_off_index,1)-1
% % % %         cut_off_index=max(cut_off_index);
% % % %     else
% % % %         cut_off_index=min(find(diff(cut_off_index)~=1));
% % % %     end
% % %     %
% % %     data_temp=data_temp(cut_off_index(1):end,:);
% % %     Total_data{ii,1}=data_temp;
% % %     %
% % %     cut_off_index=find(data_temp(:,3)==max(data_temp(:,3)));
% % %     %
% % %     Loading_data{ii,1}=data_temp(1:cut_off_index(1),:);
% % %     %
% % %     data_temp=data_temp(cut_off_index(1):end,:);
% % %     clear cut_off_index
% % %     %Reference time to 
% % %     data_temp(:,1)=data_temp(:,1)-data_temp(1,1);
% % %     Stress_relax_data{ii,1}=data_temp;
% % %     clear data_temp
% % % end
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  
% % % 
% % % hold on
% % % data_temp=Stress_relax_data{1,1};
% % % plot(data_temp(:,1),data_temp(:,3),'rx')
% % % clear data_temp
% % % data_temp=Stress_relax_data{2,1};
% % % plot(data_temp(:,1),data_temp(:,3),'bx')
% % % clear data_temp
% % % data_temp=Stress_relax_data{3,1};
% % % plot(data_temp(:,1),data_temp(:,3),'gx')
% % % clear data_temp
% % % hold off
% % % 
% % % 
% % % 
% % % clearvars -except Stress_relax_data Loading_data  Total_data
% % % 
% % % Stress_relax_data_400=Stress_relax_data;
% % % Loading_data_400=Loading_data;
% % % Total_data_400=Total_data;
% % % 
% % % clear Stress_relax_data 
% % % clear Loading_data  
% % % clear Total_data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load VE data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VE_EXP_Data.mat');
Temperatures=[600;500;400];
%Data stored in 0.1%/s, 0.01%/s, 0.001%/s format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate E0 values from fastest strain rate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0_estimates=zeros(3,1);
for ii=1:1:3
    data_name_temp=strcat('Total_data_',string(Temperatures(ii)));
    data_temp=eval(data_name_temp);
    %
    highest_rate_date_temp=data_temp{1,1};
    cut_off_strain_temp=highest_rate_date_temp(1,2)+0.06;
    cut_off_ind_temp=find(abs(highest_rate_date_temp(:,2)-cut_off_strain_temp)==min(abs(highest_rate_date_temp(:,2)-cut_off_strain_temp)));
    E0_estimates(ii,1)=(highest_rate_date_temp(cut_off_ind_temp(1),3)-highest_rate_date_temp(1,3))/(highest_rate_date_temp(cut_off_ind_temp(1),2)-highest_rate_date_temp(1,2));
    E0_estimates(ii,1)=E0_estimates(ii,1)*1.2;
    clear data_name_temp data_temp highest_rate_date_temp cut_off_strain_temp cut_off_ind_temp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot E0 values against loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:1:3
    data_name_temp=strcat('Total_data_',string(Temperatures(ii)));
    data_temp=eval(data_name_temp);
    %
    figure = figure('Color',[1 1 1]);
    title(strcat(string(Temperatures(ii)), '^oC'),'fontSize',14,'fontWeight','bold');
    hold on
    date_temp_temp=data_temp{1,1};
    plot(date_temp_temp(:,2)-date_temp_temp(1,2),date_temp_temp(:,3)-date_temp_temp(1,3), 'rx','MarkerSize', 2, 'LineWidth', 2);
    date_temp_temp=data_temp{2,1};
    plot(date_temp_temp(:,2)-date_temp_temp(1,2),date_temp_temp(:,3)-date_temp_temp(1,3), 'bx','MarkerSize', 2, 'LineWidth', 2);
    date_temp_temp=data_temp{3,1};
    plot(date_temp_temp(:,2)-date_temp_temp(1,2),date_temp_temp(:,3)-date_temp_temp(1,3), 'gx','MarkerSize', 2, 'LineWidth', 2);
    strain_temp=transpose(linspace(0,0.12,100));
    stress_temp=strain_temp.*E0_estimates(ii,1);
    plot(strain_temp,stress_temp, 'k--','MarkerSize', 2, 'LineWidth', 2);
    clear data_temp_temp
    NAMES_all={'0.1%/s';'0.01%/s';'0.001%/s';'Elastic'};
    xlabel('\epsilon (%)','fontSize',14,'fontWeight','bold');
    ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
    legend(NAMES_all,'fontSize',14,'fontWeight','bold','Location','best')
    set(gca,'fontsize',14,'fontWeight','bold')
    grid on
    xlim([0 0.15])
    ylim([0 ceil(max(stress_temp)/10)*10])
    hold off
    %
    clear data_name_temp data_temp date_temp_temp strain_temp stress_temp NAMES_all figure
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate total viscous strain components from loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate total viscous strain components from relaxation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot stress relaxation data (all temperatures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=[plotting_cols_master(1:3,:);plotting_cols_master(1:3,:);plotting_cols_master(1:3,:)];
plotting_cols=jet(9);
%
NAMES_all=cell(9,1);
kk=1;
for ii=1:1:3
    for jj=1:1:3
        NAMES_all{kk,1}=strcat(NAMES{jj,1},{' '},string(Temperatures(ii,1)),'^oC');
        kk=kk+1;
    end
end
clear ii jj kk
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
kk=1;
for ii=1:1:3
    for jj=1:1:3
        data_temp=eval(strcat('Stress_relax_data_',string(Temperatures(ii,1)),'{jj,1}'));
        plot(data_temp(:,1),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2,...
            'Color', plotting_cols(kk,:), 'MarkerEdgeColor', plotting_cols(kk,:), 'MarkerFaceColor', plotting_cols(kk,:));
%         if jj==1
%             plot(data_temp(:,1),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2,...
%                 'Color', plotting_cols(kk,:), 'MarkerEdgeColor', plotting_cols(kk,:), 'MarkerFaceColor', plotting_cols(kk,:));
%         elseif jj==2
%             plot(data_temp(:,1),data_temp(:,3), 'bo','MarkerSize', 2, 'LineWidth', 2,...
%                 'Color', plotting_cols(kk,:), 'MarkerEdgeColor', plotting_cols(kk,:), 'MarkerFaceColor', plotting_cols(kk,:));
%         elseif jj==3
%             plot(data_temp(:,1),data_temp(:,3), 'b-','MarkerSize', 2, 'LineWidth', 2,...
%                 'Color', plotting_cols(kk,:), 'MarkerEdgeColor', plotting_cols(kk,:), 'MarkerFaceColor', plotting_cols(kk,:));
%         end
        kk=kk+1;
        clear data_temp
    end
end
xlabel('t (s)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES_all,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
%xlim([0 cut_off_time])
%ylim([50 150])
hold off
clear figure plotting_cols xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Estimate E
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Use slowest data
% % loading_data_sizes=zeros(size(Loading_data,1),1);
% % for ii=1:1:size(Loading_data,1)
% %     loading_data_sizes(ii,1)=size(Loading_data{ii,1},1);
% % end
% % data_temp=Loading_data{find(loading_data_sizes==max(loading_data_sizes)),1};
% % clear loading_data_sizes ii
% % E=(data_temp(:,2)./100)\(data_temp(:,3).*1e6);
% % E=0.9*E;
% % clear data_temp
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smooth data and estimate rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting_cols=varycolor(size(Select_data_sets,1));
plotting_cols=plotting_cols_master(1:size(Select_data_sets,1),:);
%
VE_strain=cell(size(Select_data_sets,1),1);
VE_strain_rate=cell(size(Select_data_sets,1),1);
VE_stress=cell(size(Select_data_sets,1),1);
%
for ii=1:1:size(Select_data_sets,1)
    data_temp=Stress_relax_data{ii,1};
    %
    time=data_temp(:,1);
    stress=data_temp(:,3);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sample_vec=unique(floor(logspace(0,(log(size(time,1))/log(10)),500)));
%     fit_time=time(sample_vec);
%     fit_stress=stress(sample_vec);
%     %
%     figure = figure('Color',[1 1 1]);
%     hold all
%     plot(time,stress,'bx')
%     plot(time(sample_vec),stress(sample_vec),'ro')
%     hold off
%     clear figure
%     %
%     VE_function=fit(fit_time,fit_stress,'exp2');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure = figure('Color',[1 1 1]);
%     hold all
%     plot(time,stress)
%     [fit_time,fit_stress]=ginput;
%     hold off
%     clear figure
%     VE_function=fit(fit_time,fit_stress,'exp2');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VE_function=fit(time(1:floor(size(time,1)*(1/2))),stress(1:floor(size(time,1)*(1/2))),'exp2');
    VE_function=fit(time,stress,'exp2');
    %VE_function=fit(fit_time,fit_stress,'exp2');
    %
    stress_smooth=feval(VE_function,time);
    stress_smooth_rate=((VE_function.a*VE_function.b).*(exp(VE_function.b.*time)))+((VE_function.c*VE_function.d).*(exp(VE_function.d.*time)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure = figure('Color',[1 1 1]);
    hold all
    plot(time,stress, 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    plot(time,stress_smooth, 'k-','MarkerSize', 2, 'LineWidth', 2);
    NAMES_temp=cell(2,1);
    NAMES_temp{1,1}=NAMES{ii,1};
    NAMES_temp{2,1}='Fit';
    xlabel('t (s)','fontSize',14,'fontWeight','bold');
    ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
    legend(NAMES_temp,'fontSize',14,'fontWeight','bold','Location','best')
    set(gca,'fontsize',14,'fontWeight','bold')
    grid on
    hold off
    clear figure NAMES_temp 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    total_strain=(mean(data_temp(:,2)))/100;
    VE_strain_1=total_strain-((stress_smooth.*10^6)./E);
    VE_strain_rate_1=-((stress_smooth_rate.*10^6)./E);
    %
    VE_strain{ii,1}=VE_strain_1;
    VE_strain_rate{ii,1}=VE_strain_rate_1;
    VE_stress{ii,1}=stress_smooth;
    %
    clear data_temp time stress VE_function stress_smooth stress_smooth_rate total_strain VE_strain_1 VE_strain_rate_1
end
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
for ii=1:1:size(Select_data_sets,1)
    VE_strain_temp=VE_strain{ii,1};
    VE_strain_rate_temp=VE_strain_rate{ii,1};
    VE_stress_temp=VE_stress{ii,1};
    %
    x_plot=(VE_stress_temp./VE_strain_temp);
    y_plot=(VE_strain_rate_temp./VE_strain_temp);
    y_plot(find(x_plot<=0))=[];
    x_plot(find(x_plot<=0))=[];
    %
    plot(x_plot,y_plot, 'bx','MarkerSize', 2, 'LineWidth', 2,...
        'Color', plotting_cols(ii,:), 'MarkerEdgeColor', plotting_cols(ii,:), 'MarkerFaceColor', plotting_cols(ii,:));
    %    
    clear VE_strain_rate_temp VE_strain_temp VE_stress_temp x_plot y_plot
end
%
xlabel('\sigma/\epsilon_{VE} (MPa)','fontSize',14,'fontWeight','bold');
ylabel('\epsilon^._{VE}/\epsilon_{VE}','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
% xlim([0 6e5])
% ylim([0 5e-3])
grid on
hold off
clear figure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






