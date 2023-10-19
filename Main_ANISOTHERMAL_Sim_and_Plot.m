clear all
close all
clc
plotting_cols_master=[1,0,0;
                      0,0,1;
                      0,1,0;
                      0,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify isothermal parameters and loading waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VEVP_optim_400.mat','x_op')
x_op_400=x_op;
clear x_op
load('VEVP_optim_500.mat','x_op')
x_op_500=x_op;
clear x_op
load('VEVP_optim_600.mat','x_op')
x_op_600=x_op;
clear x_op
%
param_set={x_op_400;x_op_500;x_op_600};
clear x_op_400 x_op_500 x_op_600
%
param_temperatures=[400;500;600];
%
param_funcs=cell(18,1);
for ii=1:1:18
    param_temps=zeros(size(param_temperatures,1),1);
    for jj=1:1:size(param_temperatures,1)
        data_temp=param_set{jj,1};
        param_temps(jj,1)=data_temp(ii,1);
        clear data_temp
    end
    fit_temp=polyfit(param_temperatures([1,3]),param_temps([1,3]),1);
    param_funcs{ii,1}=fit_temp;
%     hold on
%     plot(param_temperatures,param_temps,'bo')
%     plot(param_temperatures,polyval(fit_temp,param_temperatures),'k-')
%     hold off
%     pause 
%     close all
    clear param_temps fit_temp
end
%
load('OP12.mat');
%
waveform_raw=readmatrix('OP1.csv');
waveform_raw(:,3)=(waveform_raw(:,3).*500)+500;
waveform_n=50;
for ii=2:1:size(waveform_raw,1)
    time=linspace(waveform_raw(ii-1,1),waveform_raw(ii,1),waveform_n);
    strain=linspace(waveform_raw(ii-1,2),waveform_raw(ii,2),waveform_n);
    temperature=linspace(waveform_raw(ii-1,3),waveform_raw(ii,3),waveform_n);
    if ii==2
        waveform=[time',strain',temperature'];
    else
        time=time(:,2:end);
        strain=strain(:,2:end);
        temperature=temperature(:,2:end);
        waveform=[waveform;[time',strain',temperature']];
    end
    clear time strain temperature
end
%
NAMES=cell(2,1);
NAMES{1,1}='Exp.';
NAMES{2,1}='Model';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



thermal_expansion=1.41e-5;
t_strain=[0;diff(waveform(:,3)).*1.41e-5].*100;
t_strain=cumsum(t_strain);

waveform_t_strain=waveform;
waveform_t_strain(:,2)=waveform_t_strain(:,2)-t_strain;

Rel_Tol=1e-4;
[state_var]=Chab_ODE_Solve_TMF(param_funcs,{waveform_t_strain},Rel_Tol);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Select_cycle=1;
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
data_temp=data_proc{Select_cycle,1};
plot(data_temp(:,2),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2);
results_temp=state_var{Select_cycle,1};
plot(waveform(:,2),results_temp(:,4), 'k--','MarkerSize', 2, 'LineWidth', 2);
clear data_temp results_temp
%
xlabel('\epsilon (%)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
xlim([-0.8 0.8])
ylim([-600 600])
hold off
clear figure  plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Select_cycle=1;
%
figure = figure('Color',[1 1 1]);
%title(['Cycle ', num2str(Select_cycle)],'fontSize',14,'fontWeight','bold');
hold all
data_temp=data_proc{Select_cycle,1};
plot(data_temp(:,1),data_temp(:,3), 'bx','MarkerSize', 2, 'LineWidth', 2);
results_temp=state_var{Select_cycle,1};
plot(waveform(:,1),results_temp(:,4), 'k--','MarkerSize', 2, 'LineWidth', 2);
clear data_temp results_temp
%
xlabel('t (s)','fontSize',14,'fontWeight','bold');
ylabel('\sigma (MPa)','fontSize',14,'fontWeight','bold');
legend(NAMES,'fontSize',14,'fontWeight','bold','Location','best')
set(gca,'fontsize',14,'fontWeight','bold')
grid on
%xlim([-0.6 0.6])
ylim([-600 600])
hold off
clear figure  plotting_cols Select_cycle xlim ylim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%