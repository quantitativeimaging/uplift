%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File name:    Datalogger_pp.m
%
%   Purpose  :    Read measurements of force and displacement as 
%                 function of time respectively, then interpolates
%                 to obtain synchronised data for plotting.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all

filenamelog = 'nov 7 type b 2d uplift h130d45 exp 6';

% Read from force and displacement output files
[time_force,value_force] = ReadDataLog([filenamelog, ' force.txt']);
[time_dis,value_dis] = ReadDataLog([filenamelog, ' disp.txt']);


% Start plotting only after both force and displacement are being measured
time_offset = max(min(time_force),min(time_dis));
time_force = time_force - time_offset;
time_dis = time_dis - time_offset;

% Stop plotting immediately after one of the measurement stops
tmax = min(max(time_force),max(time_dis));
number_point = max(length(time_force),length(time_dis))*10; % reduce for faster speed, increase for better interpolation
time = linspace(0,tmax,number_point);

% Interpolate raw data
force = InterpolateTime(time,time_force,value_force);
dis = InterpolateTime(time,time_dis,value_dis);

% Potentially some unit conversion
% Convert Displacement to mm (currently 10^-6 m)
value_dis = value_dis / 100;
dis = dis / 100;

% Find peak force displacement
force_peak = max(force(:,:));
id_Peak_F = find(force == force_peak);
disp_peak = min(dis(id_Peak_F));
disp_peak
force_peak

% Find Uplift Speed
% take displacement readings once displacement begins
zero_disp = find(value_dis <= 0.2);
id_reg_start = max(zero_disp);
disp_start = value_dis(id_reg_start);
t_start = time_dis(id_reg_start);
t_peak =time(id_Peak_F)

%disp_start
%t_start

id_reg_end = min(find(value_dis >= max(value_dis) - 0.2));
% Find location of max displacement (and steady displacement)
%id_reg_start
%id_reg_end
   
dis_move = value_dis(id_reg_start:id_reg_end);
time_move = time_dis(id_reg_start:id_reg_end);
% average_speed
% b = fitlm(time_move, dis_move);
p = polyfit(time_move, dis_move,1);
p
reg_y = polyval(p, time_move);
%forceMax = ['Peak force = ', num2str(force_peak)]
%text(disp_peak, force_peak, forceMax,'HorizontalAlignment','left')

% Plot
figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
subplot(2,2,1) % plot 2 by 1 figures
plot(time,force) % plot interpolated force v time
hold on
scatter(time_force,value_force) % plot exp force v time
legend('fitted','raw')
xlabel('time(s)')
ylabel('force / N')
title('Fitting of force measurement')

subplot(2,2,2)
plot(time,dis)
hold on
scatter(time_dis,value_dis)
plot(time_move, reg_y)
legend('fitted','raw')
xlabel('time(s)')
ylabel('displacement / mm')
title('Fitting of displacement measurement')

subplot(2,2,3)
plot(dis,force)
xlabel('displacement / mm')
ylabel('force / N')
title('Processed result')

subplot(2,2,4)
scatter(time_move,dis_move)
xlabel('time (s)')
ylabel('displacement / mm')
title('Uplift Speed Test')

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,['\bf ' ,filenamelog],'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'fontsize',14)

%% Calculating moving standard deviation
index_foo = zeros(size(force));
sigma =zeros(size(force)); 
movAverage = zeros(size(force));
sigma_fil =zeros(size(force)); 
movAverage_fil = zeros(size(force));
index_foo(1) = 0;
sigma(1) = 0; sigma_fil(1) = 0;
movAverage(1) = 0; movAverage_fil(1) = 0;;
interval = 100;

for i = 1:(size(force)-interval)
    index_foo(i) = i;
    sigma(i) = std(force(i:i+interval,1));
    movAverage(i) = mean(force(i:i+interval,1));
    
    if sigma(i) < 0.2
        sigma_fil(i) = sigma(i); 
        movAverage_fil(i) = movAverage(i); 
    end    
    
    
end 

figure(99)
% % yyaxis left
% % plot(sigma(:,1))
% % ylim([0 1]);
% % yyaxis right
% % plot(movAverage(:,1))
% % xlabel('time')
% % legend('moving sigma','moving average')

[ax h1 h2] = plotyy(time,sigma,time,movAverage,'plot'); % Used in R2015a. 
ylim(ax(1),[0 1]);
ylim(ax(2),[30 40]);
xlabel('time')
legend('moving sigma', 'moving average')

% % Select moving Averages at only the bottom 3 percentile 
figure(98)
% % yyaxis left
% % scatter(index_foo, sigma_fil)
% % ylim([0 1]);
% % yyaxis right
% % scatter(index_foo, movAverage_fil)
% % %ylim([30 40]);
% % xlabel('time')
% % legend('moving sigma','moving average')

[ax h3 h4] = plotyy(index_foo,sigma_fil,index_foo,movAverage_fil,'plot'); % Used in R2015a. 
ylim(ax(1),[0 1]);
ylim(ax(2),[30 40]);
xlabel('time')
legend('moving sigma', 'moving average')
