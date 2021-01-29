%% Hunter's HH_original

clearvars;
% close all
ode_opts = odeset('AbsTol',10^-6,'RelTol',10^-6,...
                  'MaxStep',.01);
              
Vin = -65.0557;      %mV
m = 0.0531;
h = 0.3183;
n = 0.5890;

ts = [0 40];
DC_on = [5 40];
% I_step = 250e6;
DC_step = 300e6; % 300e6 for repeated spikes


K_array = [100];

for ii = 1:1
    K  = K_array(ii);
    y0 = [Vin m n h];
    y0 = repmat(y0,1,K);

    [t,Y] = ode15s(@(t,Y) HH_axon_original(t,Y,@stepCurrentChuck,{DC_on, DC_step}), ts, y0, ode_opts);
    Vcalc = Y(:,1:4:size(Y,2));
    mArray = Y(:,2:4:size(Y,2));
    nArray = Y(:,3:4:size(Y,2));
    hArray = Y(:,4:4:size(Y,2));
%     figure(1);
%     % subplot(2,2,ii);
%     set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
%     plot(t,Vcalc);
%     str = sprintf('%d Compartments', K);
%     title(str)
%     ylabel('Membrane Potential (mV)')
%     xlabel('Time (ms)')


end
% set(groot,'defaultAxesColorOrder','remove')

set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
figure(2); clf;hold on;
for i=1:K
    plot(t,Vcalc(:,i)-i)
end
% plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1.5)
% leg = legend('1', '25','50', '75', '100');
% title(leg,'Compartment')
fig = gca;
xlabel('Time (ms)')
ylabel('Relative Membrane Voltage (mV)')
fig.FontSize = 15;
title('Direct Current Conduction Block', 'FontSize', 30)
hold off;
set(groot,'defaultAxesColorOrder','remove')


% figure(3); clf;
% subplot(3,1,1)
% plot(t,mArray(:,50))
% title('m')
% subplot(3,1,2)
% plot(t,nArray(:,50))
% title('n')
% subplot(3,1,3)
% plot(t,hArray(:,50))
% title('h')


%% No Conduction Block

clearvars;
% close all
ode_opts = odeset('AbsTol',10^-6,'RelTol',10^-6,...
                  'MaxStep',.01);
              
Vin = -65.0557;      %mV
m = 0.0531;
h = 0.3183;
n = 0.5890;

ts = [0 40];
DC_on = [10 40];
% I_step = 250e6;
DC_step = 300e6; % 300e6 for repeated spikes
AC_freq = 1; % kHz
AC_amp = 0; % uA
AC_on = 20; % ms


K_array = [100];

for ii = 1:1
    K  = K_array(ii);
    y0 = [Vin m n h];
    y0 = repmat(y0,1,K);

    [t,Y] = ode15s(@(t,Y) HH_axon_IntraMP(t,Y,@stepCurrentChuck,{DC_on, DC_step}, ...
        @sineCurrent,{AC_freq,AC_amp,AC_on}), ...
        ts, y0, ode_opts);
    Vcalc = Y(:,1:4:size(Y,2));
    mArray = Y(:,2:4:size(Y,2));
    nArray = Y(:,3:4:size(Y,2));
    hArray = Y(:,4:4:size(Y,2));
%     figure(1);
%     % subplot(2,2,ii);
%     set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
%     plot(t,Vcalc);
%     str = sprintf('%d Compartments', K);
%     title(str)
%     ylabel('Membrane Potential (mV)')
%     xlabel('Time (ms)')


end
% set(groot,'defaultAxesColorOrder','remove')

set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
figure(2); clf;hold on;
for i=1:K
    plot(t,Vcalc(:,i)+i)
end
% plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1.5)
% leg = legend('1', '25','50', '75', '100');
% title(leg,'Compartment')
fig = gca;
xlabel('Time (ms)')
ylabel('Relative Membrane Voltage (mV)')
fig.FontSize = 15;
title('No Conduction Block', 'FontSize', 20)
hold off;
set(groot,'defaultAxesColorOrder','remove')


figure(3); clf;
subplot(3,1,1)
plot(t,mArray(:,50))
title('m')
subplot(3,1,2)
plot(t,nArray(:,50))
title('n')
subplot(3,1,3)
plot(t,hArray(:,50))
title('h')

%% DC Conduction Block

clearvars;
% close all
ode_opts = odeset('AbsTol',10^-6,'RelTol',10^-6,...
                  'MaxStep',.01);
              
Vin = -65.0557;      %mV
m = 0.0531;
h = 0.3183;
n = 0.5890;

ts = [0 60];
DC_on = [10 40];
% I_step = 250e6;
DC_step = 320e6; % 300e6 for repeated spikes
AC_freq = 0; % kHz
AC_amp = -400e6; % uA
AC_on = 16; % ms


K_array = [100];

for ii = 1:1
    K  = K_array(ii);
    y0 = [Vin m n h];
    y0 = repmat(y0,1,K);

    [t,Y] = ode15s(@(t,Y) HH_axon_IntraMP(t,Y,@stepCurrentChuck,{DC_on, DC_step}, ...
        @sineCurrent,{AC_freq,AC_amp,AC_on}), ...
        ts, y0, ode_opts);
    Vcalc = Y(:,1:4:size(Y,2));
    mArray = Y(:,2:4:size(Y,2));
    nArray = Y(:,3:4:size(Y,2));
    hArray = Y(:,4:4:size(Y,2));
%     figure(1);
%     % subplot(2,2,ii);
%     set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
%     plot(t,Vcalc);
%     str = sprintf('%d Compartments', K);
%     title(str)
%     ylabel('Membrane Potential (mV)')
%     xlabel('Time (ms)')


end
% set(groot,'defaultAxesColorOrder','remove')

% set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
figure(2); clf;hold on;
% for i=1:K
%     plot(t,Vcalc(:,i)+i)
% end
plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1.5)
leg = legend('1', '25','50', '75', '100');
title(leg,'Compartment')
fig = gca;
xlabel('Time (ms)')
ylabel('Relative Membrane Voltage (mV)')
fig.FontSize = 15;
title('No Conduction Block', 'FontSize', 20)
hold off;
% set(groot,'defaultAxesColorOrder','remove')


figure(3); clf;
subplot(3,1,1)
plot(t,mArray(:,50))
title('m')
subplot(3,1,2)
plot(t,nArray(:,50))
title('n')
subplot(3,1,3)
plot(t,hArray(:,50))
title('h')