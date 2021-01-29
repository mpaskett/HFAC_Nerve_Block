%% Simulates HH axon with multiple compartments.

clearvars;
% close all
ode_opts = odeset('AbsTol',10^-6,'RelTol',10^-6,...
                  'MaxStep',.01);
              
Vin = -65;          %mV
a_m = alfa_m(Vin);
b_m = beta_m(Vin);
a_h = alfa_h(Vin);
b_h = beta_h(Vin);
a_n = alfa_n(Vin);
b_n = beta_n(Vin);

m = a_m./(a_m+b_m);
n = a_n./(a_n+b_n);
h = a_h./(a_h+b_h);

ts = [0 40];
T_start = [10 40];
% I_step = 250e6;
I_step = 250e5; % 250e5 for repeated spikes

K_array = [100];

for ii = 1:1
    K  = K_array(ii);
    y0 = [Vin m n h];
    y0 = repmat(y0,1,K);

    [t,Y] = ode15s(@(t,Y) HH_axon_MP(t,Y,@stepCurrentChuck,{T_start, I_step}), ts, y0, ode_opts);
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
    plot(t,Vcalc(:,i))
end
% plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1.5)
leg = legend('1', '25','50', '75', '100');
title(leg,'Compartment')
xlabel('Time (ms)')
ylabel('Membrane Voltage (mV)')
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