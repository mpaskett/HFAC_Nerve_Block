%% Simulates HH axon with multiple compartments.

clearvars;
tic
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

ts = [0 35];
T_start = [0 50];
I_step = 5e4;     % "synaptic" input to initiate action potential (fA)

K_array = [100];     % set to 10 for myelinated, 100 for unmyelinated
%   ts - time span
%   Y - V, m, h, n values
%   function handle
%   function arguments
%   freq - frequency (0 for DC) in kHz
%   amp - amplitude in A
%   a - radius (um)
% Optional
%   Lm - myelin length

%%%%%%%%%%%%%_Conduction block paramaters_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 4;       % kHz
amp  = 2e-5;  % amplitude of input current for blocking stimulus(A)
a    = 8;      % axon radius (um)
Lm   = 1000;    % myelin length (um); 0 for unmyelinated, 1000 for myelinated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:1
    
K  = K_array(ii);
y0 = [Vin m n h];
y0 = repmat(y0,1,K);

[t,Y] = ode15s(@(t,Y) Copy_of_HH_axon_myelin_2(t,Y,@stepCurrentChuck,{T_start, I_step},freq,amp,a,Lm), ts, y0, ode_opts);
Vcalc = Y(:,1:4:size(Y,2));
mArray = Y(:,2:4:size(Y,2));
nArray = Y(:,3:4:size(Y,2));
hArray = Y(:,4:4:size(Y,2));
% figure(1);
% % subplot(2,2,ii);
% set(groot,'defaultAxesColorOrder',copper(size(Vcalc,2)))
% plot(t,Vcalc);
% str = sprintf('%d Compartments', K);
% title(str)
% ylabel('Membrane Potential (mV)')
% xlabel('Time (ms)')


end
% set(groot,'defaultAxesColorOrder','remove')

figure(2); clf;
% plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1.5)
plot(t,Vcalc(:,[1 (K_array/2 - ceil(.25*K_array)) ceil((K_array/2)) (K_array/2 + ceil(.25*K_array)) K_array]),'LineWidth', 1)
leg = legend('1', '25','50', '75', '100');
title(leg,'Compartment')
xlabel('Time (ms)')
ylabel('Membrane Voltage (mV)')

% plot parameters of middle compartment
plotCompartment = K_array/2;

figure(3); clf;
subplot(3,1,1)
plot(t,mArray(:,plotCompartment))
title('m')
subplot(3,1,2)
plot(t,nArray(:,plotCompartment))
title('n')
subplot(3,1,3)
plot(t,hArray(:,plotCompartment))
title('h')
toc
beep