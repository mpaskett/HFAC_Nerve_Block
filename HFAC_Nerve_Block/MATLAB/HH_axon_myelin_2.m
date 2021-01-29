function dYdt = HH_axon_myelin_2(t,Y,varargin)
% inputs:
%   ts - time span
%   Y - V, m, h, n values
%   function handle
%   function arguments
%   freq - frequency (0 for DC)
%   amp - amplitude
%   a - radius (um)
%   myelinLength - myelin Length


% axon parameters
L =    10000;    % um (axon length)
% a =        5;    % um (axon radius)
% stimType = 0;    % 0 for DC, 1 for AC stim, 2 for no stim

% set membrane parameters
cm	=    10;	% fF/um2
gna	=  1900;	% pS/um2
gk	=   360;	% pS/um2
gl 	=     3;	% pS/um2
% gna	=  30e3;	% pS/um2
% gk	=   .08e4;	% pS/um2
% gl 	=   .007e4;	% pS/um2
Vna	=    50;	% mV
Vk	=   -77;	% mV
Vl	= -54.4;	% mV
ri  = 10^-6;    % TOhm um


if nargin < 7
    error('HH requires at least 7 inputs')
end
 % add frequency, amplitude, radius, internal resistance
 
if nargin > 2
    if nargin > 3
        args = varargin{2};
        Icap = feval(varargin{1},t,args{:});
    else
        Icap = feval(varargin{1},t);
    end
else
    Icap = 0;
end
%   freq - frequency (0 for DC)
%   amp - amplitude
%   a - radius (um)
%   r_int - internal resistance

freq = varargin{3};
amp  = varargin{4};
a    = varargin{5};
myelinLength = varargin{6};
nodeORLength = 1; %length of node of ranvier (um)

if myelinLength > 0  % Checks to see if myelinated or unmyelinated
     K = length(Y)/4; % number of compartments
%      L = K*myelinLength/250; % (um) for myelination 
     L = K*nodeORLength;
else
    K = length(Y)/4; % number of compartments
end

if nargin > 8
    error('Too many input arguments')
end

dYdt = zeros(size(Y)); % get size
% K = length(Y)/4; % number of compartments

%Current at each compartment
I = zeros(K,1);
I(1) = Icap(1);

x = L/K; % delta x, um
A = 2*pi*a*x; % Area, um^2


% Stim block calc
if t>30 && t<50
%     Istimulation = 500e-6;      % Stim current in A
    if myelinLength > 0
        axonL = ((myelinLength*K)+(K*nodeORLength))/1e6;
        Ktemp = K*(myelinLength+nodeORLength);
    else
        axonL = L/1e6;
        Ktemp = K;
    end
    Istimulation = amp*cos((2*pi*freq*t));      % Stim current in A
    electrodeDistance = 400e-6 ; % distance of stim electrode from membrane
    stimCompartment   = Ktemp/2;     % compartment being stimulated
    conductance = 1/(3.33);      % conductance of extracellular space (S/m)(from 333Ohm*cm)
    [Vout, ~] = fieldPotentialCalc(Istimulation, axonL, Ktemp, ...
    electrodeDistance, stimCompartment, conductance);
    if Ktemp ~= K
        Vout = Vout(myelinLength:myelinLength:end);
    end
else
    Vout = zeros(K,1);
end
Vout = -Vout;   % magically fixes everything

k0 = [0:K-1]*4; %initialize compartments
Vin = Y(k0+1) + Vout;
Vm  = Y(k0+1);
m   = Y(k0+2);
h   = Y(k0+3);
n   = Y(k0+4);


if mod(floor(t),1) == 0
    disp(t);
end

% voltage dependent alpha and beta values
am = alfa_m(Vm);
ah = alfa_h(Vm);
an = alfa_n(Vm);
bm = beta_m(Vm);
bh = beta_h(Vm);
bn = beta_n(Vm);

% voltage dependent gate time constant values
taum = 1./ (am+bm);
tauh = 1./ (ah+bh);
taun = 1./ (an+bn);

% steady-state gating variables
minf = am .* taum;
hinf = ah .* tauh;
ninf = an .* taun;

% dYdt for gating variables m,h, and n
dYdt(k0+2) = (minf-m)./taum;
dYdt(k0+3) = (hinf-h)./tauh;
dYdt(k0+4) = (ninf-n)./taun;

% axial resitance, membrane capacitance
Ri = ri * (x + myelinLength) / (pi* a^2); 
Cm = cm * A;

%conductance values
Gna = gna * A * m.^3 .*h; 
Gk = gk * A * n.^4;
Gl = gl * A;

%compartment indices
ia = 1;
ib = 2:K-1;
ic = K;
ka = k0(ia)+1;
kb = k0(ib)+1;
kc = k0(ic)+1;

% differential equation for the first, middle, and final compartments
dYdt(ka) = (Gna(ia).*(Vna - (Vin(ia)-Vout(ia))) + ...
            Gk(ia).* (Vk -  (Vin(ia)-Vout(ia))) + ...
            Gl     * (Vl -  (Vin(ia)-Vout(ia))) + ...
            (Vin(ia+1) - Vin(ia))/Ri + ...
            I(ia))/Cm;
dYdt(kb) = (Gna(ib).*(Vna - (Vin(ib)-Vout(ib))) + ...
            Gk(ib).* (Vk - (Vin(ib)-Vout(ib))) + ...
            Gl* (Vl - (Vin(ib)-Vout(ib))) + ...
            (Vin(ib+1) - Vin(ib))/Ri + ...
            (Vin(ib-1) - Vin(ib))/Ri +...
            I(ib))/Cm;
dYdt(kc) = (Gna(ic).*(Vna - (Vin(ic)-Vout(ic))) + ...
            Gk(ic).* (Vk - (Vin(ic)-Vout(ic))) + ...
            Gl* (Vl - (Vin(ic)-Vout(ic))) + ...
            (Vin(ic-1) - Vin(ic))/Ri + ...
            I(ic))/Cm;
        
% figure(5)
% plot(dYdt(1:4:end,:))
end
