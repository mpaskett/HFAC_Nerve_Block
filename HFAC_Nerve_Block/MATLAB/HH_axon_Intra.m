function dYdt = HH_axon_Intra(t,Y,varargin)

% axon parameters
L =    10000;    % um (axon length)
a =       10;    % um (axon radius)
stimType = 0;    % 0 for DC, 1 for AC stim, 2 for no stim

% set membrane parameters
cm	=    10;	% fF/um2
gna	=  1200;	% pS/um2
gk	=   360;	% pS/um2
gl 	=     3;	% pS/um2
Vna	=    50;	% mV
Vk	=   -77;	% mV
Vl	= -54.4;	% mV
ri  = 10^-6;    % TOhm um


if nargin < 2
    error('HH requires atleast two inputs ts, a time span, and V_init, an initial membranevoltage')
end
if nargin > 2
    if nargin > 3
        args = varargin{2};
        Icap = feval(varargin{1},t,args{:});
    else
        Icap = feval(varargin{1},t);
    end
else Icap = 0;
end
if nargin > 4
    error('Too many input arguments')
end


dYdt = zeros(size(Y)); % get size
K = length(Y)/4; % number of compartments

%Current at each compartment
I = zeros(K,1);
I(1) = Icap(1);

% Insert Stim Block
midIndex = K/2;

if stimType == 1
    I(midIndex) =  300e8*sin((2*pi*1*t));   %AC Stim block (1kHz)
elseif stimType == 0
    I(midIndex) =  300e7;                   %DC Stim block
end

x = L/K; % delta x, um
A = 2*pi*a*x; % Area, um^2

% Stim block calc
if t>100
    Istimulation = 1e-6;      % Stim current in A
    %Istimulation = 500e-6*sin((2*pi*1*t));      % Stim current in A
    axonL = L/1e6;               % convert to (m)
    electrodeDistance = 5e-6 ; % distance of stim electrode from membrane
    stimCompartment   = K/2;     % compartment being stimulated
    conductance = 1/(3.33);      % conductance of extracellular space (S*m)(from 333Ohm*cm)
    [Vout, ~] = fieldPotentialCalc(Istimulation, axonL, K, ...
    electrodeDistance, stimCompartment, conductance);
    
    % just for error checking
%     Vcheck = Y(1:4:end);
%     figure(5); plot(Vcheck); title('Vm')
%     xlabel('Compartment')
%     ylabel('Voltage')
else
    Vout = zeros(K,1);
end

k0 = [0:K-1]*4; %initialize compartments
Vin = Y(k0+1) + Vout;
Vm  = Y(k0+1);
m   = Y(k0+2);
h   = Y(k0+3);
n   = Y(k0+4);


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
Ri = ri * x / (pi* a^2);
Cm = cm * A;

%conductance values
Gna = gna * A * m.^3 .*h; 
Gk = gk * A * n.^4;
Gl = gl * A;

%comparment indices
ia = 1;
ib = 2:K-1;
ic = K;
ka = k0(ia)+1;
kb = k0(ib)+1;
kc = k0(ic)+1;

% differential equation for the first, middle, and final compartments
% find dVdt
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
        
% more debugging
% if t > 20
%     dV = dYdt(1:4:end);
%     figure(4); plot(dV); title('dV')
%     xlabel('Compartment')
%     ylabel('Change in Voltage')
%     check = 1;
% end
end