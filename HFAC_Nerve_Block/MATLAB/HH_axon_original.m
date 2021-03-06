function dYdt = HH_axon_original(t,Y,varargin)

% axon parameters
L =    10000;    % um (axon length)
a =      100;    % um (axon radius)
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
if t > 10
    if stimType == 1
        I(midIndex) =  550e7*sin((2*pi*1*t));   %AC Stim block (1kHz)
    elseif stimType == 0
        I(midIndex) =  315e7;                   %DC Stim block
    end
end

x = L/K; % delta x, um
A = 2*pi*a*x; % Area, um^2

k0 = [0:K-1]*4; %initialize compartments
V = Y(k0+1);
m   = Y(k0+2);
h   = Y(k0+3);
n   = Y(k0+4);

% voltage dependent alpha and beta values
am = alfa_m(V);
ah = alfa_h(V);
an = alfa_n(V);
bm = beta_m(V);
bh = beta_h(V);
bn = beta_n(V);

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
dYdt(ka) = (Gna(ia).*(Vna - V(ia)) + ...
            Gk(ia).* (Vk - V(ia)) + ...
            Gl* (Vl - V(ia)) + ...
            (V(ia+1) - V(ia))/Ri + ...
            I(ia))/Cm;
dYdt(kb) = (Gna(ib).*(Vna - V(ib)) + ...
            Gk(ib).* (Vk - V(ib)) + ...
            Gl* (Vl - V(ib)) + ...
            (V(ib+1) - V(ib))/Ri + ...
            (V(ib-1) - V(ib))/Ri +...
            I(ib))/Cm;
dYdt(kc) = (Gna(ic).*(Vna - V(ic)) + ...
            Gk(ic).* (Vk - V(ic)) + ...
            Gl* (Vl - V(ic)) + ...
            (V(ic-1) - V(ic))/Ri + ...
            I(ic))/Cm;
end