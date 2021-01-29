function [I, varargout] = sineCurrent(t, varargin)
% I = sineCurrent(T,F,Imax, T_on)
% SINECURRENT returns Imax * sine(2*pi*F*T). T may be a scalar or a time vector.
% F and Imax are optional. F defaults to 0.5. Imax defaults to 2.
if ~exist('t') error('sineCurrent requires at least one input parameter, t'); end
if nargin>1
    f = varargin{1};
else
    f = 0.1;
end
if nargin>2
    Imax = varargin{2};
else
    Imax = 1;
end
if nargin > 3
    T_on = varargin{3};
else
    T_on = 0;
end
if t > T_on
%     I = Imax * sin(2*pi*f*t);
    I = Imax * cos(2*pi*f*t);
else
    I = 0;
end
varargout{1} = f; % optional output of frequency
varargout{2} = Imax; % optional output of current amplitude