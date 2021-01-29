function [I, varargout] = stepCurrentChuck(t,varargin)

% I = stepCurrentChuck(t,T,A)
%   STEPCURRENTCHUCK returns the value of I at time t, given a step
%   function with the following parameter option:
%
%   I = stepCurrentChuck(t)
%   I returns as: 0 if t<0, 0.5 if t=0, and 1 if t>0. If t is a vector,
%   returns a vector of the same length, with I at each value.
%
%   I = stepCurrentChuck(t,T)
%   where T is a vector. I returns as:
%       0 if t<T(1), 0.5 if t=T(1), 1 if T(1)<t<T(2), 0.5 if t=T(2),
%       0 if T(2)<t<T(3), 0.5 if t=T(3), 1 if T(3)<t<T(4), 0.5 if t=T(4),
%       0 if T(4)<t<T(5), 0.5 if t=T(5), 1 if T(5)<t<T(6), and so on...
%
%   I = stepCurrentChuck(t,T,A)
%   where T is a vector and A is a scalar. I returns as:
%       0 if t<T(1), A/2 if t=T(1), A if T(1)<t<T(2), A/2 if t=T(2),
%       0 if T(2)<t<T(3), A/2 if t=T(3), A if T(3)<t<T(4), A/2 if t=T(4),
%       0 if T(4)<t<T(5), A/2 if t=T(5), A if T(5)<t<T(6), and so on...

if ~exist('t')
    error('stepCurrent requires at least one parameter, t');
end

if nargin > 1
    Tv = varargin{1};
else
    Tv = 0;
end
if nargin > 2
    Iv = varargin{2};
else
    Iv = 1;
end

Lt = length(Tv);
switch length(Iv)
    case 1
        Iv = Iv * ones(1,Lt+1);
        Iv(find(mod([1:Lt+1],2))) = 0;
    case Lt
        Iv = [0,Iv];
    otherwise
        error('Parameter 3 must be a scalar or equal the length of parameter 2');
end

for k = 1:length(t)
    inds = find(Tv<t(k));
    if ~isempty(inds)
        pre = inds(end) + 1;
    else
        pre = 1;
    end
    if ~isempty(find(Tv==t(k)))
        I(k) = (Iv(pre) + Iv(pre + 1))/2;
    else
        I(k) = Iv(pre);
    end
end