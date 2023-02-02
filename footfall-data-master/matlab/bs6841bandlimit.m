function [b,a] = bs6841bandlimit(w)
% function [b,a] = bs6841bandlimit(w)
% BS6841 frequency weighting filter design 
%
% Calculate analogue transfer function for band limiting filters
% b and a are coefficients of DESCENDING powers of s for numerator and 
% denominator respectively, starting with s^4
%  
% w is a single character string, 'b' through to 'g'
% specifying weighting to be used.
%
% The folliwing frequency weightings are available:
%
% Wb, Wc, Wd, We, Wf, Wg  
%
% PY 23/05/02. Ref BS 6841:1987 "Measurement and evaluation of human 
% exposure to whole-body mechanical vibration and repeated shock"
% 
% see also BS6841FREQWEIGHT

if nargin ~= 1
	error('Invalid number of input arguements')
end
if ~ischar(w) || length(w)~=1
	error('Invalid weighting')
end
a = zeros(1,5);
b = zeros(1,5);
if w == 'b'
	f = [0.4,100,16,16,2.5,4];
	Q = [0.71,0.55,0.9,0.95];
% 	K = 0.4;
elseif w == 'c'
	f = [0.4,100,8,8];
	Q = [0.71,0.63];
% 	K = 1.0;
elseif w == 'd'
	f = [0.4,100,2,2];
	Q = [0.71,0.63];
% 	K = 1.0;
elseif w == 'e'
	f = [0.4,100,1,1];
	Q = [0.71,0.63];
% 	K = 1.0;
elseif w == 'f'
	f = [0.08,0.63,1e6,0.25,0.0625,0.1];
	Q = [0.71,0.86,0.8,0.8];
% 	K = 0.4;
elseif w == 'g'
	f = [0.8,100,1.5,5.3];
	Q = [0.71,0.68];
% 	K = 0.42;
else
	error('Invalid weighting')
end

b(3) = 4*pi^2*f(2)^2;
c(3) = 1;
c(2) = 2*pi*f(1)/Q(1);
c(1) = 4*pi^2*f(1)^2;
d(3) = 1;
d(2) = 2*pi*f(2)/Q(1);
d(1) = 4*pi^2*f(2)^2;

a(1:3) = a(1:3) + c(1).*d;
a(2:4) = a(2:4) + c(2).*d;
a(3:5) = a(3:5) + c(3).*d;

% sort into descending rather than ascending order

a = atrflip(a);
b = atrflip(b);
