function [af,afrms] = bs6841aweight(acc,fs,w,units,do_plots)
% function [af,afrms] = bs6841aweight(acc,fs,w,units,plotopts)
%
% This function produces and plots a weighted time histories and their rolling RMS values
% WARNING: Values suspect at frequencies above 25 Hz
% 
% INPUTS
% acc = acceleration time history (required)
% fs = sampling frequency (required)
% w = curve to use 'b' through to 'g' (required)
% units    =  1 = Acceleraation (default)
%          =  2 = RF (if input is SI)
% do_plots =  1 = rms + filtered signal (default)
%             2 = rms only
%             3 = filtered signal only
%             
% OUTPUTS
% af = a-weighted filtered time history
% afrms = 1s period rolling rms of a-weighted time history
%
% NJS 12/11/08
%
% Based HEAVILY on following...
%
% function [vdv,af] = bs6841vdv(aacc,fs,w,do_plots)
%
% BS6841 filtering and evaluation of vibration doasge values
% filter a using the band limiting and frequency weighting filters 
% specified in BS6841.
%
% fs is the sampling frequency
%
% w is a single character string, 'b' through to 'g'
% specifying weighting to be used.
%
% The folliwing frequency weightings are available:
%
% Wb, Wc, Wd, We, Wf, Wg  
%                                                    
% af is the filtered acceleration, which may have a different time base
% if decimation is used
%
% Sampling frequency is checked to ensure that it is adequate
% for the band chosen. An error is issued if the Nyquist frequency 
% is less than the ~3dB point of the combined filter, and a warning is 
% issued if it is less than approx 20dB down point.
%
% Signal length is checked to ensure that there is at least 2 whole 
% period  at the lower 20dB down point. For all except the 'b'
% weighting this also implies that there are 6 whole periods at the 3dB
% down point. For the 'b' weighting, because of the step in the frequency
% response, this implies 28 periods at the lower 3dB point.
%  
% Over-sampled signals are decimated to ensure filter stability, but not
% too drastically to try and ensure that pre-warping is not necessary
% for the bilinear transform from the s domain to the z domain. The
% responses of the analogue filters (as specified in the standard) and 
% their digital realisations (as implimented by Matlab) should be compared.
% If there appear to be significant discrepancies then either a) try
% a higher sampling frequency, b) consider modifying these routines to 
% include pre-warp in the bilinear transform
%
% PY, 28/05/02
%
% See also BS6841BANDLIMIT, BS6841FREQWEIGHT

% check input arguements
if nargin < 5
  do_plots = 1;
end
if nargin < 4
	error('Invalid number of input arguements')
end
if nargin < 3
	units = 1;
end
if ~ischar(w) || length(w)~=1  % CHANGED FROM isstr 06/11/2008
	error('Invalid weighting')
end

% Acceleration to RF conversion
rffactor = 0.005; % m/s2 per RF

% check sampling frequency 
% set fminsoft, fminhard and Tmin
% frequency response of combined band limit and frequency weighting 
% filters are approx 3dB down at fminhard and 20dB down at fminsoft
% If sampling frequency is less than twice fminhard, am error occurs,
% if it is less than fminsoft, a warning is issued that sampling 
% frequency may be insufficient to represent all frequencies in that
% particular band.
% Signal length is checked to ensure that there are at least 2 periods 
% of Tmin
np = length(acc);
if w == 'b'
	fminsoft = 100;
	fminhard = 10;
	Tmin = 1/0.2;
elseif w == 'c'
	fminsoft = 50;
	fminhard = 7;
	Tmin = 1/0.1;
elseif w == 'd'
	fminsoft = 11;
	fminhard = 2;
	Tmin = 1/0.1;
elseif w == 'e'
	fminsoft = 10;
	fminhard = 1;
	Tmin = 1/0.1;
elseif w == 'f'
	fminsoft = 0.8;
	fminhard = 0.2;
	Tmin = 1/0.04;
elseif w == 'g'
	fminsoft = 50;
	fminhard = 7;
	Tmin = 1/0.4;
else
	error('Invalid weighting')
end
if fs < 2*fminhard
	fprintf('\nIncrease sampling frequency to at least %3.2f Hz\n',fs,2*fminhard)
	error('Sampling frequency (%3.2f Hz) insufficient to represent\n vibration for this weighting',fs)
elseif fs < 2*fminsoft && do_plots~=-2
	fprintf('\nWarning, sampling frequency (%3.2f Hz) insufficient\nto represent higher frequencies for this weighting',fs)
	fprintf('\nAlso performance of digital filter at higher frequencies may be compromised')
	fprintf('\nCheck by comparing analogue and digital filter characteristics') 
	fprintf('\nType help bs6841vdv for more info') 
	fprintf('\nIncrease sampling frequency to at least %3.2f Hz wil help\n',2*fminsoft)
end
% check length of signal
dt = 1/fs;
T = dt * np;
if T < Tmin 
	fprintf('\nA signal length of at least %2.1 secs is required',Tmin)
	error('Insufficient signal length to give\n representative values of rms/rmq etc.')
end

% decimate if fs > 20*fminsoft
fsold = fs;
while fs > 20*fminsoft
    for i = 1:size(acc,2)
        accdec(:,i) = decimate(acc(:,i),2);
    end
    acc = accdec;
    clear accdec;
    fs = fs/2;
end
if fsold ~= fs
	fprintf('\nN.B. Decimation used to reduce\n sampling frequency from %f to %f \n',fsold,fs)
end

dt = 1/fs; % dt may have change in decimation

% get transfer function (analogue) frequency weighting

[bfw,afw] = bs6841freqweight(w);


% get transfer function (analogue) band limiting (descending powers of s)

[bbl,abl] = bs6841bandlimit(w);

% calculate analogue frequency responses (descending powers of s)
f = logspace(-2,3,1024);

% combine
a = zeros(1,9);
b = zeros(1,9);
for jj = 1:5
	a(jj:jj+4) = a(jj:jj+4) + afw(jj).*abl(1:5);
	b(jj:jj+4) = b(jj:jj+4) + bfw(jj).*bbl(1:5);
end

% convert to z transform
[bz,az] = bilinear(b,a,fs);

% Filter the signal
af = filter(bz,az,acc);

% % % Filter with a 25 Hz low pass to ignore unreliable data
% % af = butt(af,dt,12,25,1,0);

% Calculate rolling RMS of af over 1 second intervals
afrms = rms(af,dt,1);

% Construct the time vector
tt = dt.*(1:length(af));

if units ~= 1
    af = af/rffactor;
    afrms = afrms/rffactor;
    labely = 'RF';
else
    labely = 'Acceleration, m/s^2';
end

% Plot af time history
if do_plots == 1 % RMS + filtered
    plot(tt,af,tt,afrms)
    title('Filtered data and rolling RMS')
    xlabel('Time')
    ylabel(labely)
    legend('Filtered Signal','Rolling RMS (1s)')
elseif do_plots == 2 % RMS only
    plot(tt,afrms)
    title('a-weighted rolling RMS')
    xlabel('Time')
    ylabel(labely)
elseif do_plots == 3 % filtered only
    plot(tt,af)
    title('a-weighted filtered data')
    xlabel('Time')
    ylabel(labely)
else
end



