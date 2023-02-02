function [vdv,af,dt] = bs6841vdv(acc,fs,w,do_plots)
% function [vdv,af,dt] = bs6841vdv(acc,fs,w,do_plots)
%
% BS6841 filtering and evaluation of vibration dosage values (vdv)
% filter acc using the band limiting and frequency weighting filters
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
%                                                    T
%                                                (  (   4     )1/4
% vdv is the vibration dosage value evaluated as (  |  a  dt  )
%                                                (  )         )
%                                                    0
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
% period at the lower 20dB down point. For all except the 'b'
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
% do_plots = 1 to get filtered and non-filtered traces
% do_plots = 2 to get traces and filter characteristics
% do_plots = 0 for no_plots
% do_plots = -1 to plot only filter response for given weighting and
% sampling freq (a is ignored) into CURRENT figure window
% do_plots = -2 to do all analysis with minimal output to screen
%
% PY, 28/05/02
% Nicholas Simpson 10/11/14 - dt added to output
%
% See also BS6841BANDLIMIT, BS6841FREQWEIGHT

% check input arguments
if nargin<4
    do_plots=2;
end
if nargin <3
    error('Invalid number of input arguments')
end
if ~ischar(w) || length(w)~=1
    error('Invalid weighting')
end
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
    error('Sampling frequency (%3.2f Hz) insufficient to represent\n vibration for this weighting')
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
[fp,fm] = oafft(acc,1/fs,1,0);% get original spectral content befpore decimation
fsold = fs;
while fs > 20*fminsoft
    for i = 1:size(acc,2)
        accdec(:,i) = resample(acc(:,i),1,2);
    end
    acc = accdec;
    clear accdec;
    fs = fs/2;
end
if fsold ~= fs
    fprintf('\nN.B. Decimation used to reduce\n sampling frequency from %f to %f \n',fsold,fs)
end

% get transfer function (analogue) frequency weighting
[bfw,afw] = bs6841freqweight(w);

% get transfer function (analogue) band limiting (descending powers of s)
[bbl,abl] = bs6841bandlimit(w);

% calculate analogue frequency responses (descending powers of s)
f = logspace(-2,3,1024);
hfw = freqs(bfw,afw,f);
hbl = freqs(bbl,abl,f);
% combine
a = zeros(1,9);
b = zeros(1,9);
for jj = 1:5
    a(jj:jj+4) = a(jj:jj+4) + afw(jj).*abl(1:5);
    b(jj:jj+4) = b(jj:jj+4) + bfw(jj).*bbl(1:5);
end

% convert to z transform
[bz,az] = bilinear(b,a,fs);

% filter characteristics
if do_plots==1 || do_plots==-1
    hz = freqz(bz,az,f,fs);
    hs = freqs(b,a,f);
    if do_plots==-1
        subplot(2,1,1)
        loglog(f./(2*pi),abs(hfw),f./(2*pi),abs(hbl))
        grid on
        title(['Analogue Frequency Responses: ',w,' Weighting'])
        legend('Frequency Weighting','Band Limiting')
        axis([10^-2, 10^2,10^-2,10])
        ylabel('Filter Modulus')
        subplot(2,1,2)
        loglog(f./(2*pi),abs(hs),f,abs(hz))
        grid on
        legend('Analogue Filter','Digital Filter')
        title('Combined Frequency Responses')
        axis([10^-2, 10^2,10^-2,10])
        ylabel('Filter Modulus')
        af = [];
        vdv = [];
        disp('doplot = -1, processing of acceleration is ignored, returning empty ''vdv'' and ''af'' variables.');
        disp('Type help bs6841vdv for more info.');
        return
    end
end

% filter
af = filter(bz,az,acc);

% calculate vdv
dt = 1/fs;
vdv = (sum(dt*(af.^4))).^0.25;

% calculate frequency response of digital and analogue filters:

% plot weighting and band limiting responses, digital and analogue
% characteristics
if do_plots == 2
    hz = freqz(bz,az,f,fs);
    hs = freqs(b,a,f);
    [fpf,fmf] = oafft(af,1/fs,1,0);
    figure
    subplot(3,1,1)
    loglog(f./(2*pi),abs(hfw),f./(2*pi),abs(hbl))
    grid on
    title(['Analogue Frequency Responses: ',w,' Weighting'])
    legend('Frequency Weighting','Band Limiting')
    axis([10^-2, 10^2,10^-2,10])
    ylabel('Filter Modulus')
    subplot(3,1,2)
    loglog(f./(2*pi),abs(hs),f,abs(hz))
    grid on
    legend('Analogue Filter','Digital Filter')
    title('Combined Frequency Responses')
    axis([10^-2, 10^2,10^-2,10])
    ylabel('Filter Modulus')
    subplot(3,1,3)
    loglog(fp,fm,fpf,fmf)
    grid on
    legend('Original Signal','Filtered Signal')
    xlabel('Frequency, Hz')
    ylabel('FFT Modulus')
    %axis([10^-2, 10^2,10^-2,10])
end
if do_plots == 1 || do_plots == 2
    figure
    subplot(2,1,1)
    tt = dt.*(1:length(af));
    plot(tt,acc,tt,af)
    legend('Original Signal','Filtered Signal')
    ylabel('Acceleration, m/s^2')
    title('Filtered Signal for VDV calculation')
    subplot(2,1,2)
    plot(tt,af)
    % calculate crest factor
    peak = max(af);
    rms = (mean(af.^2)).^.5;
    cf = peak/rms;
    title(['Filtered Signal, Crest Factor = ',num2str(cf)])
    ylabel('Acceleration, m/s^2')
    xlabel('Time, secs')
end