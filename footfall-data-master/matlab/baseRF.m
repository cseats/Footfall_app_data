function [rF rfTimeHistory rmsFilteredMax] = baseRF(data,dt,weighting,doPlots,newFig)

% function [rF rfTimeHistory rmsFilteredMax] = baseRF(data,dt,weighting,doPlots)
%
%   Function to calculate maximum response factor and response factor time history
%   from wighted time history
%
%   INPUTS (non-optional)
%   data                - Acceleration time history, channels in columns
%   dt                  - Time step in seconds
%   weighting           - BS6841 weighting curve as STRING
%                           = 0 for DEFAULT 'none'
%
%   INPUTS (optional)
%   doPlot       = -1   - No plots
%                = 0    - Print to screen (DEFAULT)
%                = 1    - RF time history
%
%   OUTPUTS
%   rF                  - Single value max response factor
%   rfTimeHistory       - Filtered RF timehistory of signal
%
% by Nicholas Simpson
% 18/03/2009
%
% 23/07/2012 by Nicholas Simpson
% Updated to use rms2 for increased speed

disp('**********baseRF.m***********')
disp('OLD FUNCTION PLEASE USE dbRF')
disp('**********baseRF.m***********')

% Weight the data
if nargin < 3
    weighting = 'none';
end

if nargin < 4
    doPlots = 0;
end

if nargin < 5
    newFig = 1;
end

timeEnd = (size(data,1)-1)*dt;

if strcmpi(weighting,'none') == 0
    [vDv weightedData,dt] = bs6841vdv(data,1/dt,weighting,0);
    data = weightedData;
end

% BEWARE bs6841vdv can decimate, so construct a new time vector


time = 0:dt:(size(data,1)-1)*dt;

% Calcualate RMS
[rmsFiltered] = rms2(data,dt,1); % Want max RMS in the filtered signal
rmsFilteredMax = max(rmsFiltered);

rfScaleFactor = 0.005;

rfTimeHistory = rmsFiltered ./ rfScaleFactor;
rF = rmsFilteredMax ./ rfScaleFactor; % Now have RF for each channel

if doPlots >= 0
    % Print calcualted RF to screen
    fprintf('\n-------------------------- RF --------------------------\n')
    fprintf('\n Channel\t\t  RF\t\t   Peak RMS acceleration')
    plotOneCounter = 0;
    for i = 1:size(rF,2)
        fprintf('\n\t%1.0f\t\t\t%5.3f\t\t\t\t%6.5f',(i-1),rF(i),rmsFilteredMax(i))
    end
    fprintf('\n--------------------------------------------------------\n')
end


if doPlots == 1
    % Plot rolling RF for multichannel

    % Build up multichannel plot command
    plotCommand = '';
    for i = 1:size(data,2)
        plotCommand = [plotCommand 'time,rfTimeHistory(:,' num2str(i) '),'];
    end
    plotCommand = ['plot(' plotCommand(1:(end-1)) ');'];
    if newFig == 1;
        figure
    else
        subplot(1,1,1);
    end
    eval(plotCommand);

    % Label the plot

    xlabel('Time, s')
    ylabel('Response Factor')
    title(['Response factor time history, weighted by BS6841 curve ' weighting])
    axis tight
    % Build up legend by looping through data columns
    legendText = '';
    for i = 1:size(data,2)
        legendText = [legendText '''Ch. ' num2str(i) ''','];
    end
    legendCommand = ['legend(' legendText(1:(end-1)) ');'];
    eval(legendCommand);
    set(gcf,'Name','RF time history')
end
