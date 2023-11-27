function [y,yWinDur] = rmsDC(y, winDur, DC, fs)
%*Calculate a centered root-mean-square amplitude*
%rmsDC(signal, window_duration[ms], DCoffset=1, fps)

% Authors:
% Paolo Tecchio
% Updated by Brent J. Raiteri, May 2022.
% If you find errors please email: brent.raiteri@rub.de
% Tested in R2022a.

%% Calculate window length based on window duration and frames per second
winPts = winDur/(1/fs)/1000;

if rem(winPts,2) == 0
    winPts = winPts+1;
    yWinDur = winPts*(1/fs)*1000;
end

%% Remove DC offset
if DC == 1
    y = y-mean(y,'omitnan');
end

%% Calculate centred root-mean-square amplitude
if winDur == fs
    y = sqrt(mean(y.^2));
else
y = sqrt(movmean(y.^2, winPts));
end

% Subtract minimum from signal
%y = y-min(y);