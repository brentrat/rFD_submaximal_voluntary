function [fasStretch,fasStretchX] = fas_stretch(FL,varargin)
%Calculate fascicle stretch from a fascicle length-time trace.
%[LocalMaxStretch_CumulativeStretch_MinToMax] = fas_stretch(FL)
%*Requires rcumsum function available from: 
%https://au.mathworks.com/matlabcentral/fileexchange/28685-rcumsumc

% Author:
% BJ Raiteri, May 2022, if you find errors pls email brent.raiteri@rub.de
% Tested in R2022a.

if ~isrow(FL)
    FL = FL';
end

FLdiff = diff(FL);
FLstrX = find(FLdiff>0);
FLstr = FLdiff(FLstrX);
FLstrXdiff = diff(FLstrX)==1;
FLstrXdiff = [1 FLstrXdiff];

if length(FLstrXdiff) == 1
    FLstrLoop = 1;
else
    FLstrLocal = rcumsum(FLstrXdiff);
    FLstrLoop = find(FLstrLocal==0);
end

if isempty(FLstrLoop)

    FLstrLoop = 1;

end

for bb = 1:length(FLstrLoop)

    if bb == 1

        fasStr(bb,1) = sum(FLstr(1:FLstrLoop(bb)-1));

    elseif bb == length(FLstrLoop)

        fasStr(bb,1) = sum(FLstr(FLstrLoop(bb):end));

    else

        fasStr(bb,1) = sum(FLstr(FLstrLoop(bb):FLstrLoop(bb+1)-1));

    end

end

fasStretch = max(fasStr);
fasStretch(1,2) = sum(FLstr);

[FLmin,FLminX] = min(FL);
[FLmax,FLmaxX] = max(FL(FLminX:end));
FLmaxX = FLmaxX+FLminX-1;
fasStretchX = [FLminX FLmaxX];

fasStretch(1,3) = FLmax-FLmin;

if nargin > 1
    time = varargin{1};
    fasStretch(1,4) = fasStretch(1,3)/(time(FLmaxX)-time(FLminX)); 
end