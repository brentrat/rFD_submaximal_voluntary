function [netWork,posWork,negWork] = work_calc(x, y)
%% Calculate work by determining area ~= 0 using the trapezoidal method
%   Positive + Negative work = Net work

if length(x) > length(y)
    x = x(1:length(y));
elseif length(x) < length(y)
    y = y(1:length(x));
end

netWork = trapz(x,y);
[pos, posNum] = bwlabel(y > 0);
[neg, negNum] = bwlabel(y < 0);

if netWork > 0
    posWork = trapz(x(pos==posNum), y(pos==posNum));
    negWork = posWork-netWork;
elseif netWork < 0
    negWork = trapz(x(pos==posNum), y(pos==posNum));
    posWork = negWork-netWork;
end

end