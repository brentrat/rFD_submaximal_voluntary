function [cu_corrected, order_corrected]=Wcu(order, fc, type, fs)
%*Calculate correct cut-off frequency for n-pass Butterworth filter*
%Wcu(filter_order_following_n-passes, cut-off_frequency_in_Hz, filter_type[low, high, stop, bandpass], sampling_frequency)

% -------- Inputs: --------

% order         -> filter order (higher order is a more selective filter).
% The order that is provided as input to this function will be honoured when
% possible. Because the filter order doubles with each filter pass, a
% 2nd order filter will be generated by dual passing a 1st order filter.
% Uneven filter order are therefore not possible, if a uneven filter order
% is provided as input to this function a filter will be designed to the
% next even number. For example, if order = 3, then a 2nd order filter will
% be dual passed, resulting in a 4th order.

% type          -> is the desired cut-off frequency of the filter.

% filtertype -> can be 'low', 'high', 'stop', 'bandpass', corrosponding to
% a low-pass, high-pass, stop-band, pass-band filter type. However, Winter
% correction is not adapted yet for bandpass or stopband filter types.

% fs            -> sampling rate in samples/s

% -------- Output: --------

% cu_corrected    -> corrected cut-off (cu or fc) frequency.
% order_corrected -> corrected filter order for dual passing.

% -------- Description: --------
% This function corrects the cut-off frequency and order for passing the
% data through a Butterworth filter twice. Dual pass filtering is performed for the
% following reason: The Buterworth  
% filter creates time delays or phase lags that are not consistent across
% the frequency range after the data is passed through the filter.
% Therefore, a simple shift of the filtered data in time would not suffice
% as a correction. To compensate for these delays, it is common to pass the
% data through the filter again but in reverse, also referred to as
% dual-pass zero-lag or zero-phase shift filtering. However, because the
% data is filtered twice (once in each direction), the selectivity (order)
% of the filter changes and additional power is removed from the signal.
% With each filter pass, the cut-off frequency (the frequency at which the
% amplitude of the data is reduced to ~70.7% or -3dB) is pushed towards
% lower or higher frequencies, depending if the filter is designed to
% filter out high (low-pass filter) or low (high-pass filter) frequencies,
% respectively.

% Author:
% W. van den Hoorn (Feb 2020, for any errors please mail: w.vandenhoorn@uq.edu.au)
% Function tested using R2014b matlab version.

% References:
% Biomechanics and motor control of human movement, David A. Winter, 4th Ed

% check if uneven filter order is requested
if order/2~=round(order/2)
    disp(['NB: dual pass filtering cannot have uneven filter order. The inputted filter order: ',num2str(order),' is adapted to the following filter order: ',num2str(order+1),'.'])
    order=order+1;
else
    order=order;
end

% correct filter order for dual passing:
order_corrected=order/2;

% The Winter correction factor, NB Winter correction (equation on page 69
% of Winter book) is only valid for 2nd order low pass filter type
npasses=2;
% Order correction derived from MOtoNMS function files (Griffith
% University) or see also:
% https://www.codeproject.com/Articles/1267916/Multi-pass-Filter-Cutoff-Correction

C=(((2^(1/npasses))-1)^(.5/order_corrected));

% This correction factor needs to be applied to the angular cut-off frequency
% (equation 3.8 Winter book)

% check for filter type:
switch type
    case 'low'
        if length(fc)>1
            disp(['please enter one cut-off frequency for ',type, ' type Butterworth filter.'])
            return
        end
        % equation 3.8 Winter book
        % wc reflects the angular cut-off frequency normalised by the
        % sampling rate. (factor 2 is cancelled out), and is prewarped to
        % compensate for non-linear shift of cut-off frequency when using
        % the bi-linear transformation
        cu_angular=tan(pi*fc/fs);   
        cu_angular_corrected=cu_angular/C;
        % transform the cut-off frequency back to frequency
        cu_corrected=atan(cu_angular_corrected)*fs/pi;
    case 'bandpass'
        disp('Winter correction is not suitable for bandpass and bandstop filters...  no correction applied to cut-off and filter order. Note that amplitude will be attenuated by 50% instead of 70.71% at the cut-of frequencies')
        cu_corrected=fc;
        order_corrected=order;
    case 'high'
        if length(fc)>1
            disp(['please enter one cut-off frequency for ',type, ' type Butterworth filter.'])
            return
        end
        cu_angular=tan(pi*fc/fs);
        % Because it is a high-pass filter the inverse of the Winter
        % correction needs to be applied
        cu_angular_corrected=cu_angular*C;
        % transform the cut-off frequency back to frequency
        cu_corrected=atan(cu_angular_corrected)*fs/pi;
    case 'stop'
        disp('Winter correction is not suitable for bandpass and bandstop filters...  no correction applied to cut-off and filter order. Note that amplitude will be attenuated by 50% instead of 70.71% at the cut-of frequencies')
        cu_corrected=fc;
        order_corrected=order;
    otherwise
        disp('Please enter one of the following filter tyes: ''high'', or ''low''.');
        return
end