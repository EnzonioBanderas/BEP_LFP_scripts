function [LPfiltered_data] = LowPassfilter(order, cut_off, fs, dataset)
%
% [LPfiltered_data] = lowpassfilter(infoLP, dataset)
%
%       Function for performing Low-pass filter
%
%       Inputs:
%           
%           order:      order of the filter
%           cut_off:    cut-off frequency
%           fs:         sampling frequency
%           dataset:    array containing a signal in mV in each line
%
%       Outputs:
%           LPfiltered_data: array containing a filtered signal in mV in each line
%

cutoffrad = cut_off / (fs/2);
[b,a] = butter(order,cutoffrad,'low');

LPfiltered_data = filtfilt(b, a, dataset);

end
