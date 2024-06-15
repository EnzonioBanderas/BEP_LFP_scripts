function [HPfiltered_data] = HighPassfilter(order, cut_off, fs, dataset)
%
% [HPfiltered_data] = highpassfilter(infoHP, dataset)
%
%       Function for performing High-pass filter
%
%       Inputs:
%           
%           order:      order of the filter
%           cut_off:    cut-off frequency
%           fs:         sampling frequency
%           dataset:    array containing a signal in mV in each line
%
%       Outputs:
%           HPfiltered_data: array containing a filtered signal in mV in each line
%

cutoffrad = cut_off / (fs/2);
[b,a] = butter(order,cutoffrad,'high');

HPfiltered_data = filtfilt(b, a, dataset);

end

    

