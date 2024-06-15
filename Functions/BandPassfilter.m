function [BANDfiltered_data] = BandPassfilter(order, cut_off, fs, dataset)
%
% [BANDfiltered_data] = bandpassfilter(infoBAND, dataset)
%
%       Function for performing Band-pass filter
%
%       Inputs:
%           
%           order:      order of the filter
%           cut_off:    cut-off frequencies
%           fs:         sampling frequency
%           dataset:    array containing a signal in mV in each line
%
%       Outputs:
%           BANDfiltered_data: array containing a filtered signal in mV in each line
%
if any(size(cut_off)~=[1,2])
    error('cut_off input should be 1x2 double')
end

cutoffrad = cut_off / (fs/2);
[b,a] = butter(order,cutoffrad);

BANDfiltered_data = filtfilt(b, a, dataset);

end

