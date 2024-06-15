function [HighLowPassfiltered_data] = HighLowPassfilter(order, cut_off, fs, dataset)
%
% [HighLowPassfilter] = bandpassfilter(infoBAND, dataset)
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
if any(size(cut_off)~=[1,2])||~isnumeric(cut_off)
    error('cut_off input should be 1x2 double')
end
if cut_off(2)-cut_off(1)<=0
    error('cut_off input should a frequency band')
end

dataset=HighPassfilter(order,cut_off(1),fs,dataset);
HighLowPassfiltered_data=LowPassfilter(order,cut_off(2),fs,dataset);

end

