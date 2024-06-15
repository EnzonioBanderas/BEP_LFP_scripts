function [band_data] = BandExtraction(minfreq, maxfreq, freqres, dataset)
%
% [band_data] = BandExtraction(minfreq, maxfreq, dataset)
%
%       Function extracting data in the 'minfreq' to 'maxfreq' band. 
%
%       Inputs:
%           
%           minfreq:	lower limit of the band (degrees)
%           maxfreq:    higher limit of the band (degrees)
%           freqres:    frequency resolution of FFT
%           dataset:    array containing FFT results 
%
%       Outputs:
%           band_data:  array containing the FFT values in the selected frequency band
%

min_index = (minfreq / freqres)+ 1; 
max_index = (maxfreq / freqres)+ 1;

band_data = dataset([min_index:max_index],:);

end


