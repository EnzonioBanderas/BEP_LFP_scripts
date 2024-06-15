function data_linearfit = LineFit(data,n_w,n_o)
% Fits line with least squares method on data divided in overlapping windows 
% and averages overlapping fitted lines. 
% Input: Data vector, Number of points of window, Number of points of overlap
% Output: Linear fit of data
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

    [data_buffered,data_truncation]=buffer(data,n_w,n_o,'nodelay'); %overlap segments
    data_buffered_linearfit=data_buffered-detrend(data_buffered); %get linear fit
    data_linearfit=buffer_inv(data_buffered_linearfit,n_o,'mean'); %mean linear fit for overlapping segments
    if ~isempty(data_truncation)
        data_truncation_linearfit=data_truncation-detrend(data_truncation);
        data_linearfit=[data_linearfit;data_truncation_linearfit];
    end
    
end

