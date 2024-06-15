function [PSD_matrix,t,f] = spectrogram_dis(data,points,nWindowLength,nOverlap,nFFT,fs)
% Use start and end point matrix to compute PSD for each overlapping
% segment of selected data while avoiding discontinuities.
% Input: [data vector, start and end point matrix, number of points of
%         window, number of points of overlap, number of points of FFT, 
%         sampling rate]
% Output: matrix with PSD for each window in each segment and time (x) and frequency (y) vectors 
% Example: [PSD_matrix,t,f]=spectrogram_dis(data_D_HP_LP,select_points,...
%          paramsspect.window,paramsspect.noverlap,paramsspect.nfft,Settings.fs);
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% Get number of segments
nSegment=size(points,1);

% Preallocate PSD_matrix
PSD_matrix=cell(nSegment,1);
t=PSD_matrix;

% Compute PSD of each window within each segment
for i=1:nSegment
    [~,f,t{i},PSD_matrix{i}]=spectrogram(data(points(i,1):points(i,2)),...
    nWindowLength, nOverlap, nFFT, fs);
    t{i}=t{i}(:);
end

% Concatenate PSD matrices of each segment into one PSD_matrix
PSD_matrix=[PSD_matrix{:}];

if nargout>1
    for i=1:nSegment
    t{i}=t{i}+points(i,1)/fs;
    end
    t=cat(1,t{:});
end

end

