function [data,data_t,data_lin] = LFP_Process_EF(EF_Name,Settings,Channels,varargin)
% LFP_PROCESS converts, resamples, detrends and filters data in an
% Experiment Folder (EF) created by OpenEphys.
% [DATA,DATA_T,DATA_LIN] = LFP_PROCESS_EF(EF_NAME,SETTINGS,CHANNELS,RESAMPLE)
% Processes EF with the name EF_NAME and SETTINGS and CHANNELS defined at the start of
% LFP_Analyzer_EFF. RESAMPLE input can be added optionally. If RESAMPLE is false no 
% resampling is done on LFP traces. Possibly returns 3 outputs: DATA is
% LFP trace data of all channels in a cell array, DATA_T is the corresponding time vector
% and DATA_LIN is the linear fit of all channels in a cell array.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
%% Go into folder, check inputs and preallocate
EphysFiles=dir([EF_Name,filesep,'*_CH*.continuous']); % get all .continuous files (channels)
nChannel=size(Channels,1);

% Check and parse inputs
narginchk(3,4)
if nargin<4
    Resample=true;
else
    Resample=varargin{1};
end

if nChannel~=length(EphysFiles)
    error('Defined number of channels does not match number of channels in experiment folder')
end

% Preallocate data cells
data=cell(nChannel,1);
info=data;
data_lin=data;

%% Load and process data
for i=1:nChannel
    % Load OpenEhys .continuous files
    [data{i},~,info{i}] = load_open_ephys_data_faster([EF_Name,filesep,EphysFiles(i).name]);
end

% get ends of experiment
data_t_end=zeros(nChannel,1);
for i=1:nChannel
    data_t_end(i)=length(data{i});
end
data_t_end_min=min(data_t_end);

% truncate data to make sure that they have the same lengths
for i=1:nChannel
    data_t_end_trunc=data_t_end(i)-data_t_end_min;
    data{i}(end-data_t_end_trunc+1:end)=[];
end

for i=1:nChannel
    
    if Resample
        % Resample data 
        data{i} = resample(data{i}, Settings.fs, info{i}.header.sampleRate);
    else
        Settings.fs=info{i}.header.sampleRate;
    end

    if Settings.Detrend
    % Detrend data for resampled data
    n_w=round(Settings.window_detrend*Settings.fs); %number of points of window
    n_o=round(n_w*Settings.overlap_detrend); %number of points of overlap
    data_lin{i}=LineFit(data{i},n_w,n_o); % get linear fit
    data{i}=data{i}-data_lin{i}; %subtract fit
    end
    
    % Filter
    data{i} = HighPassfilter(2, Settings.HP, Settings.fs, data{i});
    data{i} = LowPassfilter(2, Settings.LP, Settings.fs, data{i});
    
end

% One time vector for both channels
data_t=(0:1/Settings.fs:(length(data{1})-1)/Settings.fs)';

end