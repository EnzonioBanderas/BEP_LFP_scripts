% Use to re-analyze the data.
% Place Experiment Folders (EFs) in one Experiment Folder Folder (EFF), select that folder.
% Settings can be changed by altering answers to input windows.
% Same as LFP_Anlayzer_EFF but no selection of LFP trace is required.
% Optionally, LFP trace can be reselected if under the time threshold by
% answering the input window question 'Reselect if selection time is under threshold?'
% with 'yes'.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

clearvars, close all

%% Initial input
    
% Select experiment folder folder and create list of experiment folder names
EFF = uigetdir('','Select folder containing experiments');
cd (EFF)

% Define or load previously created Settings and Channels 
if exist([pwd,filesep,'Settings.mat'],'file')&&exist([pwd,filesep,'Channels.mat'],'file') % both Settings.mat and Channels.mat exist
    DefineSettings = strcmp(inputdlg('Define Settings','Define Settings?',1,"no"),'yes');
else
    DefineSettings=true;
end

% Define Channels
if DefineSettings % if yes input was given for Define Settings input window or if there was no Settings.mat or Channels.mat in EFF
    
% Define number of channels
nChannel = str2double(cell2mat(inputdlg('Enter number of channels:','Number of channels',1,{'1'})));
prompt_channels=cell(nChannel,1);

% prompt names
for i=1:nChannel
    prompt_channels(i)={['Channel ',num2str(i)]};
end

% Define channel names
Channels = inputdlg(prompt_channels,'Names of channels',1);

% Save channel names
save('Channels','Channels')
else % load previously defined Channels
    load('Channels')
    nChannel=size(Channels,1);
end

% Define Settings
if DefineSettings % if yes input was given for Define Settings input window
% Settings: prompt window
SettingsCell={'Sampling rate (Samples/s)','High-pass cutoff frequency (Hz)',...
              'Low-pass cutoff frequency (Hz)','Window size (s)',...
              'Window overlap (ratio)','Time threshold for data acceptance (min)',...
              'Detrend data?','Detrend window overlap (s)','Detrend overlap (ratio)',...
              'Plot linear fit to data when selecting?',...
              'Select traces for multiple channels simultaneously?',...
              'Select experiment folders manually?',...
              'Randomize order of appearance of mice?'...
              'Reselect if selection time is under threshold?'};
SettingsDef={'1000','1','100','2','0.5','5','no','1','0.5','no','yes','no','no','no'};
Settings.input=inputdlg(SettingsCell,'Enter settings',1,SettingsDef);

% Settings: Assign input
Settings.fs = str2double(Settings.input{1});
Settings.HP = str2double(Settings.input{2});
Settings.LP = str2double(Settings.input{3});
Settings.window = str2double(Settings.input{4});
Settings.overlap = str2double(Settings.input{5});
Settings.time_threshold = str2double(Settings.input{6}); 
Settings.Detrend = strcmp(Settings.input{7},'yes');
Settings.window_detrend = str2double(Settings.input{8}); 
Settings.overlap_detrend = str2double(Settings.input{9});
Settings.ShowLinearFit = strcmp(Settings.input{10},'yes');
Settings.SelectMultipleChannels = strcmp(Settings.input{11},'yes');
Settings.SelectExperimentFolders = strcmp(Settings.input{12},'yes');
Settings.Randomize = strcmp(Settings.input{13},'yes');
Settings.Reselect = strcmp(Settings.input{14},'yes');

% Save settings
save('Settings','Settings')

else % load previously defined Settings
    load('Settings')
end

if Settings.SelectExperimentFolders
% Select experiment folders
EF_Select = uiselect([],'Select experiment folders to analyse');
nEF=length(EF_Select);
EF_Name=cell(nEF,1);
for i=1:nEF
    % get filenames
    temp.FileName=strsplit(EF_Select{i},filesep);
    EF_Name(i)=temp.FileName(end);
end

else
    
EF = dir('*_*_*'); %added so that only mouse-day folder names are included in directory
nEF=length(EF);
EF_Name={EF(:).name}';

end

% Preallocate
data_EF=cell(nEF,1);
data_lin_EF=data_EF;
select_points_EF=cell(nEF,1);

%% Process data
bar = waitbar(0, 'Processing data, please wait');
for i=1:nEF
waitbar((i-1)/nEF, bar, 'Processing data, please wait');
[data_EF{i},~,data_lin_EF{i}]=...
    LFP_Process_EF(EF_Name{i},Settings,Channels);
end
close(bar)

%% Load old selection data
select_logical=cell(nChannel,1);
nWindow=Settings.fs*Settings.window;
for i=1:nEF
    cd(EF_Name{i})
    load('EF_data','select_points')
    select_points_EF{i}=select_points;
    for ii=1:nChannel
    select_logical{ii}=points2logical(select_points{ii});
    end
    % truncate
    select_logical_length=zeros(nChannel,1);
    for ii=1:nChannel
    select_logical_length(ii)=length(select_logical{ii});
    end
    for ii=1:nChannel
    truncation_length=select_logical_length(ii)-min(select_logical_length);
    select_logical{ii}(end-truncation_length+1:end)=[];
    end
    % combine
    for ii=2:nChannel
    select_logical{ii}=select_logical{ii}&select_logical{ii-1};
    end
    select_points_temp=logical2points(select_logical{nChannel});
    select_points_temp=select_points_temp(select_points_temp(:,2)-select_points_temp(:,1)+1>=nWindow,:);
    select_time=sum((select_points_temp(:,2)-select_points_temp(:,1)+1)/Settings.fs);
    
    if (isempty(select_points_temp)||select_time<Settings.time_threshold*60)&&Settings.Reselect
        
        data_t=(0:1/Settings.fs:(length(data_EF{i}{1})-1)/Settings.fs)';
        select_points=LFP_Select_EF(data_EF{i},data_t,data_lin_EF{i},...
                            EF_Name{i},Settings,Channels);
    else
        select_points=cell(nChannel,1);
        for ii=1:nChannel
            select_points{ii}=select_points_temp;
        end
        
    end
    select_points_EF{i}=select_points;
    cd ..
end
close(gcf)

%% Save data
for i=1:nEF
    cd(EF_Name{i})
    data=data_EF{i};
    select_points=select_points_EF{i};
    save('EF_data','data','select_points')
    cd ..
end

msgbox('Re-Analysis Complete')