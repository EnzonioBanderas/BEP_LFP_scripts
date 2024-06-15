% Use to analyze the data.
% Place Experiment Folders (EFs) in one Experiment Folder Folder (EFF), select that folder.
% Settings can be changed by altering answers to input windows.
%
% Marcel van Velze 19.09.2017
% marcelbvv@gmail.com
% Enzo Nio Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Initial input

% Select experiment folder folder and create list of experiment folder names
EF = uigetdir('','Select folder containing experiments');
cd (EF)

% Define or load previously created Settings and Channels 
if exist([pwd,filesep,'Settings.mat'],'file')&&exist([pwd,filesep,'Channels.mat'],'file') % both Settings.mat and Channels.mat exist
    DefineSettings = strcmp(inputdlg('Define Settings','Define Settings?',1,"no"),'yes');
else
    DefineSettings=true;
end

% Define Channels
if DefineSettings % if yes input was given for Define Settings input window
    
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
              'Randomize order of appearance of mice?'};
SettingsDef={'1000','1','100','2','0.5','5','yes','1','0.5','no','yes','no','no'};
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
data_experiment=cell(nEF,1);
data_t_experiment=data_experiment;
data_lin_experiment=data_experiment;
select_points_experiment=cell(nEF,1);

%% Process data
bar = waitbar(0, 'Processing data, please wait');
for i=1:nEF
waitbar((i-1)/nEF, bar, 'Processing data, please wait');
[data_experiment{i},data_t_experiment{i},data_lin_experiment{i}]=LFP_Process_EF(EF_Name{i},Settings,Channels);
end
close(bar)

if Settings.Randomize
%% Randomize experiment order
ExperimentNameOld=EF_Name;
% get names
Name=cell(1,nEF);
for i = 1:nEF
    tempsplit = strsplit(EF_Name{i}, '_');
    Name{i} = tempsplit{1};
end

% get unique names
[unique_Name,~,unique_index]=unique(Name);
nName=length(unique_Name);

% define random permutation
rand_switch=randperm(nName); 

% merge data according to name of experiment 
EF_Name_name=cell(nName,1);
data_name=EF_Name_name;
data_t_name=EF_Name_name;
data_lin_name=EF_Name_name;
index_name=EF_Name_name;
for i=1:nName
    EF_Name_name{i}=EF_Name(unique_index==i);
    data_name{i}=data_experiment(unique_index==i);
    data_t_name{i}=data_t_experiment(unique_index==i);
    data_lin_name{i}=data_lin_experiment(unique_index==i);
    index_name{i}=ones(sum(unique_index==rand_switch(i)),1)*i;
end

% randomize merged data
EF_Name_name=EF_Name_name(rand_switch);
data_name=data_name(rand_switch);
data_t_name=data_t_name(rand_switch);
data_lin_name=data_lin_name(rand_switch);

% unpack data
EF_Name=cat(1,EF_Name_name{:});
data_experiment=cat(1,data_name{:});
data_t_experiment=cat(1,data_t_name{:});
data_lin_experiment=cat(1,data_lin_name{:});
index_experiment=cat(1,index_name{:});

% The named that are showed during selection
EF_Name_forselect=cell(nEF,1);
for i=1:nEF
    EF_Name_forselect{i}=['Mouse ',num2str(index_experiment(i))]; 
end

else
    EF_Name_forselect=EF_Name;
end

%% Select data
for i=1:nEF
    select_points_experiment{i} = LFP_Select_EF(data_experiment{i},data_t_experiment{i},data_lin_experiment{i},...
                            EF_Name_forselect{i},Settings,Channels);
end
close(gcf)

%% Save data
for i=1:nEF
    cd(EF_Name{i})
    data=data_experiment{i};
    data_t=data_t_experiment{i};
    select_points=select_points_experiment{i};
    save('EF_data','data','select_points')
    cd ..
end

msgbox('Analysis Complete')