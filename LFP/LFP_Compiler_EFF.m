% Use to compile all the analyzed data in Experiment Folders (EFs)
% into one file (DataTable.mat) in the Experiment Folder Folder (EFF).
% Compiled data includes spectral information of LFP traces in the form of
% Power Spectral Densities (PSDs).
%
% Marcel van Velze 19.09.2017
% marcelbvv@gmail.com
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all

warning('off', 'all')

% Select experiment folder and create list of experiment names
EFF = uigetdir('','Select folder containing experiments');
cd (EFF)
EF = dir('*_*_*');
nEF=length(EF);

% Get measurement type (end of name of F)
Measurement=strsplit(EFF,'_');
Measurement_Method=convertCharsToStrings(Measurement{end});

% Load Channels and Settings defined by Analyzer
load('Channels.mat')
nChannel=size(Channels,1);
load('Settings.mat')

% Settings
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

% Pre-allocation
nTableRow=nEF*nChannel;
temp.string_pre=strings(nTableRow,1);
temp.cell_pre=cell(nTableRow,1);
Measurement_Method_string=temp.string_pre;
Channel_string=temp.string_pre;
Mouse_Name_string=temp.string_pre;
Mouse_Day_string=temp.string_pre;
data_cell=temp.cell_pre;
select_logical_cell=temp.cell_pre;
time_cell=temp.cell_pre;
Power_cell=temp.cell_pre;
Power_Alt_cell=temp.cell_pre;
Settings_cell=temp.cell_pre;
fff_cell=temp.cell_pre;
            
% FFT pre-processing
nWindow=Settings.fs*Settings.window;
nOverlap=round(nWindow*Settings.overlap);
nFFT=nWindow;

bar = waitbar(0, 'Compiling Data, please wait');
count=0;
for i = 1:nEF
    
    waitbar((i-1)/nEF, bar, 'Compiling Data, please wait');
    cd (EF(i).name)
    
    % Load Experiment Folder (EF) data
    load ('EF_data')
    
    % Get name and day of measurement of mouse
    Mouse = strsplit(EF(i).name, '_');
    Mouse_Day=Mouse{4};
    Mouse_Name=Mouse{1};
    Mouse = [Mouse_Name, '_', Mouse_Day];
    
    for ii = 1:nChannel
        
        count=count+1;
    
        % Calculate Power Spectrum
        [PSD_eachwindow,~,fff]=spectrogram_dis(data{ii},select_points{ii},...
                                       nWindow, nOverlap, nFFT, Settings.fs);
        PSD=mean(PSD_eachwindow,2);
        
        % Code for saving excel file which contains only the relevant data (1-100Hz)
        fff_logical=fff>=1&fff<=100;
        Frequency = ['Frequency'; num2cell(fff(fff_logical))];
        Power_Channel = [Mouse; num2cell(PSD(fff_logical))];
        xlswrite(Channels{ii,1}, [Power_Channel,Frequency]); %Save 1-100Hz bins to excel file
    
        % Code for DataTable which contains all data
        Measurement_Method_string(count)=Measurement_Method;
        Channel_string(count)=convertCharsToStrings(Channels{ii,1});
        Mouse_Name_string(count)=convertCharsToStrings(Mouse_Name);
        Mouse_Day_string(count)=convertCharsToStrings(Mouse_Day);
        time_cell{count}=sum((select_points{ii}(:,2)-select_points{ii}(:,1)+1))/Settings.fs;
        Power_cell{count}=PSD;
        Settings_cell{count}=Settings;
        fff_cell{count}=fff;
        
    end
    cd ..
end

% Mouse Name
uniq.Name=unique(Mouse_Name_string);
nName=length(uniq.Name);

% Default genotype answer
Genotype_def=strings(nName,1);
if exist('DataTable.mat','file')
    load('DataTable')
    for i=1:nName
        Genotype_def(i)=DataTable.Genotype(find(strcmp(DataTable.Name,uniq.Name{i}),1,'first'));
    end
else
  
% Genotype_definput="";
% Genotype_definput="unknown";
% Genotype_definput="WT_129";
Genotype_definput="WT_C57";
Genotype_def(:)=Genotype_definput;

end

% Genotype input
uniq.Genotype=inputdlg({uniq.Name{:}}','Enter Genotype',[1,50],{Genotype_def{:}}',options);
Genotype_string=strings(nTableRow,1);
for i=1:nName
    Genotype_string(strcmp(Mouse_Name_string,uniq.Name{i}))...
    =convertCharsToStrings(uniq.Genotype{i});
end

% Default targeted answer
Targeted_def=strings(nName,1);
Targeted_def_temp=strings(nName,1);
Targeted_definput='no';
Targeted_def_temp{1}=Targeted_definput;
nChannel_full=length(Channels);
for i=2:nChannel_full
    Targeted_def_temp{i}=['_',Targeted_definput];
end
Targeted_def(:)=convertCharsToStrings([Targeted_def_temp{:}]);

channelInputWindowTitle=cell(nChannel,1);
channelInputWindowTitle(1)=Channels(1);
for i=2:nChannel
    channelInputWindowTitle{i}=['_',Channels{i}];
end
channelInputWindowTitle=cat(2,channelInputWindowTitle{:});
% Define Targeted per Mouse Name (Channels{1}_Channels{2} == yes_no)
uniq.Targeted=inputdlg({uniq.Name{:}}',['Enter Targeted "',channelInputWindowTitle,'"'],[1,100],{Targeted_def{:}}',options);

Targeted_string=strings(nTableRow,1);
for i=1:nName
    tempstring=strsplit(uniq.Targeted{i},'_');
    logical.name=strcmp(Mouse_Name_string,uniq.Name(i));
    uniq.Channel=unique(Channel_string(logical.name));
    nChannel=length(uniq.Channel);
    for ii=1:nChannel
        tempindex=find(strcmp(Channels,uniq.Channel(ii)));
        logical.name_channel=logical.name&strcmp(Channel_string,uniq.Channel(ii));
        Targeted_string(logical.name_channel)...
        =convertCharsToStrings(tempstring{tempindex});
    end
end

% Save DataTable to .mat file in experiment folder
TableNames={'Genotype','Method','Channel','Name','Day','Targeted','Time_Length','Power','Settings','Power_Bins'};
DataTable=table(Genotype_string,Measurement_Method_string,Channel_string,Mouse_Name_string,Mouse_Day_string,Targeted_string,...
                time_cell,Power_cell,Settings_cell,fff_cell,...
                'VariableNames',TableNames);
save('DataTable.mat', 'DataTable');

close(bar);
msgbox('Compiling Complete')