% Format converter for converting Experiment Folders (EFs) in an Experiment
% Folder Folder (EFF) to a .mat format, inside of which there is a
% matformat structure with the names of channels as fields. These files can
% be used by Parade for analysis. Can only be run after LFP_Analyzer_EFF
% has been used, and therefore does processing (local detrending,
% filtering) by this script. Optionally LFP_ParadeConverter_unprocessed_EFF only
% applies resampling, so that it processing can be done by Parade instead
% of LFP_Analyzer_EFF.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all

% EFF and Parade folder with Output folder select
EFF = uigetdir('','Select folder containing experiment folders');
paradeFolder = uigetdir('','Select Parade folder to which converted files are saved');
cd (EFF)
EF = dir('*_*_*');
nEF=length(EF);

% Channel names
load('Channels.mat')
load('Settings.mat')

bar = waitbar(0, 'Converting data, please wait');
for i = 1:nEF

    waitbar((i-1)/length(EF), bar, 'Converting data, please wait'); 
    cd (EF(i).name)
    
    % Mouse code
    Mouse = strsplit(EF(i).name, '_');
    Mouse_Day=Mouse{4};
    Mouse_Name=['MI',Mouse{1}];
    Mouse_Name(Mouse_Name=='-')='_';
    Mouse = [Mouse_Name, '_', Mouse_Day];
    
    load('EF_data')
    
    for ii = 1:size(Channels,1)
        
        % Create input file
        matformat.(Channels{ii,1}) = data{ii};
        
    end
    
    % Go to Parade output folder and save data
    save([paradeFolder,filesep,Mouse,'_filtered'], 'matformat')
    clear matformat
    cd ..
    
end
close(bar);
msgbox('Conversion Complete')