% Plot traces, spectrograms and PSDs per day for an Experiment Folder
% Folder (EFF) which is a folder containing experiment folders (e.g. '18-07-15_Group9_Wheel').
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Select and go to Experiment Folder Folder (EFF)
EFF = uigetdir('','Select folder containing experiment folders');
cd(EFF)

%% Plot settings
% plot property 'Color'
c=colormap(lines); % colors to most contrasting color map
% c(1:3,:)=c([2,3,1],:); %switch around colors depending on genotype
% c(1:2,:)=c([2,1],:);
close(gcf)

% plot property 'LineStyle'
ls{1}='-';
ls{2}='--';
ls{3}=':';
ls{4}='-.';

% plot properties 'LineWidth' and 'MarkerSize'
lw=1.5;

% misc plot settings
maxFreq=15;
normFreq=50;

%% Unselected data and normalized spectrogram plots per EFF
% method
temp.split=strsplit(EFF,filesep);
EFF=temp.split{end};
temp.split=strsplit(EFF,'_');
temp.method=temp.split{end};

% load settings and channel names for experiment folder folder
load('Channels')
nChannel=size(Channels,1);
load('Settings')

% get names of Experiment Folders (EF)
EF=dir('*_*_*');
nEF=length(EF);

% get names and days for each experiment folder
experiment.Name=cell(nEF,1);
experiment.Day=cell(nEF,1);
for ii=1:nEF
    temp.split=strsplit(EF(ii).name,'_');
    experiment.Name{ii}=temp.split{1};
    experiment.Day{ii}=temp.split{end};
end
uniq.Name=unique(experiment.Name);
uniq.Day=unique(experiment.Day);
nName=length(uniq.Name);
nDay=length(uniq.Day);

% FFT pre-processing
nWindow = Settings.fs * Settings.window;
nOverlap = nWindow * Settings.overlap;
nFFT = nWindow;

% predefine figures
for iii=1:nChannel
    fig.LFPs{iii}=figure('Name',['Unselected data ',EFF]);
    fig.norm_spec{iii}=figure('Name',['Normalized spectrogram ',EFF]);
end
for iii=1:nName
    fig.name_plot{iii}=figure('Name',['Normalized spectrogram ',EFF]);
    h{iii}={};
    legendcell{iii}={};
end

% preallocation
PSD=cell(nChannel,1);

for ii=1:nEF
    % get indices
    ind.Name=find(strcmp(uniq.Name,experiment.Name{ii}));
    ind.Day=find(strcmp(uniq.Day,experiment.Day{ii}));

    % go into experiment folder and load data 
    cd(EF(ii).name)
    load('EF_data')

    % time stamps for each channel should be the same
    data_t=0:1/Settings.fs:(length(data{1})-1)/Settings.fs;

    for iii=1:nChannel

        %%% Beginning of plot of unselected/selected LFPs traces %%%
        % get number of selected segments
        nSegment=size(select_points{iii},1);

        % plot
        set(0,'CurrentFigure',fig.LFPs{iii})
        subplot(nName,nDay,(ind.Name-1)*nDay+ind.Day), hold on
        temp.h(1)=plot(data_t,data{iii},'Color',c(2,:)); %unselected
        for iiii=1:nSegment
            temp.index=select_points{iii}(iiii,1):select_points{iii}(iiii,2);
            plot(data_t(temp.index),data{iii}(temp.index),...
                'Color',c(1,:)) %selected
        end

        % plot title for accepted data
        temp.time_length=sum(select_points{iii}(:,2)-select_points{iii}(:,1)+1)/Settings.fs;
        if temp.time_length>=Settings.time_threshold*60
        title([experiment.Name{ii},'\_',experiment.Day{ii}])
        % plot title for rejected data
        else
        title([experiment.Name{ii},'\_',experiment.Day{ii}],'Color',[1,0,0])
        end

%             xlabel('Time (s)'),ylabel('LFPs (\muV)')
        temp.h(2)=plot(NaN,NaN,'Color',c(1,:));
%             legend(temp.h,["unselected","selected"])
        ylim([-500,500])
        
        %%% Beginning of plot of normalized (Settings.HP-maxFreqHz) spectrogram %%%
        % normalized spectrogram calculations
        [~,temp.fff,temp.ttt,temp.ps]=spectrogram(data{iii},...
            nWindow, nOverlap,nFFT, Settings.fs);

        % plot
        set(0,'CurrentFigure',fig.norm_spec{iii})
        subplot(nName,nDay,(ind.Name-1)*nDay+ind.Day)  
        temp.logical=temp.fff>Settings.HP&temp.fff<=maxFreq;
        temp.ps_norm=(temp.ps(temp.logical,:)./sum(temp.ps(temp.logical,:)))*100; %normalized power spectrum with frequencies of interest (%)
        imagesc(temp.ttt,temp.fff(temp.logical),temp.ps_norm)
        set(gca,'YDir','normal')
        caxis([0,mean(max(temp.ps_norm))])
        colorbar              
        
        % plot title for accepted data
        temp.time_length=sum(select_points{iii}(:,2)-select_points{iii}(:,1)+1)/Settings.fs;
        if temp.time_length>=Settings.time_threshold*60
        title([experiment.Name{ii},'\_',experiment.Day{ii}])
        % plot title for rejected data
        else
        title([experiment.Name{ii},'\_',experiment.Day{ii}],'Color',[1,0,0])
        end
            
%       xlabel('Time (s)'),ylabel('LFPs PSDs (%)')

    PSD{iii}=mean(temp.ps,2);

    end
    
    %%% Beginning of plot of LFP PSDs per day %%%
    % Figure for each name, subplots for each channel
    set(0,'CurrentFigure',fig.name_plot{ind.Name})
    set(fig.name_plot{ind.Name},'Name',uniq.Name{ind.Name})
    for iii=1:nChannel
        subplot(1,nChannel,iii), hold on
        plot(temp.fff,PSD{iii},'Color',c(ind.Day,:),'LineStyle',ls{mod(ind.Day,4)+1},'LineWidth',lw);
        
    if iii==nChannel
        if temp.time_length<Settings.time_threshold*60
            hh=plot(NaN,'Color',c(ind.Day,:),'LineStyle',ls{mod(ind.Day,4)+1},'LineWidth',lw,...
                        'Marker','x','MarkerSize',5*lw,'MarkerEdgeColor',[1,0,0]);
        else
            hh=plot(NaN,'Color',c(ind.Day,:),'LineStyle',ls{mod(ind.Day,4)+1},'LineWidth',lw);
        end
        h{ind.Name}=[h{ind.Name},hh];
        legendcell{ind.Name}=[legendcell{ind.Name},experiment.Day{ii}];
    end
        
    end
    
    
    cd ..
end

for ii=1:nName
    set(0,'CurrentFigure',fig.name_plot{ii})
    for iii=1:nChannel
        subplot(1,nChannel,iii)
        grid on
        title(Channels{iii,1})
        xlabel('Frequency (Hz)'),ylabel('PSD (\muV^2)')
        xlim([Settings.HP,maxFreq])
        if iii==nChannel
            legend(h{ii},legendcell{ii})
        end
    end
end