% Plot various graphs for an Experiment Folder Folder Folder (EFFF) which 
% is a folder containing folders containing experiment folders
% Assumption: settings not changed over EFFs
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Data import
% Select folder which contains folders which contain experiment folders (EFFF) and load DataTables
EFFF = uigetdir('','Select folder containing folders containing experiment folders');
cd (EFFF)
EFF = dir('*_*_*');
nEFF=length(EFF);

DataTable_cell=cell(nEFF,1);
for i=1:nEFF
    cd(EFF(i).name)
    load('DataTable')
    cd ..
    DataTable_cell{i}=DataTable;
end
DataTable=cat(1,DataTable_cell{:});

%%% Assume settings and therefore frequency bins to be constant over EFF's %%%
Settings=DataTable.Settings{1};
fff=DataTable.Power_Bins{1};

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

% plot less (somewhat unecessary) plots
plotlessCell={'Logarithmic plots?',...
              'Seperate Mean Name over Day plots?',...
              'Extra mean plots?',...
              'LFP plots?',...
              'Spectrograms?',...
              'Plot PSDs for each day',...
              'Targeted?'};
plotlessDef={'no','no','no','no','no','no','no'};
plotless=inputdlg(plotlessCell,'Enter yes or no',1,plotlessDef);
if strcmp(plotless{1},'no')
    plot_logarithmic=2;
else
    plot_logarithmic=4;
end
if strcmp(plotless{2},'no')
    plot_seperatemean=false;
else
    plot_seperatemean=true;
end
if strcmp(plotless{3},'no')
    plot_extramean=false;
else
    plot_extramean=true;
end
if strcmp(plotless{4},'no')
    plot_data=false;
else
    plot_data=true;
end
if strcmp(plotless{5},'no')
    plot_spectrogram=false;
else
    plot_spectrogram=true;
end
if strcmp(plotless{6},'no')
    plot_PSDeachday=false;
else
    plot_PSDeachday=true;
end
if strcmp(plotless{7},'no')
    plot_targeted=false;
else
    plot_targeted=true;
end

%% Brain wave definitions and misc preallocation
temp.Logical_Time=cell2mat(DataTable.Time_Length)>DataTable.Settings{1}.time_threshold*60;

delta=[2,4];
theta=[5,10];
beta=[13,30];
gamma=[30,50];
bandstr=["delta";"theta";"beta";"gamma";"total"];
nBand=length(bandstr);

%% Unselected data and normalized spectrogram plots per EFF
for i=1:nEFF
    % method
    temp.split=strsplit(EFF(i).name,'_');
    temp.method=temp.split{end};
    
    % go into folder with experiment folders and get names of experiment folders
    cd(EFF(i).name)
    EF=dir('*_*_*');
    nEF=length(EF);
    
    % load settings and channel names for experiment folder folder
    load('Channels')
    nChannel=size(Channels,1);
    load('Settings')
    
    % FFT pre-processing
    nWindow = Settings.fs * Settings.window;
    nOverlap = nWindow * Settings.overlap;
    nFFT = nWindow;
    
    % predefine figures
    for iii=1:nChannel
        if plot_data
        fig.LFPs{iii}=figure('Name',['Unselected data ',EFF(i).name]);
        end
        if plot_spectrogram
        fig.norm_spec{iii}=figure('Name',['Normalized spectrogram ',EFF(i).name]);
        end
    end
    
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
            if plot_data
            
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
            
            end

            if plot_spectrogram
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
            
%             xlabel('Time (s)'),ylabel('LFPs PSDs (%)')
            end
            
        end
        cd ..
    end
    cd ..
end

%% Name each Day plots for each method and each channel
if plot_PSDeachday

uniq.Day_Full=unique(DataTable.Day);
uniq.Method=unique(DataTable.Method);
for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
uniq.Channel=unique(DataTable.Channel(temp.Logical_i));
    for ii=1:length(uniq.Channel)
    temp.Logical_ii=strcmp(DataTable.Channel,uniq.Channel{ii})&temp.Logical_i;
    uniq.Name=unique(DataTable.Name(temp.Logical_ii));
    
        % for loop for 4 different figures
        for j=1:plot_logarithmic
        if j==1
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Each Day for Name -Absolute'])
        maxFreq=10;
        end
        if j==2
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Each Day for Name -Normalized'])
        maxFreq=10;
        end
        if j==3
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Each Day for Name -Absolute -Logarithmic'])
        maxFreq=50;
        end
        if j==4
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Each Day for Name -Normalized -Logarithmic'])
        maxFreq=50;
        end
    
        for iii=1:length(uniq.Name)
        subplot(1,length(uniq.Name),iii), hold on
        temp.Logical_iii=strcmp(DataTable.Name,uniq.Name{iii})&temp.Logical_ii;
        temp.Day=DataTable.Day(temp.Logical_iii);
        uniq.Day=unique(temp.Day);
        temp.Power=DataTable.Power(temp.Logical_iii);
        temp.Time_Length=DataTable.Time_Length(temp.Logical_iii);
    
        % Plot if
        if sum(temp.Logical_iii)~=0
    
        hhh=zeros(sum(temp.Logical_iii),1);
        for iiii=1:sum(temp.Logical_iii)
            
        % for normalized plots
        if (j==2||j==4)
            temp.Power{iiii}=(temp.Power{iiii}/bandpower2(temp.Power{iiii},fff,[Settings.HP,normFreq]))*100;
        end
        
        c_index=find(strcmp(uniq.Day_Full,temp.Day{iiii}));
        while c_index>size(c,1)
            c_index=c_index-size(c,1);
        end
        ls_index=find(strcmp(uniq.Day_Full,temp.Day{iiii}));
        while ls_index>length(ls)
            ls_index=ls_index-length(ls);
        end
        
        if temp.Time_Length{iiii}>=Settings.time_threshold*60
        plot(fff,temp.Power{iiii},'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw)
        hhh(iiii)=plot(NaN,NaN,'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw);
        else
        plot(fff,temp.Power{iiii},'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw)
        hhh(iiii)=plot(NaN,NaN,'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw,...
                             'Marker','x','MarkerSize',5*lw,'MarkerEdgeColor',[1,0,0]);
        end
            
        end
    
        % Figure settings
        legend(hhh, {temp.Day{:}},'Interpreter','none');
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        title(uniq.Name{iii})
        xlabel('Frequency (Hz)')
        if j==3||j==4 % for logarithmic plots
            set(gca, 'YScale', 'log')
        end
        if j==1||j==3 % for non-normalized plots
            ylabel('LFPs PSD (\muV^2)')
        else % for normalized plots
            ylabel('LFPs PSD (%)')
        end
        
    
        end
        end % end of 1 plot
        end % end of 4 plots
    
    end
end

end
%% Mean Name over Days plots
if plot_seperatemean

SEM_transparency=0.2;

uniq.Method=unique(DataTable.Method);
for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
uniq.Channel=unique(DataTable.Channel(temp.Logical_i));
    for ii=1:length(uniq.Channel)
    temp.Logical_ii=strcmp(DataTable.Channel,uniq.Channel{ii})&temp.Logical_i;
    uniq.Name=unique(DataTable.Name(temp.Logical_ii));
    
        % for loop for 4 different figures
        for j=1:plot_logarithmic
        if j==1
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Mean Name over Days -Absolute'])
        maxFreq=10;
        end
        if j==2
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Mean Name over Days -Normalized'])
        maxFreq=10;
        end
        if j==3
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Mean Name over Days -Absolute -Logarithmic'])
        maxFreq=50;
        end
        if j==4
        figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},...
               ' Mean Name over Days -Normalized -Logarithmic'])
        maxFreq=50;
        end
    
%         hh=zeros(length(uniq.Name),1);
        count=0;
        for iii=1:length(uniq.Name)
            
            temp.Logical_iii=strcmp(DataTable.Name,uniq.Name{iii})&temp.Logical_ii;
            temp.Logical_iii_Time=temp.Logical_iii&temp.Logical_Time;
        
            temp.Name=DataTable.Name(temp.Logical_iii_Time);
            temp.Power=DataTable.Power(temp.Logical_iii_Time);
            
            temp.Power_arr=cell2mat(temp.Power');
            temp.Power_Mean=mean(temp.Power_arr,2);
            temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
        
            if sum(temp.Logical_iii_Time)~=0
                count=count+1;
            
            % for normalized plots
            if (j==2||j==4)
                temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
                temp.Power_Mean=mean(temp.Power_arr,2);
                temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
            end
        
            % for the consistency of the color and style of lines
            c_index=find(strcmp(uniq.Name,temp.Name{1}));
            while c_index>size(c,1)
                c_index=c_index-size(c,1);
            end
            ls_index=find(strcmp(uniq.Name,temp.Name{1}));
            while ls_index>length(ls)
                ls_index=ls_index-length(ls);
            end
        
            % First figure in subplot with mean of each mouse
            subplot(1,length(uniq.Name)+1,1), hold on 
            hh(count)=plot(fff,temp.Power_Mean,...
                        'Color',c(c_index,:),...
                        'LineStyle',ls{ls_index},...
                        'LineWidth',lw);
            legend_name{count}=uniq.Name{iii};
            if j==1||j==2 %show SEM if non-logarithmic
            patch([fff;...
                   fff(end:-1:1);...
                   fff(1)],...
                  [temp.Power_Mean-temp.SEM;...
                   temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
                   temp.Power_Mean(1)-temp.SEM(1)],...
                   c(c_index,:),...
                   'FaceAlpha',SEM_transparency,...
                   'EdgeAlpha',SEM_transparency)
            end
        
            % Last figures (which put together make the first figure)
            subplot(1,length(uniq.Name)+1,iii+1), hold on
            h_single=plot(fff,temp.Power_Mean,...
                'Color',c(c_index,:),...
                'LineStyle',ls{ls_index},...
                'LineWidth',lw);
            if j==1||j==2 %show SEM if non-logarithmic
            patch([fff;...
                   fff(end:-1:1);...
                   fff(1)],...
                  [temp.Power_Mean-temp.SEM;...
                   temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
                   temp.Power_Mean(1)-temp.SEM(1)],...
                   c(c_index,:),...
                   'FaceAlpha',SEM_transparency,...
                   'EdgeAlpha',SEM_transparency)
            end
            
            % figure settings (not first)
            legend(h_single,{temp.Name{1}},'AutoUpdate','off','Interpreter','none')
            xlim([Settings.HP,maxFreq])
            xticks(linspace(0,maxFreq,6))
            title(uniq.Name{iii})
            xlabel('Frequency (Hz)')
            % if logarithmic plot
            if j==3||j==4
                set(gca, 'YScale', 'log')
            end    
            if j==1||j==3 % for non-normalized plots
                ylabel('LFPs PSD \muV^2')
            else % for normalized plots
                ylabel('LFPs PSD (%)')
            end
        
            end
        
        end % end of 1 plot
        
        % first figure settings
        subplot(1,length(uniq.Name)+1,1)
        legend(hh,legend_name,'Interpreter','none')
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        title('Mean Name over Days')
        xlabel('Frequency (Hz)')
        % if logarithmic plot
        if j==3||j==4
            set(gca, 'YScale', 'log')
        end
        if j==1||j==3 % for non-normalized plots
            ylabel('LFPs PSD (\muV^2)')
        else % for normalized plots
            ylabel('LFPs PSD (%)')
        end
    
        end % end of 4 plots 
    end 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean Each Name over Days (Channel-Method) 
SEM_transparency=0.2;
maxFreq=10;
count=0;

uniq.Name_Full=unique(DataTable.Name);
uniq.Channel_Full=unique(DataTable.Channel);
uniq.Method=unique(DataTable.Method);

for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
uniq.Channel=unique(DataTable.Channel(temp.Logical_i));
    for ii=1:length(uniq.Channel)
    temp.Logical_ii=strcmp(DataTable.Channel,uniq.Channel{ii})&temp.Logical_i;
    temp.Logical_ii_Time=temp.Logical_ii&temp.Logical_Time;
    uniq.Name=unique(DataTable.Name(temp.Logical_ii_Time));
    figure('Name',[uniq.Method{i},' ',uniq.Channel{ii},' Mean each Name over Days'])
        hh=zeros(length(uniq.Name),2);
        for iii=1:length(uniq.Name)
            
            temp.Logical_iii_Time=strcmp(DataTable.Name,uniq.Name{iii})&temp.Logical_ii_Time;
        
            temp.Name=DataTable.Name(temp.Logical_iii_Time);
            temp.Power=DataTable.Power(temp.Logical_iii_Time);
            
            temp.Power_arr=cell2mat(temp.Power');
            temp.Power_Mean=mean(temp.Power_arr,2); % save this variable
            temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
            
            % save every needed variable in a temporary 'table'
            count=count+1;
            temptable.Method(count)=uniq.Method(i);
            temptable.Channel(count)=uniq.Channel(ii);
            temptable.Name(count)=uniq.Name(iii);
            temptable.Power{count}=temp.Power_Mean;
            temptable.Genotype(count)=DataTable.Genotype(find(temp.Logical_iii_Time,1,'first'));
            temptable.Targeted(count)=DataTable.Targeted(find(temp.Logical_iii_Time,1,'first'));
            
            if sum(temp.Logical_iii_Time)~=0   
        
            % for the consistency of the color and style of lines
            c_index=find(strcmp(uniq.Name_Full,temp.Name{1}));
            while c_index>size(c,1)
                c_index=c_index-size(c,1);
            end
            ls_index=find(strcmp(uniq.Name_Full,temp.Name{1}));
            while ls_index>length(ls)
                ls_index=ls_index-length(ls);
            end
               
            for jj=1:2
                
            if jj==2 % Normalized plot
                temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
                temp.Power_Mean=mean(temp.Power_arr,2);
                temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
            end
                
            subplot(1,2,jj), hold on
            hh(iii,jj)=plot(fff,temp.Power_Mean,...
                'Color',c(c_index,:),...
                'LineStyle',ls{ls_index},...
                'LineWidth',lw);
            patch([fff;...
                   fff(end:-1:1);...
                   fff(1)],...
                  [temp.Power_Mean-temp.SEM;...
                   temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
                   temp.Power_Mean(1)-temp.SEM(1)],...
                   c(c_index,:),...
                   'FaceAlpha',SEM_transparency,...
                   'EdgeAlpha',SEM_transparency)
            end
            
            end
                
        end
            
        subplot(1,2,1)
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
        title([uniq.Channel{ii},' ',uniq.Method{i},' -Absolute'])
        legend(hh(:,1),uniq.Name,'Interpreter','none')
xticks('auto')
% set(gca,'Yscale','log')
grid on
        subplot(1,2,2)
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        xlabel('Frequency (Hz)'),ylabel('LFPs PSD(%)')
        title([uniq.Channel{ii},' ',uniq.Method{i},' -Normalized'])
        legend(hh(:,2),uniq.Name,'Interpreter','none')
xticks('auto')
% set(gca,'Yscale','log')
grid on
            
    end
    
end 

%% Mean Genotype-Channel-Method over Names and Days (norm and abs)
SEM_transparency=0.2;
maxFreq=20;

% Create 1 figure for each method and a cell within a cell for each subplot
uniq.Method_Full=unique(temptable.Method);
nMethod_Full=length(uniq.Method);
hh=cell(nMethod_Full,1);
legendCell=cell(nMethod_Full,1);
bandpower_cell=cell(nMethod_Full,1);
for i=1:nMethod_Full
fig.method(i)=figure('Name',['Mean Genotype-Channel for ',uniq.Method_Full{i}]);
hh{i}=cell(2,1);
legendCell{i}=cell(2,1);
bandpower_cell{i}=cell(2,1);
end

index=0;
count=0;
c_index=index;
ls_index=index;

uniq.Genotype=unique(temptable.Genotype);
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv"];
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv","STOPrescue-129/Sv"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:length(uniq.Genotype)
temp.Logical_g=strcmp(temptable.Genotype,uniq.Genotype{g});
uniq.Method=unique(temptable.Method(temp.Logical_g));
for i=1:length(uniq.Method)
    
% make appropiate method figure current
ind.method=find(strcmp(uniq.Method_Full,uniq.Method(i)));
set(0,'CurrentFigure',fig.method(ind.method))

temp.Logical_g_i=temp.Logical_g&strcmp(temptable.Method,uniq.Method{i});
uniq.Channel=unique(temptable.Channel(temp.Logical_g_i));
for ii=1:length(uniq.Channel)
temp.Logical_g_ii=temp.Logical_g_i&strcmp(temptable.Channel,uniq.Channel{ii});

temp.Power=temptable.Power(temp.Logical_g_ii);

temp.Power_arr=cell2mat(temp.Power);
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
if sum(temp.Logical_g_ii)~=0   
index=index+1;
for jj=1:2
        
    % no consistency because of 3 variables
    c_index=g;
    while c_index>size(c,1)
    c_index=c_index-size(c,1);
    end
    ls_index=ii;
    while ls_index>length(ls)
    ls_index=ls_index-length(ls);
    end
         
    if jj==2 % Normalized plot
        temp.Power_arr_norm=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr_norm,2);
        temp.SEM=std(temp.Power_arr_norm,0,2)/sqrt(size(temp.Power_arr_norm,2));
    end
         
    subplot(1,2,jj), hold on
    temp.hh=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    hh{i}{jj}=[hh{i}{jj},temp.hh];
    temp.legend={[uniq.Genotype{g},' ',uniq.Channel{ii},' ',uniq.Method{i}]};
    legendCell{i}{jj}=[legendCell{i}{jj},temp.legend];
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
 
end   
                     
% temporary 'table'
count=count+1;
bandtable.Genotype(count)=uniq.Genotype(g);
bandtable.Method(count)=uniq.Method(i);
bandtable.Channel(count)=uniq.Channel(ii);
bandtable.Power{count}=[bandpower2(temp.Power_arr,fff,delta);...
                        bandpower2(temp.Power_arr,fff,theta);...
                        bandpower2(temp.Power_arr,fff,beta);...
                        bandpower2(temp.Power_arr,fff,gamma);...
                        bandpower2(temp.Power_arr,fff)];
         
end
end

end
end

for i=1:nMethod_Full
if exist('hh')~=0
    
set(0,'CurrentFigure',fig.method(i))

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Genotype-Channel-Method -Absolute')
legend(hh{i}{1},legendCell{i}{1},'Interpreter','none')
grid on

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Genotype-Channel-Method -Normalized')
legend(hh{i}{2},legendCell{i}{2},'Interpreter','none')
grid on

end
end

%% Histograms of frequency band powers
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% figure for each method-channel
% subplots with each genotype-band for both absolute-normalized
uniq.Channel_Full=unique(bandtable.Channel);
nChannel_full=length(uniq.Channel_Full);
for ii=1:nChannel_full
    fig.Channel(ii)=figure('Name',['Band powers for ',uniq.Channel_Full{ii}]);
end

Genotype_bandtable=unique(bandtable.Genotype);
nGenotype_bandtable=length(Genotype_bandtable);

uniq.Method=unique(bandtable.Method);
nMethod=length(uniq.Method);
power_cellmatrix=cell(nMethod,nChannel_full,nGenotype_bandtable,nGenotype_bandtable,nBand);
for i=1:nMethod
logical.i=strcmp(bandtable.Method,uniq.Method(i));
uniq.Channel=unique(bandtable.Channel(logical.i));
nChannel=length(uniq.Channel);
for ii=1:nChannel
set(0,'CurrentFigure',fig.Channel(ii))
logical.ii=logical.i&strcmp(bandtable.Channel,uniq.Channel(ii));
uniq.Genotype=unique(bandtable.Genotype(logical.ii));
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv"]; 
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv","STOPrescue-129/Sv"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGenotype=length(uniq.Genotype);

%%% define bin location
temp.step2=bbd+bw;
temp.step1=bw+(nGenotype-1)*(temp.step2)+bgd;
select_points=(bgd:temp.step1:bgd+(nBand-1)*temp.step1)';
bin_matrix=zeros(nBand,nGenotype);
for j=1:nBand
    bin_matrix(j,:)=select_points(j)+bw/2:temp.step2:select_points(j)+bw/2+(nGenotype-1)*temp.step2;
end

for g=1:nGenotype
logical.ii_g=logical.ii&strcmp(bandtable.Genotype,uniq.Genotype(g));
for jj=1:2
    
%%% start plot
subplot(nMethod,2,(i-1)*nMethod+jj), hold on, grid on
for j=1:nBand
    if jj==1
        temp.power=bandtable.Power{logical.ii_g}(j,:);
    else %normalize
        temp.power=(bandtable.Power{logical.ii_g}(j,:)./bandtable.Power{logical.ii_g}(end,:))*100;
    end
    ind.ii=find(strcmp(uniq.Channel_Full,uniq.Channel(ii)));
    ind.g=find(strcmp(uniq.Genotype(g),Genotype_bandtable));
    power_cellmatrix{i,ind.ii,ind.g,j,jj}=temp.power;
    temp.power_mean=mean(temp.power);
    temp.SEM=std(temp.power)/sqrt(length(temp.power));
    bar(bin_matrix(j,g),temp.power_mean,...
    'FaceColor',c(g,:),'BarWidth',bw);
    errorbar(bin_matrix(j,g),temp.power_mean,temp.SEM,...
    'Color',c(g,:),'LineStyle','none','LineWidth',lw)
end

%xticks
xticks(select_points+(bw+(nGenotype-1)*temp.step2)/2)
xticklabels(bandstr)

%legend
for gg=1:nGenotype
    h_bin(gg)=bar(NaN,'FaceColor',c(gg,:),'BarWidth',bw);
end

if j==1
legend(h_bin,uniq.Genotype,'Interpreter','none','AutoUpdate','off')
end

%title and ylabel
if jj==1
title([uniq.Channel{ii},' ',uniq.Method{i},' -absolute'])
ylabel('Power (\muV^2)')
else
title([uniq.Channel{ii},' ',uniq.Method{i},' -normalized'])
ylabel('Power (%)')
end
%%%% end plot

end
end
end
end

%Save statistics to table
count=0;
for i=1:nMethod 
for ii=1:nChannel_full
for j=1:nBand
for jj=1:2
    
    if jj==1
        normstr="Absolute";
    else
        normstr="Normalized";
    end
    
H=cell(nGenotype_bandtable); 
P=H;
CI=H;
STATS=H;
for g=1:nGenotype_bandtable
for gg=1:nGenotype_bandtable
    
    if ~isempty(power_cellmatrix{i,ii,g,j,jj})&&~isempty(power_cellmatrix{i,ii,gg,j,jj})
    [H{g,gg},P{g,gg},CI{g,gg},STATS{g,gg}]...
    =ttest2(power_cellmatrix{i,ii,g,j,jj},power_cellmatrix{i,ii,gg,j,jj},...
    'vartype','unequal');
    end
    
end
end

    count=count+1;
    
    bandtable2.method(count,1)=uniq.Method(i);
    bandtable2.channel(count,1)=uniq.Channel_Full(ii);
    bandtable2.normalization(count,1)=normstr;
    bandtable2.band(count,1)=bandstr(j);
    
    temp.charH='H(:,1)';
    temp.charP='P(:,1)';
    temp.charCI='CI(:,1)';
    temp.charSTATS='STATS(:,1)';
    for ggg=2:nGenotype_bandtable
        temp.charH=[temp.charH,',H(:,',num2str(ggg),')'];
        temp.charP=[temp.charP,',P(:,',num2str(ggg),')'];
        temp.charCI=[temp.charCI,',CI(:,',num2str(ggg),')'];
        temp.charSTATS=[temp.charSTATS,',STATS(:,',num2str(ggg),')'];
    end
    bandtable2.H{count,1}=eval(['table(',temp.charH,')']);
    bandtable2.H{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.H{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.P{count,1}=eval(['table(',temp.charP,')']);
    bandtable2.P{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.P{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.CI{count,1}=eval(['table(',temp.charCI,')']);
    bandtable2.CI{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.CI{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.STATS{count,1}=eval(['table(',temp.charSTATS,')']);
    bandtable2.STATS{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.STATS{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.P_matrix{count,1}=cell2mat(P);
    
end
end
end
end
TableNames={'Method','Channel','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix'};
BandTable=table(bandtable2.method,bandtable2.channel,bandtable2.normalization,bandtable2.band,...
                bandtable2.H,bandtable2.P,bandtable2.CI,bandtable2.STATS,bandtable2.P_matrix,...
                'VariableNames',TableNames);
save('BandTable.mat', 'BandTable');

if plot_targeted
%% Targeted for each mice
clearvars temptable
SEM_transparency=0.2;
maxFreq=10;
count=0;

uniq.Name_Full=unique(DataTable.Name);
uniq.Targeted_Full=unique(DataTable.Targeted);
uniq.Method=unique(DataTable.Method);

for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
uniq.Targeted=unique(DataTable.Targeted(temp.Logical_i));
    for ii=1:length(uniq.Targeted)
    temp.Logical_ii=strcmp(DataTable.Targeted,uniq.Targeted{ii})&temp.Logical_i;
    temp.Logical_ii_Time=temp.Logical_ii&temp.Logical_Time;
    uniq.Channel=unique(DataTable.Channel(temp.Logical_ii_Time));
    for ch=1:length(uniq.Channel)
    figure('Name',[uniq.Method{i},' ',uniq.Targeted{ii},' ',uniq.Channel{ch},' Mean each Name over Days'])
    temp.Logical_ii_Time_ch=temp.Logical_ii_Time&strcmp(DataTable.Channel,uniq.Channel{ch});
        uniq.Name=unique(DataTable.Name(temp.Logical_ii_Time_ch));
        hh=zeros(length(uniq.Name),2);
        for iii=1:length(uniq.Name)
            
            temp.Logical_iii_Time_ch=strcmp(DataTable.Name,uniq.Name{iii})&temp.Logical_ii_Time_ch;
        
            temp.Name=DataTable.Name(temp.Logical_iii_Time_ch);
            temp.Power=DataTable.Power(temp.Logical_iii_Time_ch);
            
            temp.Power_arr=cell2mat(temp.Power');
            temp.Power_Mean=mean(temp.Power_arr,2); % save this variable
            temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
            
            if sum(temp.Logical_iii_Time_ch)~=0   
                
            % save every needed variable in a temporary 'table'
            count=count+1;
            temptable.Method(count)=uniq.Method(i);
            temptable.Name(count)=uniq.Name(iii);
            temptable.Power{count}=temp.Power_Mean;
            temptable.Genotype(count)=DataTable.Genotype(find(temp.Logical_iii_Time_ch,1,'first'));
            temptable.Targeted(count)=DataTable.Targeted(find(temp.Logical_iii_Time_ch,1,'first'));
            temptable.Channel(count)=DataTable.Channel(find(temp.Logical_iii_Time_ch,1,'first'));  
        
            % for the consistency of the color and style of lines
            c_index=find(strcmp(uniq.Name_Full,temp.Name{1}));
            while c_index>size(c,1)
                c_index=c_index-size(c,1);
            end
            ls_index=find(strcmp(uniq.Name_Full,temp.Name{1}));
            while ls_index>length(ls)
                ls_index=ls_index-length(ls);
            end
               
            for jj=1:2
                
            if jj==2 % Normalized plot
                temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
                temp.Power_Mean=mean(temp.Power_arr,2);
                temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
            end
                
            subplot(1,2,jj), hold on
            hh(iii,jj)=plot(fff,temp.Power_Mean,...
                'Color',c(c_index,:),...
                'LineStyle',ls{ls_index},...
                'LineWidth',lw);
            patch([fff;...
                   fff(end:-1:1);...
                   fff(1)],...
                  [temp.Power_Mean-temp.SEM;...
                   temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
                   temp.Power_Mean(1)-temp.SEM(1)],...
                   c(c_index,:),...
                   'FaceAlpha',SEM_transparency,...
                   'EdgeAlpha',SEM_transparency)
            end
            
            end
                
        end
            
        subplot(1,2,1)
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
        title([uniq.Targeted{ii},' ',uniq.Method{i},' -Absolute'])
        legend(hh(:,1),{uniq.Name{:}}','Interpreter','none')
xticks('auto')
% set(gca,'Yscale','log')
grid on
        subplot(1,2,2)
        xlim([Settings.HP,maxFreq])
        xticks(linspace(0,maxFreq,6))
        xlabel('Frequency (Hz)'),ylabel('LFPs PSD(%)')
        title([uniq.Targeted{ii},' ',uniq.Method{i},' -Normalized'])
        legend(hh(:,2),{uniq.Name{:}}','Interpreter','none')
xticks('auto')
% set(gca,'Yscale','log')
grid on
            
    end
    end
end

%% Targeted mean in figure for each method, genotypes plotted together, with channels separated
SEM_transparency=0.2;
maxFreq=20;

% Create 1 figure for each method and a cell within a cell for each subplot
uniq.Method_Full=unique(temptable.Method);
nMethod_Full=length(uniq.Method_Full);
hh_targ=cell(nMethod_Full,2);
legend_targ=cell(nMethod_Full,2);
bandpower_cell=cell(nMethod_Full,1);
for i=1:nMethod_Full
fig.method(i)=figure('Name',['Mean Genotype-Targeted for ',uniq.Method_Full{i},' with channels separated']);
bandpower_cell{i}=cell(2,1);
end

index=0;
count=0;
c_index=index;
ls_index=index;

uniq.Genotype=unique(temptable.Genotype);
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv"];
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv","STOPrescue-129/Sv"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:length(uniq.Genotype)
temp.Logical_g=strcmp(temptable.Genotype,uniq.Genotype{g});
uniq.Method=unique(temptable.Method(temp.Logical_g));
for i=1:length(uniq.Method)
    
% make appropiate method figure current
ind.method=find(strcmp(uniq.Method_Full,uniq.Method(i)));
set(0,'CurrentFigure',fig.method(ind.method))

temp.Logical_g_i=temp.Logical_g&strcmp(temptable.Method,uniq.Method{i});
uniq.Targeted=unique(temptable.Targeted(temp.Logical_g_i));
for ii=1:length(uniq.Targeted)
temp.Logical_g_ii=temp.Logical_g_i&strcmp(temptable.Targeted,uniq.Targeted{ii});
uniq.Channel=unique(temptable.Channel(temp.Logical_g_ii));
for iii=1:length(uniq.Channel)
temp.Logical_g_iii=temp.Logical_g_ii&strcmp(temptable.Channel,uniq.Channel{iii});

temp.Power=temptable.Power(temp.Logical_g_iii);

temp.Power_arr=cell2mat(temp.Power);
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));

if sum(temp.Logical_g_iii)~=0   
    
    % no consistency because of 3 variables
    c_index=ii;
    while c_index>size(c,1)
    c_index=c_index-size(c,1);
    end
    ls_index=g;
    while ls_index>length(ls)
    ls_index=ls_index-length(ls);
    end
    
for jj=1:2
        
         
    if jj==2 % Normalized plot
        temp.Power_arr_norm=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr_norm,2);
        temp.SEM=std(temp.Power_arr_norm,0,2)/sqrt(size(temp.Power_arr_norm,2));
    end
         
    subplot(1,2,jj), hold on
    hh_targ_temp=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    hh_targ{ind.method,jj}=[hh_targ{ind.method,jj};hh_targ_temp];
    legend_targ_temp=[uniq.Genotype{g},' ',uniq.Targeted{ii},' ',uniq.Method{i}];
    legend_targ{ind.method,jj}=[legend_targ{ind.method,jj};{legend_targ_temp}];
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
 
end   
                     
% temporary 'table'
count=count+1;
bandtable3.Genotype(count)=uniq.Genotype(g);
bandtable3.Method(count)=uniq.Method(i);
bandtable3.Targeted(count)=uniq.Targeted(ii);
bandtable3.Channel(count)=uniq.Channel(iii);
bandtable3.Power{count}=[bandpower2(temp.Power_arr,fff,delta);...
                        bandpower2(temp.Power_arr,fff,theta);...
                        bandpower2(temp.Power_arr,fff,beta);...
                        bandpower2(temp.Power_arr,fff,gamma);...
                        bandpower2(temp.Power_arr,fff)];
         
end
end
end
end
end

for i=1:nMethod_Full
if exist('hh_targ')~=0
    
set(0,'CurrentFigure',fig.method(i))

%absolute settings
subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Genotype-Targeted-Method -Absolute')
legend(hh_targ{i,1},legend_targ{i,1},'Interpreter','none')

%normalized settings
subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Genotype-Targeted-Method -Normalized')
legend(hh_targ{i,2},legend_targ{i,2},'Interpreter','none')

end
end

%% Histogram Targeted without channels separated
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% clear h_bin variable possibly used in previous sections
clearvars h_bin

% figure for each method-Genotype
% subplots with each genotype-band for both absolute-normalized
uniq.Genotype_Full=unique(bandtable3.Genotype);
nGenotype_Full=length(uniq.Genotype_Full);
uniq.Channel_Full=unique(bandtable3.Channel);
nChannel_Full=length(uniq.Channel_Full);
for ii=1:nGenotype_Full
    for iii=1:nChannel_Full
        fig.Genotype_Channel(ii,iii)=figure('Name',['Targeted band powers for ',uniq.Genotype_Full{ii},' ',uniq.Channel_Full{iii}]);
    end
end

Targeted_bandtable=unique(bandtable3.Targeted);
nTargeted_bandtable=length(Targeted_bandtable);

uniq.Method=unique(bandtable3.Method);
nMethod=length(uniq.Method);
power_cellmatrix=cell(nMethod,nGenotype_Full,nChannel_Full,nTargeted_bandtable,nTargeted_bandtable,nBand);
for i=1:nMethod
logical.i=strcmp(bandtable3.Method,uniq.Method(i));
uniq.Genotype=unique(bandtable3.Genotype(logical.i));
nGenotype=length(uniq.Genotype);
for ii=1:nGenotype
logical.ii=logical.i&strcmp(bandtable3.Genotype,uniq.Genotype(ii));
uniq.Channel=unique(bandtable3.Channel);
for iii=1:nChannel
logical.iii=logical.ii&strcmp(bandtable3.Channel,uniq.Channel(iii));
set(0,'CurrentFigure',fig.Genotype_Channel(ii,iii))

uniq.Targeted=unique(bandtable3.Targeted(logical.iii));
nTargeted=length(uniq.Targeted);

%%% define bin location
temp.step2=bbd+bw;
temp.step1=bw+(nTargeted-1)*(temp.step2)+bgd;
select_points=(bgd:temp.step1:bgd+(nBand-1)*temp.step1)';
bin_matrix=zeros(nBand,nTargeted);
for j=1:nBand
    bin_matrix(j,:)=select_points(j)+bw/2:temp.step2:select_points(j)+bw/2+(nTargeted-1)*temp.step2;
end

for t=1:nTargeted
logical.iii_t=logical.iii&strcmp(bandtable3.Targeted,uniq.Targeted(t));
for jj=1:2
    
%%% start plot
subplot(nMethod,2,(i-1)*nMethod+jj), hold on
for j=1:nBand
    if jj==1
        temp.power=bandtable3.Power{find(logical.iii_t,1,'first')}(j,:);
    else %normalize
        temp.power=(bandtable3.Power{find(logical.iii_t,1,'first')}(j,:)./bandtable3.Power{find(logical.iii_t,1,'first')}(end,:))*100;
    end
    ind.ii=find(strcmp(uniq.Genotype_Full,uniq.Genotype(ii)));
    ind.iii=find(strcmp(uniq.Channel_Full,uniq.Channel(iii)));
    ind.g=find(strcmp(uniq.Targeted(t),Targeted_bandtable));
    power_cellmatrix{i,ind.ii,ind.iii,ind.g,j,jj}=temp.power;
    temp.power_mean=mean(temp.power);
    temp.SEM=std(temp.power)/sqrt(length(temp.power));
    bar(bin_matrix(j,t),temp.power_mean,...
    'FaceColor',c(t,:),'BarWidth',bw);
    errorbar(bin_matrix(j,t),temp.power_mean,temp.SEM,...
    'Color',c(t,:),'LineStyle','none','LineWidth',lw)
end

%xticks
xticks(select_points+(bw+(nTargeted-1)*temp.step2)/2)
xticklabels(bandstr)

%legend
h_bin=gobjects(nTargeted,1);
for tt=1:nTargeted
    h_bin(tt)=bar(NaN,'FaceColor',c(tt,:),'BarWidth',bw);
end
legend(h_bin,cellstr(uniq.Targeted),'Interpreter','none')

%title and ylabel
if jj==1
title([uniq.Genotype{ii},' ',uniq.Method{i},' -absolute'])
ylabel('Power (\muV^2)')
else
title([uniq.Genotype{ii},' ',uniq.Method{i},' -normalized'])
ylabel('Power (%)')
end
%%%% end plot

end
end
end
end
end

%Save statistics to table
count=0;
for i=1:nMethod 
for ii=1:nGenotype_Full
for iii=1:nChannel_Full
for j=1:nBand
for jj=1:2
    
    if jj==1
        normstr="Absolute";
    else
        normstr="Normalized";
    end

H=cell(nTargeted_bandtable); 
P=H;
CI=H;
STATS=H;
for t=1:nTargeted_bandtable
for tt=1:nTargeted_bandtable
    
    if ~isempty(power_cellmatrix{i,ii,iii,t,j,jj})&&~isempty(power_cellmatrix{i,ii,iii,tt,j,jj})
    [H{t,tt},P{t,tt},CI{t,tt},STATS{t,tt}]...
    =ttest2(power_cellmatrix{i,ii,iii,t,j,jj},power_cellmatrix{i,ii,iii,tt,j,jj},...
    'vartype','unequal');
    end
    
end
end

    count=count+1;
    
    bandtable2_Targeted.method(count,1)=uniq.Method(i);
    bandtable2_Targeted.Genotype(count,1)=uniq.Genotype_Full(ii);
    bandtable2_Targeted.Channel(count,1)=uniq.Channel_Full(iii);
    bandtable2_Targeted.normalization(count,1)=normstr;
    bandtable2_Targeted.band(count,1)=bandstr(j);
    
    temp.charH='H(:,1)';
    temp.charP='P(:,1)';
    temp.charCI='CI(:,1)';
    temp.charSTATS='STATS(:,1)';
    for ttt=2:nTargeted_bandtable
        temp.charH=[temp.charH,',H(:,',num2str(ttt),')'];
        temp.charP=[temp.charP,',P(:,',num2str(ttt),')'];
        temp.charCI=[temp.charCI,',CI(:,',num2str(ttt),')'];
        temp.charSTATS=[temp.charSTATS,',STATS(:,',num2str(ttt),')'];
    end
    bandtable2_Targeted.H{count,1}=eval(['table(',temp.charH,')']);
    bandtable2_Targeted.H{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.H{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P{count,1}=eval(['table(',temp.charP,')']);
    bandtable2_Targeted.P{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.CI{count,1}=eval(['table(',temp.charCI,')']);
    bandtable2_Targeted.CI{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.CI{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.STATS{count,1}=eval(['table(',temp.charSTATS,')']);
    bandtable2_Targeted.STATS{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.STATS{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P_matrix{count,1}=cell2mat(P);
    
end
end
end
end
end
TableNames={'Method','Genotype','Channel','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix'};
BandTable_Targeted_WithChannelSeparation=table(bandtable2_Targeted.method,bandtable2_Targeted.Genotype,bandtable2_Targeted.Channel,bandtable2_Targeted.normalization,bandtable2_Targeted.band,...
                bandtable2_Targeted.H,bandtable2_Targeted.P,bandtable2_Targeted.CI,bandtable2_Targeted.STATS,bandtable2_Targeted.P_matrix,...
                'VariableNames',TableNames);
save('BandTable_Targeted.mat', 'BandTable_Targeted_WithChannelSeparation');

%% Targeted mean in figure for each method, genotypes plotted together, without channels separated
clearvars bandtable3 bandtable2_Targeted
SEM_transparency=0.2;
maxFreq=20;

% Create 1 figure for each method and a cell within a cell for each subplot
uniq.Method_Full=unique(temptable.Method);
nMethod_Full=length(uniq.Method_Full);
hh_targ=cell(nMethod_Full,2);
legend_targ=cell(nMethod_Full,2);
bandpower_cell=cell(nMethod_Full,1);
for i=1:nMethod_Full
fig.method(i)=figure('Name',['Mean Genotype-Targeted for ',uniq.Method_Full{i}]);
bandpower_cell{i}=cell(2,1);
end

index=0;
count=0;
c_index=index;
ls_index=index;

uniq.Genotype=unique(temptable.Genotype);
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv"];
% uniq.Genotype=["WT-129/Sv","STOP-129/Sv","STOPrescue-129/Sv"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:length(uniq.Genotype)
temp.Logical_g=strcmp(temptable.Genotype,uniq.Genotype{g});
uniq.Method=unique(temptable.Method(temp.Logical_g));
for i=1:length(uniq.Method)
    
% make appropiate method figure current
ind.method=find(strcmp(uniq.Method_Full,uniq.Method(i)));
set(0,'CurrentFigure',fig.method(ind.method))

temp.Logical_g_i=temp.Logical_g&strcmp(temptable.Method,uniq.Method{i});
uniq.Targeted=unique(temptable.Targeted(temp.Logical_g_i));
for ii=1:length(uniq.Targeted)
temp.Logical_g_ii=temp.Logical_g_i&strcmp(temptable.Targeted,uniq.Targeted{ii});

temp.Power=temptable.Power(temp.Logical_g_ii);

temp.Power_arr=cell2mat(temp.Power);
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));

if sum(temp.Logical_g_ii)~=0   
    
    % no consistency because of 3 variables
    c_index=ii;
    while c_index>size(c,1)
    c_index=c_index-size(c,1);
    end
    ls_index=g;
    while ls_index>length(ls)
    ls_index=ls_index-length(ls);
    end
    
for jj=1:2
        
         
    if jj==2 % Normalized plot
        temp.Power_arr_norm=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr_norm,2);
        temp.SEM=std(temp.Power_arr_norm,0,2)/sqrt(size(temp.Power_arr_norm,2));
    end
         
    subplot(1,2,jj), hold on
    hh_targ_temp=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    hh_targ{ind.method,jj}=[hh_targ{ind.method,jj};hh_targ_temp];
    legend_targ_temp=[uniq.Genotype{g},' ',uniq.Targeted{ii},' ',uniq.Method{i}];
    legend_targ{ind.method,jj}=[legend_targ{ind.method,jj};{legend_targ_temp}];
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
 
end   
                     
% temporary 'table'
count=count+1;
bandtable3.Genotype(count)=uniq.Genotype(g);
bandtable3.Method(count)=uniq.Method(i);
bandtable3.Targeted(count)=uniq.Targeted(ii);
bandtable3.Power{count}=[bandpower2(temp.Power_arr,fff,delta);...
                        bandpower2(temp.Power_arr,fff,theta);...
                        bandpower2(temp.Power_arr,fff,beta);...
                        bandpower2(temp.Power_arr,fff,gamma);...
                        bandpower2(temp.Power_arr,fff)];
         
end
end

end
end

for i=1:nMethod_Full
if exist('hh_targ')~=0
    
set(0,'CurrentFigure',fig.method(i))

%absolute settings
subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Genotype-Targeted-Method -Absolute')
legend(hh_targ{i,1},legend_targ{i,1},'Interpreter','none')

%normalized settings
subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Genotype-Targeted-Method -Normalized')
legend(hh_targ{i,2},legend_targ{i,2},'Interpreter','none')

end
end

%% Histogram Targeted without channels separated
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% clear h_bin variable possibly used in previous sections
clearvars h_bin

% figure for each method-Genotype
% subplots with each genotype-band for both absolute-normalized
uniq.Genotype_Full=unique(bandtable3.Genotype);
nGenotype_Full=length(uniq.Genotype);
for ii=1:nGenotype_Full
    fig.Genotype(ii)=figure('Name',['Targeted band powers for ',uniq.Genotype_Full{ii}]);
end

Targeted_bandtable=unique(bandtable3.Targeted);
nTargeted_bandtable=length(Targeted_bandtable);

uniq.Method=unique(bandtable3.Method);
nMethod=length(uniq.Method);
power_cellmatrix=cell(nMethod,nGenotype_Full,nTargeted_bandtable,nTargeted_bandtable,nBand);
for i=1:nMethod
logical.i=strcmp(bandtable3.Method,uniq.Method(i));
uniq.Genotype=unique(bandtable3.Genotype(logical.i));
nGenotype=length(uniq.Genotype);
for ii=1:nGenotype
set(0,'CurrentFigure',fig.Genotype(ii))
logical.ii=logical.i&strcmp(bandtable3.Genotype,uniq.Genotype(ii));
uniq.Targeted=unique(bandtable3.Targeted(logical.ii));
nTargeted=length(uniq.Targeted);

%%% define bin location
temp.step2=bbd+bw;
temp.step1=bw+(nTargeted-1)*(temp.step2)+bgd;
select_points=(bgd:temp.step1:bgd+(nBand-1)*temp.step1)';
bin_matrix=zeros(nBand,nTargeted);
for j=1:nBand
    bin_matrix(j,:)=select_points(j)+bw/2:temp.step2:select_points(j)+bw/2+(nTargeted-1)*temp.step2;
end

for t=1:nTargeted
logical.ii_t=logical.ii&strcmp(bandtable3.Targeted,uniq.Targeted(t));
for jj=1:2
    
%%% start plot
subplot(nMethod,2,(i-1)*nMethod+jj), hold on
for j=1:nBand
    if jj==1
        temp.power=bandtable3.Power{find(logical.ii_t,1,'first')}(j,:);
    else %normalize
        temp.power=(bandtable3.Power{find(logical.ii_t,1,'first')}(j,:)./bandtable3.Power{find(logical.ii_t,1,'first')}(end,:))*100;
    end
    ind.ii=find(strcmp(uniq.Genotype_Full,uniq.Genotype(ii)));
    ind.g=find(strcmp(uniq.Targeted(t),Targeted_bandtable));
    power_cellmatrix{i,ind.ii,ind.g,j,jj}=temp.power;
    temp.power_mean=mean(temp.power);
    temp.SEM=std(temp.power)/sqrt(length(temp.power));
    bar(bin_matrix(j,t),temp.power_mean,...
    'FaceColor',c(t,:),'BarWidth',bw);
    errorbar(bin_matrix(j,t),temp.power_mean,temp.SEM,...
    'Color',c(t,:),'LineStyle','none','LineWidth',lw)
end

%xticks
xticks(select_points+(bw+(nTargeted-1)*temp.step2)/2)
xticklabels(bandstr)

%legend
for tt=1:nTargeted
    h_bin(tt)=bar(NaN,'FaceColor',c(tt,:),'BarWidth',bw);
end
legend(h_bin,cellstr(uniq.Targeted),'Interpreter','none')

%title and ylabel
if jj==1
title([uniq.Genotype{ii},' ',uniq.Method{i},' -absolute'])
ylabel('Power (\muV^2)')
else
title([uniq.Genotype{ii},' ',uniq.Method{i},' -normalized'])
ylabel('Power (%)')
end
%%%% end plot

end
end
end
end

%Save statistics to table
count=0;
for i=1:nMethod 
for ii=1:nGenotype_Full
for j=1:nBand
for jj=1:2
    
    if jj==1
        normstr="Absolute";
    else
        normstr="Normalized";
    end

H=cell(nTargeted_bandtable); 
P=H;
CI=H;
STATS=H;
for t=1:nTargeted_bandtable
for tt=1:nTargeted_bandtable
    
    if ~isempty(power_cellmatrix{i,ii,t,j,jj})&&~isempty(power_cellmatrix{i,ii,tt,j,jj})
    [H{t,tt},P{t,tt},CI{t,tt},STATS{t,tt}]...
    =ttest2(power_cellmatrix{i,ii,t,j,jj},power_cellmatrix{i,ii,tt,j,jj},...
    'vartype','unequal');
    end
    
end
end

    count=count+1;
    
    bandtable2_Targeted.method(count,1)=uniq.Method(i);
    bandtable2_Targeted.Genotype(count,1)=uniq.Genotype_Full(ii);
    bandtable2_Targeted.normalization(count,1)=normstr;
    bandtable2_Targeted.band(count,1)=bandstr(j);
    
    temp.charH='H(:,1)';
    temp.charP='P(:,1)';
    temp.charCI='CI(:,1)';
    temp.charSTATS='STATS(:,1)';
    for ttt=2:nTargeted_bandtable
        temp.charH=[temp.charH,',H(:,',num2str(ttt),')'];
        temp.charP=[temp.charP,',P(:,',num2str(ttt),')'];
        temp.charCI=[temp.charCI,',CI(:,',num2str(ttt),')'];
        temp.charSTATS=[temp.charSTATS,',STATS(:,',num2str(ttt),')'];
    end
    bandtable2_Targeted.H{count,1}=eval(['table(',temp.charH,')']);
    bandtable2_Targeted.H{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.H{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P{count,1}=eval(['table(',temp.charP,')']);
    bandtable2_Targeted.P{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.CI{count,1}=eval(['table(',temp.charCI,')']);
    bandtable2_Targeted.CI{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.CI{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.STATS{count,1}=eval(['table(',temp.charSTATS,')']);
    bandtable2_Targeted.STATS{count,1}.Properties.RowNames={Targeted_bandtable{:}};
    bandtable2_Targeted.STATS{count,1}.Properties.VariableNames={Targeted_bandtable{:}};
    bandtable2_Targeted.P_matrix{count,1}=cell2mat(P);
    
end
end
end
end
TableNames={'Method','Genotype','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix'};
BandTable_Targeted_WithoutChannelSeparation=table(bandtable2_Targeted.method,bandtable2_Targeted.Genotype,bandtable2_Targeted.normalization,bandtable2_Targeted.band,...
                bandtable2_Targeted.H,bandtable2_Targeted.P,bandtable2_Targeted.CI,bandtable2_Targeted.STATS,bandtable2_Targeted.P_matrix,...
                'VariableNames',TableNames);
save('BandTable_Targeted.mat', 'BandTable_Targeted_WithoutChannelSeparation');

end % end of if plot_targeted statement

if plot_extramean
%% Mean Channel-Method over Names and over Days
    
SEM_transparency=0.2;
maxFreq=20;

figure('Name','Mean Channel-Method over Names and over Days')
uniq.Channel_Full=unique(DataTable.Channel);
uniq.Method=unique(DataTable.Method);
count=0;
legendCell={};
clearvars hh

for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
uniq.Channel=unique(DataTable.Channel(temp.Logical_i));
    for ii=1:length(uniq.Channel)
    temp.Logical_ii=strcmp(DataTable.Channel,uniq.Channel{ii})&temp.Logical_i;
    temp.Logical_ii_Time=temp.Logical_ii&temp.Logical_Time;
    
    temp.Power=DataTable.Power(temp.Logical_ii_Time);

    temp.Power_arr=cell2mat(temp.Power');
    temp.Power_Mean=mean(temp.Power_arr,2);
    temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
    if sum(temp.Logical_ii_Time)~=0   
    legendCell=[legendCell,[uniq.Channel{ii},' ',uniq.Method{i}]];
    count=count+1;
    for jj=1:2
        
        % for the consistency of the color and style of lines
        c_index=find(strcmp(uniq.Channel_Full,uniq.Channel{ii})); % Color dependent on Channel
        while c_index>size(c,1)
        c_index=c_index-size(c,1);
        end
        ls_index=find(strcmp(uniq.Method,uniq.Method{i})); % Line Style dependent on Method
        while ls_index>length(ls)
        ls_index=ls_index-length(ls);
        end
         
        if jj==2 % Normalized plot
            temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
            temp.Power_Mean=mean(temp.Power_arr,2);
            temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
        end
         
        subplot(1,2,jj), hold on
        hh(count,jj)=plot(fff,temp.Power_Mean,...
            'Color',c(c_index,:),...
            'LineStyle',ls{ls_index},...
            'LineWidth',lw);
        patch([fff;...
               fff(end:-1:1);...
               fff(1)],...
              [temp.Power_Mean-temp.SEM;...
               temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
               temp.Power_Mean(1)-temp.SEM(1)],...
               c(c_index,:),...
               'FaceAlpha',SEM_transparency,...
               'EdgeAlpha',SEM_transparency)
           
    end
         
    end
    end
end

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Channel-Method -Absolute')
legend(hh(:,1),legendCell,'Interpreter','none')

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Channel-Method -Normalized')
legend(hh(:,2),legendCell,'Interpreter','none')

%% Mean Channel over Methods, Names and Days
SEM_transparency=0.2;
maxFreq=20;

figure('Name','Mean Channel over Methods, Names and Days')
uniq.Channel=unique(DataTable.Channel);
hh=zeros(length(uniq.Channel),2);

for i=1:length(uniq.Channel)
temp.Logical_i=strcmp(DataTable.Channel,uniq.Channel{i});
temp.Logical_i_Time=temp.Logical_i&temp.Logical_Time;

temp.Power=DataTable.Power(temp.Logical_i_Time);

temp.Power_arr=cell2mat(temp.Power');
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
if sum(temp.Logical_i_Time)~=0   
for jj=1:2
        
    % for the consistency of the color and style of lines
    c_index=find(strcmp(uniq.Channel,uniq.Channel{i})); % Color dependent on Channel
    while c_index>size(c,1)
    c_index=c_index-size(c,1);
    end
    ls_index=find(strcmp(uniq.Channel,uniq.Channel{i})); % Line Style dependent on Method
    while ls_index>length(ls)
    ls_index=ls_index-length(ls);
    end
         
    if jj==2 % Normalized plot
        temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr,2);
        temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    end
         
    subplot(1,2,jj), hold on
    hh(i,jj)=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
           
end
         
end
end

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Channel -Absolute')
legend(hh(:,1),{uniq.Channel{:}},'Interpreter','none')

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Channel -Normalized')
legend(hh(:,2),{uniq.Channel{:}},'Interpreter','none')

%% Mean Method over Channels, Names and Days
SEM_transparency=0.2;
maxFreq=20;

figure('Name','Mean Method over Channels, Names and Days')
uniq.Method=unique(DataTable.Method);
hh=zeros(length(uniq.Method),2);

for i=1:length(uniq.Method)
temp.Logical_i=strcmp(DataTable.Method,uniq.Method{i});
temp.Logical_i_Time=temp.Logical_i&temp.Logical_Time;

temp.Power=DataTable.Power(temp.Logical_i_Time);

temp.Power_arr=cell2mat(temp.Power');
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
if sum(temp.Logical_i_Time)~=0   
for jj=1:2
        
    % for the consistency of the color and style of lines
    c_index=find(strcmp(uniq.Method,uniq.Method{i})); % Color dependent on Channel
    while c_index>size(c,1)
    c_index=c_index-size(c,1);
    end
    ls_index=find(strcmp(uniq.Method,uniq.Method{i})); % Line Style dependent on Method
    while ls_index>length(ls)
    ls_index=ls_index-length(ls);
    end
         
    if jj==2 % Normalized plot
        temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr,2);
        temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    end
         
    subplot(1,2,jj), hold on
    hh(i,jj)=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
           
end
         
end
end

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Method -Absolute')
legend(hh(:,1),{uniq.Method{:}},'Interpreter','none')

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Method -Normalized')
legend(hh(:,2),{uniq.Method{:}},'Interpreter','none')

%% Mean Genotype-Channel over Methods, Names and Days (norm and abs)
SEM_transparency=0.2;
maxFreq=20;

figure('Name','Mean Genotype-Channel')
legendCell={};
uniq.Channel_Full=unique(DataTable.Channel);
index=0;
clearvars hh

uniq.Genotype=unique(DataTable.Genotype);
for g=1:length(uniq.Genotype)
temp.Logical_g=strcmp(DataTable.Genotype,uniq.Genotype{g});
uniq.Method=unique(DataTable.Method); 
for i=1:length(uniq.Channel)
temp.Logical_g_i=temp.Logical_g&strcmp(DataTable.Channel,uniq.Channel{i});
temp.Logical_g_i_Time=temp.Logical_g_i&temp.Logical_Time;

temp.Power=DataTable.Power(temp.Logical_g_i_Time);

temp.Power_arr=cell2mat(temp.Power');
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
if sum(temp.Logical_g_i_Time)~=0   
legendCell=[legendCell,[uniq.Genotype{g},' ',uniq.Channel{i}]];
index=index+1;

% no consistency because of 3 variables
c_index=find(strcmp(uniq.Genotype,uniq.Genotype{g}));
while c_index>size(c,1)
c_index=c_index-size(c,1);
end
ls_index=find(strcmp(uniq.Channel_Full,uniq.Channel{i}));
while ls_index>length(ls)
ls_index=ls_index-length(ls);
end
    
for jj=1:2
       
         
    if jj==2 % Normalized plot
        temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr,2);
        temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    end
         
    subplot(1,2,jj), hold on
    hh(index,jj)=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
           
end
         
end
end

end

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Genotype-Channel -Absolute')
legend(hh(:,1),legendCell,'Interpreter','none')

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Genotype-Channel -Normalized')
legend(hh(:,2),legendCell,'Interpreter','none')

%% Mean Genotype over Channels, Methods, Names and Days (norm and abs)
clearvars hh
SEM_transparency=0.2;
maxFreq=20;

figure('Name','Mean Genotype')
uniq.Channel_Full=unique(DataTable.Channel);
index=0;

legendCell={};
uniq.Genotype=unique(DataTable.Genotype);
nGenotype=length(uniq.Genotype);
for g=1:nGenotype
temp.Logical_g=strcmp(DataTable.Genotype,uniq.Genotype{g});
temp.Logical_g_Time=temp.Logical_g&temp.Logical_Time;

temp.Power=DataTable.Power(temp.Logical_g_Time);

temp.Power_arr=cell2mat(temp.Power');
temp.Power_Mean=mean(temp.Power_arr,2);
temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    
if sum(temp.Logical_g_Time)~=0   
index=index+1;
legendCell=[legendCell,uniq.Genotype{g}];

% no consistency because of 3 variables
c_index=find(strcmp(uniq.Genotype,uniq.Genotype{g}));
while c_index>size(c,1)
c_index=c_index-size(c,1);
end
ls_index=find(strcmp(uniq.Genotype,uniq.Genotype{g}));
while ls_index>length(ls)
ls_index=ls_index-length(ls);
end
    
for jj=1:2
       
         
    if jj==2 % Normalized plot
        temp.Power_arr=(temp.Power_arr./bandpower2(temp.Power_arr,fff,[Settings.HP,normFreq]))*100;
        temp.Power_Mean=mean(temp.Power_arr,2);
        temp.SEM=std(temp.Power_arr,0,2)/sqrt(size(temp.Power_arr,2));
    end
         
    subplot(1,2,jj), hold on
    hh(index,jj)=plot(fff,temp.Power_Mean,...
        'Color',c(c_index,:),...
        'LineStyle',ls{ls_index},...
        'LineWidth',lw);
    patch([fff;...
           fff(end:-1:1);...
           fff(1)],...
          [temp.Power_Mean-temp.SEM;...
           temp.Power_Mean(end:-1:1)+temp.SEM(end:-1:1);...
           temp.Power_Mean(1)-temp.SEM(1)],...
           c(c_index,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
           
end
         
end
end

subplot(1,2,1)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (\muV^2)')
title('Mean Genotype -Absolute')
legend(hh(:,1),cellstr(legendCell),'Interpreter','none')

subplot(1,2,2)
xlim([Settings.HP,maxFreq])
xlabel('Frequency (Hz)'),ylabel('LFPs PSD (%)')
title('Mean Genotype -Normalized')
legend(hh(:,2),cellstr(legendCell),'Interpreter','none')

end % end of if statement extra mean