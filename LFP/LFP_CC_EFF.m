% Select Experiment Folder Folder (EFF).
% Calculates cross-correlation for all channels in Experiment Folders
% (EFs) in EFF.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Initial input

% Select Experiment Folder Folder and create list of Experiment Folder names
EFF = uigetdir('','Select folder containing experiments');
cd (EFF)
EF=dir('*_*_*');
nEF=length(EF);
EF_Name={EF(:).name}';

% Load Channel names and settings of analysis
load('Channels')
nChannel=size(Channels,1);
load('Settings')

% Use settings to define cross correlation settings
nMaxLag=.1*Settings.fs; % maximum lag of .1 s

% Plot settings
col=colormap(lines);
close(gcf)

% Get names and days
Name=cell(nEF,1);
Day=Name;
for i=1:nEF
    tempsplit=strsplit(EF_Name{i},'_');
    Name{i}=tempsplit{1};
    Day{i}=tempsplit{end};
end
uniq.Name=unique(Name);
nName=length(uniq.Name);
uniq.Day_Full=unique(Day);
nDay_Full=length(uniq.Day_Full);

% preallocate
c_EF=cell(nEF,1);
t_EF=c_EF;
l_EF=c_EF;
CC_dis_EF=c_EF;
lag_dis_EF=c_EF;

% predefine figures
check_matrix=triu(true(nChannel),1);
for ii=1:nChannel
    for iii=1:nChannel
        if check_matrix(ii,iii)
            fig.cross1{ii,iii}=figure('Name',[Channels{ii,1},'-',Channels{iii,1},' cross-correlation']);
        end
    end
end

for i=1:nEF
    [c_EF{i},t_EF{i},l_EF{i},CC_dis_EF{i},lag_dis_EF{i}]=LFP_CC_EF(EF_Name{i},Settings,Channels,nMaxLag);
    for ii=1:nChannel
        for iii=1:nChannel
            if check_matrix(ii,iii)
                set(0,'CurrentFigure',fig.cross1{ii,iii})
                ind.Name=find(strcmp(uniq.Name,Name{i}));
                ind.Day=find(strcmp(uniq.Day_Full,Day{i}));
                subplot(nName,nDay_Full,ind.Day+(ind.Name-1)*nDay_Full)
                imagesc(t_EF{i}{ii,iii},l_EF{i}{ii,iii},c_EF{i}{ii,iii})
                grid on
                set(gca,'YDir','normal')
                title([Name{i},'_',Day{i}],'Interpreter','none')
            end
        end
    end
end

CC_mean=cell(nChannel);
lag_mean=CC_mean;
for ROW=1:nChannel
for COL=1:nChannel    
if check_matrix(ROW,COL)
    
for i=1:nName
    fig.CC_lag{i}=figure('Name',[uniq.Name{i},'_CC/lag fig:name holdon:day']);
    logical.Name=strcmp(uniq.Name{i},Name);
    uniq.Day=unique(Day(logical.Name));
    nDay=length(uniq.Day);
    clearvars CCtemp lagtemp
    CCtemp=zeros(nDay,1); lagtemp=zeros(nDay,1);
    for ii=1:nDay
        logical.Name_Day=logical.Name&strcmp(uniq.Day{ii},Day);
        ind.Day=find(strcmp(uniq.Day{ii},uniq.Day_Full));
        set(0,'CurrentFigure',fig.CC_lag{i})
        subplot(1,2,1)
        hold on
        histogram(CC_dis_EF{logical.Name_Day}{ROW,COL},...
            'FaceColor',col(ind.Day,:),'BinWidth',.05)
        CCtemp(ii)=mean(CC_dis_EF{logical.Name_Day}{ROW,COL});
        subplot(1,2,2)
        hold on
        histogram(lag_dis_EF{logical.Name_Day}{ROW,COL},...
            'FaceColor',col(ind.Day,:),'BinWidth',1/Settings.fs)
        lagtemp(ii)=mean(lag_dis_EF{logical.Name_Day}{ROW,COL});
    end
    subplot(1,2,1)
    xlabel('CrossCorrelation coefficient (CC)')
    ylabel('Number of windows')
    title(['Crosscorrelation distribution, mean CC is ',num2str(mean(CCtemp))])
    legend(uniq.Day)
    subplot(1,2,2)
    legend(uniq.Day)
    xlabel('Lag (s)')
    ylabel('Number of windows')
    title(['Lag distribution, mean lag is ',num2str(mean(lagtemp)*1000),' ms'])
end

% Calculate mean for CC and lag distributions
CC_mean{ROW,COL}=zeros(nEF,1);
lag_mean{ROW,COL}=CC_mean{ROW,COL};
for i=1:nEF
CC_mean{ROW,COL}(i)=mean(CC_dis_EF{i}{ROW,COL});
lag_mean{ROW,COL}(i)=mean(lag_dis_EF{i}{ROW,COL});
end

end
end
end

% CC_mean_name=zeros(nName,1);
% lag_mean_name=CC_mean_name;
% for i=1:nName
%     logical.Name=strcmp(uniq.Name{i},Name);
%     lag_mean_name(i)=mean(lag_mean(logical.Name));
%     CC_mean_name(i)=mean(CC_mean(logical.Name));
% end

% % Genotypes
% Genotype=cell(nName,1);
% Genotype(:)={'WT_C57'};
% % Genotype(end-4:end-1)={'het_C57'}; % G8
% Genotype([5;7;8;9])={'het_C57'}; % G7
% 
% figure('Name','CC/lag WT/het')
% uniq.Genotype=unique(Genotype);
% nGenotype=length(uniq.Genotype);
% for i=1:nGenotype
%     logical.Genotype=strcmp(uniq.Genotype{i},Genotype);
%     nGenotype_spec=sum(logical.Genotype);
%     subplot(1,2,1)
%     hold on
%     grid on
%     histogram(CC_mean_name(logical.Genotype),'BinWidth',.1,'FaceColor',col(i,:))
%     subplot(1,2,2)
%     hold on
%     grid on
%     histogram(lag_mean_name(logical.Genotype),'BinWidth',1/Settings.fs,'FaceColor',col(i,:))
% end

% Create cell arrays for table
load('DataTable')
CC_table=cell(nName,1);
genotype_table=CC_table;
targeted_table=CC_table;
method_table=CC_table;
CC_table(:)={cell(nChannel)};
lag_table=CC_table;
tempsplit=strsplit(EFF,filesep);
tempsplit=strsplit(tempsplit{end},'_');
method_table(:)=tempsplit(end);



for i=1:nName

    iTargeted=find(strcmp(uniq.Name{i},DataTable.Name)&strcmp(Channels{1},DataTable.Channel));
    targetedString=cell(nChannel,1);
    targetedString{1}=DataTable.Targeted{iTargeted};
    for ii=2:nChannel
        iTargeted=find(strcmp(uniq.Name{i},DataTable.Name)&strcmp(Channels{ii},DataTable.Channel),1,'first');
        targetedString{ii}=['_',DataTable.Targeted{iTargeted}];
    end
    
    genotype_table(i)={DataTable.Genotype{find(strcmp(uniq.Name{i},DataTable.Name),1,'first')}};
    targeted_table{i}=cat(2,targetedString{:});
    
for ROW=1:nChannel
for COL=1:nChannel
if check_matrix(ROW,COL)
    
    CC_table{i}{ROW,COL}=CC_mean{ROW,COL}(strcmp(uniq.Name{i},Name));
    lag_table{i}{ROW,COL}=lag_mean{ROW,COL}(strcmp(uniq.Name{i},Name));
    
end
end
end

end

% Save DataTable_CC to .mat file in EFF
TableNames={'Name','Genotype','Targeted','Method','CC','lag'};
DataTable_CC=table(uniq.Name,genotype_table,targeted_table,method_table,CC_table,lag_table,...
                'VariableNames',TableNames);
save('DataTable_CC.mat', 'DataTable_CC');


% for ROW=1:nChannel
% for COL=1:nChannel    
% if check_matrix(ROW,COL)
%     
%     figure
%     subplot(1,2,1)
%     hold on
%     for i=1:nEF
%         plot(t_EF{i}{ROW,COL}(1:length(CC_dis_EF{i}{ROW,COL})),CC_dis_EF{i}{ROW,COL})
%     end
%     
%     subplot(1,2,2)
%     hold on
%     for i=1:nEF
%         plot(t_EF{i}{ROW,COL}(1:length(lag_dis_EF{i}{ROW,COL})),lag_dis_EF{i}{ROW,COL})
%     end
%     
% end
% end
% end