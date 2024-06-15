% Select Experiment Folder Folder Folder (EFFF) containing Experiment Folder 
% Folders (EFFs) containing Experiment Folders (EFs).
% Plots cross-correlation coefficients and lags calculated with LFP_CC_EFF.
% Saves a .mat file with all mean cross-correlation and lag information of EFFF.
% Four figures for each method:
% 1) Mean cross-correlation and delay/lag per day for each mouse
% 2) Mean cross-correlation and delay/lag for each measurement seperated according to genotype
% 3) Mean cross-correlation and delay/lag per mouse seperated according to genotype
% 4) Mean cross-correlation and delay/lag per mouse seperated according to whether 
%    one of the sites of electrode placement was targeted or not (both sides targeted are not included)
% Assumption: Settings.mat and Channels.mat are the same across EFFs.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all

% Select Experiment Folder Folder Folder (EFFF) and 
% create list of Experiment Folder Folder (EFF) names
EFFF = uigetdir('','Select folder containing experiments');
cd (EFFF)
EFF=dir('*_*_*');
nEFF=length(EFF);
EFF_Name={EFF(:).name}';

% Load Settings and Channels assuming that they are constant across EFFs
cd(EFF_Name{1})
load('Settings')
load('Channels')
nChannel=size(Channels,1);
cd ..

% Plot settings
col=colormap(lines);
close(gcf)

% Load all DataTable_CC's in EFFs and concatenate
DataTable_CC_cell=cell(nEFF,1);
DataTable_cell=cell(nEFF,1);
for i=1:nEFF
    cd(EFF_Name{i})
    
    load('DataTable_CC')
    DataTable_CC_cell{i}=DataTable_CC;
    
    load('DataTable')
    DataTable_cell{i}=DataTable;
    
    cd ..
end
DataTable_CC=cat(1,DataTable_CC_cell{:});
DataTable=cat(1,DataTable_cell{:});

% % Targeted input for each channel
% % Define Targeted per Mouse Name (Channels{1}_Channels{2} == yes_no)
% uniq.Name=unique(DataTable.Name);
% nName=length(uniq.Name);
% strcmp(DataTable.Channel,Channels{1});
% for i=1:nName
%     logical.Name=strcmp(DataTable.Name,uniq.Name{i});
% for ii=1:nChannel
%     
%     
%     
% end
% end
% 
% uniq.Channel_full=unique(Channel_string);
% nChannel_full=length(uniq.Channel_full);
% 
% Targeted_def_temp=strings(nName,1);
% Targeted_definput='no';
% Targeted_def_temp{1}=Targeted_definput;
% for i=2:nChannel_full
%     Targeted_def_temp{i}=['_',Targeted_definput];
% end
% Targeted_def=strings(length(uniq.Name),1);
% Targeted_def(:)=convertCharsToStrings([Targeted_def_temp{:}]);
% 
% uniq.Targeted=inputdlg(uniq.Name,['Enter Targeted ',[uniq.Channel_full{:}]],1,Targeted_def);
% 
% Targeted_string=strings(nTableRow,1);
% for i=1:nName
%     tempstring=strsplit(uniq.Targeted{i},'_');
%     logical.name=strcmp(Mouse_Name_string,uniq.Name(i));
%     uniq.Channel=unique(Channel_string(logical.name));
%     nChannel=length(uniq.Channel);
%     for ii=1:nChannel
%         tempindex=find(strcmp(uniq.Channel_full,uniq.Channel(ii)));
%         logical.name_channel=logical.name&strcmp(Channel_string,uniq.Channel(ii));
%         Targeted_string(logical.name_channel)...
%         =convertCharsToStrings(tempstring{tempindex});
%     end
% end

% Predefine figures
uniq.Method_Full=unique(DataTable_CC.Method);
nMethod_Full=length(uniq.Method_Full);
for i=1:nMethod_Full
    fig1.Method{i}=figure('Name',['PerNameEachMeasurement_',uniq.Method_Full{i}]);
    fig2.Method{i}=figure('Name',['PerGenotypeEachMeasurement_',uniq.Method_Full{i}]);
    fig3.Method{i}=figure('Name',['PerGenotype_',uniq.Method_Full{i}]);
    fig4.Method{i}=figure('Name',['For Targeted and Non-Targeted ',uniq.Method_Full{i}]);
end

% simple check matrix so that no cross-correlations are computed
% twice and that no autocorrelations are computed
check_matrix=triu(true(nChannel),1);

% Concatenate possible duplicate name-method's for both CC and lag -day
% distributions
% Plot distribution CC and lag for different colors for each name
legend_cell=cell(nMethod_Full,1);
count=0;
uniq.Name_Full=unique(DataTable_CC.Name);
nName_Full=length(uniq.Name_Full);
for i=1:nName_Full
    logical.Name=strcmp(uniq.Name_Full{i},DataTable_CC.Name);
    uniq.Method=unique(DataTable_CC.Method(logical.Name));
    nMethod=length(uniq.Method);
    for ii=1:nMethod
        logical.Name_Method=logical.Name&strcmp(uniq.Method{ii},DataTable_CC.Method);
        CC=DataTable_CC.CC(logical.Name_Method);
        lag=DataTable_CC.lag(logical.Name_Method);
        nDuplicate=sum(logical.Name_Method);
        
        ind.Name=i;
        ind.Method=find(strcmp(uniq.Method{ii},uniq.Method_Full));
        
        for ROW=1:nChannel
        for COL=1:nChannel
        if check_matrix(ROW,COL)
        count=count+1;
            
            % Merge
            CC_cell=cell(nDuplicate,1);
            lag_cell=cell(nDuplicate,1);
            for iii=1:nDuplicate
                CC_cell{iii}=CC{iii}{ROW,COL};
                lag_cell{iii}=lag{iii}{ROW,COL};
            end
            CC_hist=cat(1,CC_cell{:});
            lag_hist=cat(1,lag_cell{:});
                
                % plot
                set(0,'CurrentFigure',fig1.Method{ind.Method})
                % plot CC
                subplot(1,2,1)
                hold on
                grid on
                histogram(CC_hist,'FaceColor',col(ind.Name,:),'BinWidth',.05)
                title('cross-correlation PerNameEachMeasurement')
                xlabel('cross-correlation'),ylabel('Count of windowed cross-correlation closest to 0 peaks')
                % plot lag
                subplot(1,2,2)
                hold on
                grid on
                histogram(lag_hist,'FaceColor',col(ind.Name,:),'BinWidth',1/Settings.fs)
                title('lag PerNameEachMeasurement')
                xlabel('lag (s)'),ylabel('Count of windowed cross-correlation closest to 0 peaks')
                % legend
                legend_cell{ind.Method}=[legend_cell{ind.Method},uniq.Name_Full(i)];
                
            % Mean and mode method-names for later use for PerGenotype
            % figure
            Genotype_cell{count}=DataTable_CC.Genotype{find(logical.Name_Method,1,'first')};
            targeted_split=strsplit(DataTable_CC.Targeted{find(logical.Name_Method,1,'first')},'_');
            Targeted_cell{count}=[targeted_split{ROW},'_',targeted_split{COL}];
            Method_cell{count}=uniq.Method{ii};
            
            CC_mean_cell{count}{ROW,COL}=mean(CC_hist);
            lag_mean_cell{count}{ROW,COL}=mean(lag_hist);
        
        end
        end
        end
        
    end
end
for i=1:nMethod_Full
    set(0,'CurrentFigure',fig1.Method{i})
    subplot(1,2,1)
    legend(legend_cell{i},'Interpreter','none')
    subplot(1,2,2)
    legend(legend_cell{i},'Interpreter','none')
end

% Plot same thing by concatentating genotype-methods (each measurement will be visible)
% different colors here for each genotype
legend_cell=cell(nMethod_Full,1);
uniq.Genotype_Full=unique(DataTable_CC.Genotype);
nGenotype_Full=length(uniq.Genotype_Full);
for i=1:nGenotype_Full
    logical.Genotype=strcmp(uniq.Genotype_Full{i},DataTable_CC.Genotype);
    uniq.Method=unique(DataTable_CC.Method(logical.Genotype));
    nMethod=length(uniq.Method);
    for ii=1:nMethod
        logical.Genotype_Method=logical.Genotype&strcmp(uniq.Method{ii},DataTable_CC.Method);
        CC=DataTable_CC.CC(logical.Genotype_Method);
        lag=DataTable_CC.lag(logical.Genotype_Method);
        nDuplicate=sum(logical.Genotype_Method);
        
        count=count+1;
        ind.Genotype=i;
        ind.Method=find(strcmp(uniq.Method{ii},uniq.Method_Full));
        
        for ROW=1:nChannel
        for COL=1:nChannel
        if check_matrix(ROW,COL)
            
            % Merge
            CC_cell=cell(nDuplicate,1);
            lag_cell=cell(nDuplicate,1);
            for iii=1:nDuplicate
                CC_cell{iii}=CC{iii}{ROW,COL};
                lag_cell{iii}=lag{iii}{ROW,COL};
            end
            CC_hist=cat(1,CC_cell{:});
            lag_hist=cat(1,lag_cell{:});
                
                % plot
                set(0,'CurrentFigure',fig2.Method{ind.Method})
                % plot CC
                subplot(1,2,1)
                hold on
                grid on
                histogram(CC_hist,'FaceColor',col(ind.Genotype,:),'BinWidth',.05)
                title('cross-correlation PerGenotypeEachMeasurement')
                xlabel('cross-correlation'),ylabel('Count of windowed cross-correlation closest to 0 peaks')
                % plot lag
                subplot(1,2,2)
                hold on
                grid on
                histogram(lag_hist,'FaceColor',col(ind.Genotype,:),'BinWidth',1/Settings.fs)
                title('lag PerGenotypeEachMeasurement')
                xlabel('lag (s)'),ylabel('Count of windowed cross-correlation closest to 0 peaks')
                % legend
                legend_cell{ind.Method}=[legend_cell{ind.Method},uniq.Genotype_Full(i)];
        
        end
        end
        end
        
    end
end
for i=1:nMethod_Full
    set(0,'CurrentFigure',fig2.Method{i})
    subplot(1,2,1)
    legend(legend_cell{i},'Interpreter','none')
    subplot(1,2,2)
    legend(legend_cell{i},'Interpreter','none')
end


% Take mean and mode of CC and lag respectively for each name and plot in
% similar way
legend_cell=cell(nMethod_Full,1);
uniq.Genotype_Full=unique(Genotype_cell);
nGenotype_Full=length(uniq.Genotype_Full);
for i=1:nGenotype_Full
    logical.Genotype=strcmp(uniq.Genotype_Full{i},Genotype_cell);
    uniq.Method=unique(Method_cell(logical.Genotype));
    nMethod=length(uniq.Method);
    for ii=1:nMethod
        logical.Genotype_Method=logical.Genotype&strcmp(uniq.Method{ii},Method_cell);
        CC=CC_mean_cell(logical.Genotype_Method);
        lag=lag_mean_cell(logical.Genotype_Method);
        nDuplicate=sum(logical.Genotype_Method);
        
        count=count+1;
        ind.Genotype=i;
        ind.Method=find(strcmp(uniq.Method{ii},uniq.Method_Full));
        
        for ROW=1:nChannel
        for COL=1:nChannel
        if check_matrix(ROW,COL)
            
            % Merge
            CC_cell=cell(nDuplicate,1);
            lag_cell=cell(nDuplicate,1);
            for iii=1:nDuplicate
                CC_cell{iii}=CC{iii}{ROW,COL};
                lag_cell{iii}=lag{iii}{ROW,COL};
            end
            CC_hist=cat(1,CC_cell{:});
            lag_hist=cat(1,lag_cell{:});
                
                % plot
                set(0,'CurrentFigure',fig3.Method{ind.Method})
                % plot CC
                subplot(1,2,1)
                hold on
                grid on
                histogram(CC_hist,'FaceColor',col(ind.Genotype,:),'BinWidth',.05)
                title('cross-correlation PerGenotype')
                xlabel('cross-correlation'),ylabel('Count of mice')
                % plot lag
                subplot(1,2,2)
                hold on
                grid on
                histogram(lag_hist,'FaceColor',col(ind.Genotype,:),'BinWidth',1/Settings.fs)
                title('lag PerGenotype')
                xlabel('lag (s)'),ylabel('Count of mice')
                % legend
                legend_cell{ind.Method}=[legend_cell{ind.Method},uniq.Genotype_Full(i)];
                
        end
        end
        end
        
    end
end
for i=1:nMethod_Full
    set(0,'CurrentFigure',fig3.Method{i})
    subplot(1,2,1)
    legend(legend_cell{i},'Interpreter','none')
    subplot(1,2,2)
    legend(legend_cell{i},'Interpreter','none')
end

%% Histogram plots Targeted-Nontargeted combinations
% Take mean and mode of CC and lag respectively for each name and plot in
% similar way
legend_cell=cell(nMethod_Full,1);
uniq.Targeted={'no_no';'yes_no or no_yes'};
nTargeted=length(uniq.Targeted);
for i=1:nTargeted
    if i==1
        logical.Targeted=strcmp('no_no',Targeted_cell);
    end
    if i==2
        logical.Targeted=strcmp('yes_no',Targeted_cell)|strcmp('no_yes',Targeted_cell);
    end
    uniq.Method=unique(Method_cell(logical.Targeted));
    nMethod=length(uniq.Method);
    for ii=1:nMethod
        logical.Targeted_Method=logical.Targeted&strcmp(uniq.Method{ii},Method_cell);
        CC=CC_mean_cell(logical.Targeted_Method);
        lag=lag_mean_cell(logical.Targeted_Method);
        nDuplicate=sum(logical.Targeted_Method);
        
        count=count+1;
        ind.Targeted=i;
        ind.Method=find(strcmp(uniq.Method{ii},uniq.Method_Full));
        
        for ROW=1:nChannel
        for COL=1:nChannel
        if check_matrix(ROW,COL)
            
            % Merge
            CC_cell=cell(nDuplicate,1);
            lag_cell=cell(nDuplicate,1);
            for iii=1:nDuplicate
                CC_cell{iii}=CC{iii}{ROW,COL};
                lag_cell{iii}=lag{iii}{ROW,COL};
            end
            CC_hist=cat(1,CC_cell{:});
            lag_hist=cat(1,lag_cell{:});
                
                % plot
                set(0,'CurrentFigure',fig4.Method{ind.Method})
                % plot CC
                subplot(1,2,1)
                hold on
                grid on
                histogram(CC_hist,'FaceColor',col(ind.Targeted,:),'BinWidth',.05)
                title('cross-correlation For Targeted and Non-Targeted')
                xlabel('cross-correlation'),ylabel('Count of unique cross correlation pairs')
                % plot lag
                subplot(1,2,2)
                hold on
                grid on
                histogram(lag_hist,'FaceColor',col(ind.Targeted,:),'BinWidth',1/Settings.fs)
                title('lag For Targeted and Non-Targeted')
                xlabel('lag (s)'),ylabel('Count of unique cross correlation pairs')
                % legend
                legend_cell{ind.Method}=[legend_cell{ind.Method},uniq.Targeted(i)];
                
        end
        end
        end
        
    end
end
for i=1:nMethod_Full
    set(0,'CurrentFigure',fig4.Method{i})
    subplot(1,2,1)
    legend(legend_cell{i},'Interpreter','none')
    subplot(1,2,2)
    legend(legend_cell{i},'Interpreter','none')
end

%% Save DataTable_CC
EFFF=strsplit(EFFF,filesep);
EFFF=EFFF{end};
save(['DataTableCC_',EFFF],'DataTable_CC')