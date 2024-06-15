% Plots selected .continuous files without asking for selection of data. 
% No previous script has to be run to use this script.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Get LFP recording (...CH.continuous file)
[EphysFiles,PathName] = uigetfile('*.continuous','Select all the files to analyse','MultiSelect', 'on');
%if-statement for the case that there is only one channel
if ischar(EphysFiles) 
    EphysFiles={EphysFiles}; 
end 
nEphysFiles=length(EphysFiles);
cd(PathName)

%% Setttings
Settings.fs=1000; 
Settings.Detrend=true; 
Settings.window_detrend=2; 
Settings.overlap_detrend=0.5; 
Settings.HP=1; 
Settings.LP=100; 

Settings.SelectMultipleChannels=true;
Settings.ShowLinearFit=false;
Settings.Randomize=false;
Settings.time_threshold=5;
Settings.window=2;
Settings.overlap=0.5;

Resample=true;

EphysFilesNames=cell(nEphysFiles,1);
for i=1:nEphysFiles
    EphysFilesNames{i}=EphysFiles{i}(1:end-11);
end

maxFrequency=20;

% Colors for plotting c=[blue;red;green]
c=[0     ,0.4470,0.7410;...
   0.8500,0.3250,0.0980;...
   0.4660,0.6740,0.1880];

%% Load and process data
% Preallocate data cells
data=cell(nEphysFiles,1);
info=data;
data_lin=data;

for i=1:nEphysFiles
    % Load OpenEhys .continuous files
    [data{i},~,info{i}] = load_open_ephys_data_faster(EphysFiles{i});
end

% get ends of experiment
data_t_end=zeros(nEphysFiles,1);
for i=1:nEphysFiles
    data_t_end(i)=length(data{i});
end
data_t_end_min=min(data_t_end);

% truncate data to make sure that they have the same lengths
for i=1:nEphysFiles
    data_t_end_trunc=data_t_end(i)-data_t_end_min;
    data{i}(end-data_t_end_trunc+1:end)=[];
end

for i=1:nEphysFiles
    
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

%% No selection of data
select_points=[1,length(data{1})];
nonselect_points=[];

nSegment_select=size(select_points,1);
nSegment_nonselect=size(nonselect_points,1);

%% Plot figure Traces-Spectrograms-Corrgrams
CC_check=triu(true(nEphysFiles),1);
nCC=sum(CC_check(:));
figure('Units','normalized','OuterPosition',[0,0,0.5,1])
for i=1:nEphysFiles
    subplot(2*nEphysFiles+nCC,1,i), hold on, grid on
    for ii=1:nSegment_select
        plot(data_t(select_points(ii,1):select_points(ii,2)),...
             data{i}(select_points(ii,1):select_points(ii,2)),'Color',c(3,:))
    end
    for ii=1:nSegment_nonselect
        plot(data_t(nonselect_points(ii,1):nonselect_points(ii,2)),...
             data{i}(nonselect_points(ii,1):nonselect_points(ii,2)),'Color',c(1,:))
    end
    title(['LFP trace Channel ',EphysFilesNames{i}],'Interpreter','none')
    xlabel('Time (s)'), ylabel('LFP (\muV)')
    xlim([0,data_t(end)])
end
PSD_total=cell(nEphysFiles,1);
for i=1:nEphysFiles
    subplot(2*nEphysFiles+nCC,1,i+nEphysFiles), grid on
    
    [~,FF,TT,PSD_matrix]=spectrogram(data{i},round(Settings.window*Settings.fs),round(Settings.window*Settings.fs*Settings.overlap),round(Settings.window*Settings.fs),Settings.fs);
    
    % Calculate PSD of total LFP trace
    PSD_total{i}=mean(PSD_matrix,2);
    
    imagesc(TT,FF,PSD_matrix)
    set(gca,'YDir','normal')
    caxis([0,max(PSD_total{i})])
%     colorbar   
    
    title(['LFP PSDs/spectrogram Channel ',EphysFilesNames{i}],'Interpreter','none')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    ylim([0,maxFrequency])    

end
count=0;
for i=1:nEphysFiles
    for ii=1:nEphysFiles
        if CC_check(i,ii)
            count=count+1;
            subplot(2*nEphysFiles+nCC,1,2*nEphysFiles+count), grid on
            [CC,TT,LL]=corrgram2(data{i},data{ii},round(Settings.fs*0.1),round(Settings.fs*Settings.window),round(Settings.fs*Settings.window*0.5),Settings.fs);
            imagesc(TT, LL, CC)
            xlabel('Time (s)'), ylabel('Delay (s)')
            axis xy % probably same as "set(gca,'YDir','normal')"
%             colorbar
            set(gca,'CLim',[-1,1])
            
            [CC,TT,LL]=corrgram2_dis(data{i},data{ii},select_points,round(Settings.fs*0.1),round(Settings.fs*Settings.window),round(Settings.fs*Settings.window*0.5),Settings.fs);
            nWindow=length(TT); PKS=zeros(nWindow,1); LOCS=zeros(nWindow,1);
            for iii=1:nWindow
                [pks,locs]=findpeaks(CC(:,iii),LL,'WidthReference','halfheight');
                if isempty(pks)
                    [PKS(iii),index]=max(CC(:,iii));
                    LOCS(iii)=LL(index);
                else
                    [~,index]=min(abs(locs));
                    LOCS(iii)=locs(index);
                    PKS(iii)=pks(index);
                end
            end
            title(['Windowed cross correlation between channels ',EphysFilesNames{i},'-',EphysFilesNames{ii},', mean CC = ',num2str(mean(PKS)),' and mean delay = ',num2str(mean(LOCS)*1e3),' ms'], 'fontweight', 'bold','Interpreter','none')
            
        end
    end
end

%% Plot PSDs
figure('Units','normalized','OuterPosition',[0.5,0,0.5,1]), hold on, grid on
for i=1:nEphysFiles
    plot(FF,PSD_total{i})
end

title('PSD computed without selection')
xlim([Settings.HP,maxFrequency])
xlabel('Frequency (Hz)'),ylabel('PSD (\muV^2)')
legend(EphysFilesNames,'Interpreter','none')