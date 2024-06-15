% Plots an Experiment Folder (EF). No previous script has to be run to use 
% this script.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
clearvars, close all
%% Get EF info
EF=uigetdir('','Select Experiment Folder (EF)');
EphysFiles=dir([EF,filesep,'*_CH*.continuous']);
EF_split=strsplit(EF,filesep);
EF_split=strsplit(EF_split{end},'_');
EF_name=EF_split{1};
EF_date=EF_split{2};
EF_time=EF_split{3};
EF_day=EF_split{4};
nChannel=length(EphysFiles);

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

Channels=cell(nChannel,1);
for i=1:nChannel
    Channels{i}=EphysFiles(i).name(1:end-11);
end

maxFrequency=20;

% Colors for plotting c=[blue;red;green]
c=[0     ,0.4470,0.7410;...
   0.8500,0.3250,0.0980;...
   0.4660,0.6740,0.1880];

%% Process EF data
[data,data_t,data_lin]=LFP_Process_EF(EF,Settings,Channels);
nData=length(data_t);

%% Select EF data

select_points=LFP_Select_EF(data,data_t,data_lin,EF,Settings,Channels);
close(gcf)
% select_points=select_points{1};

nonselect_points=cell(nChannel,1);
nSegment_select=cell(nChannel,1);
nSegment_nonselect=cell(nChannel,1);
for i=1:nChannel
nonselect_points{i}=logical2points((conv(~points2logical(select_points{i},nData),true(3,1),'same'))~=0);

nSegment_select{i}=size(select_points{i},1);
nSegment_nonselect{i}=size(nonselect_points{i},1);
end

%% Plot figure Traces-Spectrograms-Corrgrams
CC_check=triu(true(nChannel),1);
nCC=sum(CC_check(:));
figure('Units','normalized','OuterPosition',[0,0,0.5,1])
for i=1:nChannel
    subplot(2*nChannel+nCC,1,i), hold on, grid on
    for ii=1:nSegment_select{i}
        plot(data_t(select_points{i}(ii,1):select_points{i}(ii,2)),...
             data{i}(select_points{i}(ii,1):select_points{i}(ii,2)),'Color',c(3,:))
    end
    for ii=1:nSegment_nonselect{i}
        plot(data_t(nonselect_points{i}(ii,1):nonselect_points{i}(ii,2)),...
             data{i}(nonselect_points{i}(ii,1):nonselect_points{i}(ii,2)),'Color',c(1,:))
    end
    title(['LFP trace Channel ',Channels{i}],'Interpreter','none')
    xlabel('Time (s)'), ylabel('LFP (\muV)')
    xlim([0,data_t(end)])
end
PSD_total=cell(nChannel,1);
PSD_select=cell(nChannel,1);
for i=1:nChannel
    subplot(2*nChannel+nCC,1,i+nChannel), grid on
    
    [~,FF,TT,PSD_matrix]=spectrogram(data{i},round(Settings.window*Settings.fs),round(Settings.window*Settings.fs*Settings.overlap),round(Settings.window*Settings.fs),Settings.fs);
    
    % Calculate PSD of selected PSD trace
    PSD_select{i}=mean(spectrogram_dis(data{i},select_points{i},round(Settings.window*Settings.fs),round(Settings.window*Settings.fs*Settings.overlap),round(Settings.window*Settings.fs),Settings.fs),2);
    
    imagesc(TT,FF,PSD_matrix)
    set(gca,'YDir','normal')
    caxis([0,max(PSD_select{i})])
%     colorbar   
    
    title(['LFP PSDs/spectrogram Channel ',Channels{i}],'Interpreter','none')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    ylim([0,maxFrequency])    
    
    % Calculate PSD of total LFP trace
    PSD_total{i}=mean(PSD_matrix,2);

end
count=0;
for i=1:nChannel
    for ii=1:nChannel
        if CC_check(i,ii)
            count=count+1;
            subplot(2*nChannel+nCC,1,2*nChannel+count), grid on
            [CC,TT,LL]=corrgram2(data{i},data{ii},round(Settings.fs*0.1),round(Settings.fs*Settings.window),round(Settings.fs*Settings.window*0.5),Settings.fs);
            imagesc(TT, LL, CC)
            xlabel('Time (s)'), ylabel('Delay (s)')
            axis xy % probably same as "set(gca,'YDir','normal')"
%             colorbar
            set(gca,'CLim',[-1,1])
            
            [CC,TT,LL]=corrgram2_dis(data{i},data{ii},select_points{i},round(Settings.fs*0.1),round(Settings.fs*Settings.window),round(Settings.fs*Settings.window*0.5),Settings.fs);
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
            title(['Windowed cross correlation between channels ',Channels{i},'-',Channels{ii},', mean CC = ',num2str(mean(PKS)),' and mean delay = ',num2str(mean(LOCS)*1e3),' ms'], 'fontweight', 'bold','Interpreter','none')
            
        end
    end
end

%% Plot PSDs
figure('Units','normalized','OuterPosition',[0.5,0,0.5,1])
for i=1:nChannel
    subplot(1,2,1), hold on
    plot(FF,PSD_total{i})
    subplot(1,2,2), hold on
    plot(FF,PSD_select{i})
end

subplot(1,2,1), grid on
title('PSD computed without selection')
xlim([Settings.HP,maxFrequency])
xlabel('Frequency (Hz)'),ylabel('PSD (\muV^2)')

subplot(1,2,2), grid on
legend(Channels,'Interpreter','none')
title('PSD computed with selection')
xlim([Settings.HP,maxFrequency])
xlabel('Frequency (Hz)'),ylabel('PSD (\muV^2)')