function [c,t,l,CC_dis,lag_dis] = LFP_CC_EF(EF_Name,Settings,Channels,nMaxLag,varargin)
% LFP_CC_EF calculates cross-correlation coefficients for different time 
% delays between channels in an Experiment Folder (EF).
% [C,T,L,C_DIS,LAG_DIS] = LFP_CC_EF(EF_NAME,SETTINGS,CHANNELS,NMAXLAG,INCLUDEALL,PLOT,FULLSELECTION)
% LFP_CC_EF takes as input EF_NAME which is the name of the EF, 
% SETTINGS and CHANNELS which are the same as defined in
% LFP_Analyzer_EFF, NMAXLAG which is the maximum number of points of delay 
% over which the cross-correlation is computed, INCLUDEALL which is a
% logical which if true includes autocorrelations and 'double'
% cross-correlations in calculations, PLOT which if true plots corrgram's
% and FULLSELECTION which if true changes the selection made with
% LFP_Analyzer_EFF to a full selection of the whole trace. Outputs are
% the cross-correlation coefficients C (c), the central time of the window in seconds T (x),
% the delay vector in seconds L (y), the cross correlation of each window of the peak
% closest to 0 (CC_DIS) and the corresponding delay/lag for this peak (LAG_DIS).
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% Check and parse inputs
narginchk(4,7)
if nargin<5
    IncludeAll=false;
else
    IncludeAll=varargin{1};
end
if nargin<6
    Plot=false;
else
    Plot=varargin{2};
end
if nargin<7
    FullSelection=false;
else
    FullSelection=varargin{3};
end

% Use settings to define parameters
nWindowLength=Settings.window*Settings.fs;
nOverlap=nWindowLength*Settings.overlap;

% Get number of channels
nChannel=size(Channels,1);
if nChannel<2
    error('cross-correlation requires at least 2 channels')
end

% Go into Experiment Folder
cd(EF_Name)

% Load data
data=cell(nChannel,1);
select_points=data;
load('EF_data','data','select_points')
if length(data)~=nChannel
    error('Channel number in EF_data does not match defined channel number')
end

% select all
if FullSelection
    for i=1:nChannel
    select_points{i}=[1,length(data{i})];
    end
end

CHECK=false;
for i=2:nChannel
    check=select_points{i-1}~=select_points{i};
    if any(check(:))
        CHECK=true;
        disp('Selection is not equal for all channels')
    end
end
if CHECK
    select_logical=cell(nChannel,1);
    for i=1:nChannel
        select_logical{i}=points2logical(select_points{i},length(data{i}));
    end
    for i=2:nChannel
        select_logical{i}=select_logical{i}&select_logical{i-1};
    end
    for i=1:nChannel
    select_points{i}=logical2points(select_logical{nChannel});
    select_points{i}=select_points{i}(select_points{i}(:,2)-select_points{i}(:,1)+1>=nWindowLength+2*nMaxLag,:);
    if isempty(select_points{i})
        error('Recombined selection does not contain enough points for cross-correlation')
    end
    end
end

% preallocation
c=cell(nChannel);
t=c;
l=c;
CC_dis=c;
lag_dis=c;
if IncludeAll
    cross_matrix_logical=true(nChannel);
else
    cross_matrix_logical=triu(true(nChannel),1);
end

if Plot
    figure('Name',[EF_Name,'_corrgram']);
end

for i=1:nChannel
    for ii=1:nChannel
        if cross_matrix_logical(i,ii)
            
            [c{i,ii},t{i,ii},l{i,ii}]=corrgram2_dis(data{i},data{ii},select_points{1},nMaxLag,nWindowLength,nOverlap,Settings.fs);
            
            if Plot
                if IncludeAll
                    subplot(nChannel,nChannel,ii+(i-1)*nChannel)
                else
                    subplot(nChannel-1,nChannel-1,(ii-1)+(i-1)*(nChannel-1))
                end
                imagesc(t{i,ii},l{i,ii},c{i,ii})
                grid on
                title([Channels{i,1},'-',Channels{ii,1}])
                hcb=colorbar;
                hcb.Label.String='cross-correlation';
                set(gca,'CLim',[-1,1])
                xlabel('Time (s)'), ylabel('Lag (s)')
            end
            
            nWindow=size(c{i,ii},2);
            lag=zeros(nWindow,1);
            CC=zeros(nWindow,1);
            for iii=1:nWindow
                [pks,locs,~,~]=findpeaks(c{i,ii}(:,iii),l{i,ii},'WidthReference','halfheight');
                if ~isempty(pks)
                [~,peak_closestto0_index]=min(abs(locs));
                lag(iii)=locs(peak_closestto0_index);
                CC(iii)=pks(peak_closestto0_index);
                else
                lag(iii)=l{i,ii}(nMaxLag+1);
                CC(iii)=c{i,ii}(nMaxLag+1,iii);
                end
            end
            lag_dis{i,ii}=lag;
            CC_dis{i,ii}=CC;
            
        end
    end
end

cd ..

end

