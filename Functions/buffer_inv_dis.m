function data_combined = buffer_inv_dis(data_buffered,nOverlap,Method,points,data_length,fill)
% BUFFER_INV_DIS: buffer_inv but for discontinuous segments described by
% points. Optional data_length parameter which makes sure that the vector
% is assigned to a vector of length data_length filled with 0's or NaN's
% depending on the fill parameter.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% input check
narginchk(4,6)
if nargin<5||isempty(data_length)
    data_length=points(end,2);
end
if nargin<6
    fill=0;
end
if fill==0
    data_combined=zeros(data_length,1);
elseif isnan(fill)
    data_combined=NaN(data_length,1);
end

% preallocation
nSegment=size(points,1);

% points that describe which windows belong to which segments
nWindowLength=size(data_buffered,1);
nWindow=floor((points(:,2)-points(:,1)+1-nOverlap)/(nWindowLength-nOverlap));
WindowCount=0;
points_window=zeros(nSegment,2);
for i=1:nSegment
    WindowCount=WindowCount+nWindow(i);
    points_window(i,2)=WindowCount;
end
points_window(1,1)=1;
for i=2:nSegment
    points_window(i,1)=points_window(i-1,2)+1;
end

% assign buffered data after buffer_inv into combined data
for i=1:nSegment
    index_buffered=points_window(i,1):points_window(i,2);
    data_temp=buffer_inv(data_buffered(:,index_buffered),nOverlap,Method);
    index=points(i,1):points(i,1)+length(data_temp)-1;
    data_combined(index)=data_temp;
end

end

