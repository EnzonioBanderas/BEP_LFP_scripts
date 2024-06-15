function data_buffered = buffer_dis(data,nWindowLength,nOverlap,points)
% BUFFER_DIS does the same thing as buffer but for segments
% defined by points input.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% Check number of input and output arguments
narginchk(4,4)

% preallocation
nSegment=size(points,1);
data_buffered=cell(nSegment,1);

% buffer segments of data
for i=1:nSegment
    index=points(i,1):points(i,2);
    [data_buffered{i},~]=buffer(data(index),nWindowLength,nOverlap,'nodelay');
end

% concatenate cells
data_buffered=cat(2,data_buffered{:});

end

