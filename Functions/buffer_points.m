function points = buffer_points(varargin)
% Convert vector after using buffer to overlap segments into a start and 
% end point matrix.
% Input: [buffer matrix of vector, overlap length in points]
% Output: [points]
% Example: adjusted_vector=buffer_inv(rand(100,20),10);
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% assign buffered matrix or give error
if nargin<1
    error('Not enough input arguments')
else
    buffered=varargin{1};
end

% assign default number of overlap points
if nargin<2
    n_o=0;
else
    n_o=varargin{2};
end

% buffered constants
n_w=size(buffered,1); % number of points for segment
nSegment=size(buffered,2); % number of segments

if n_o<0||n_o~=floor(n_o)
    error(['Number of points of overlap (n_o>=0) should be a positive integer and '...
           'smaller than the number of points per segment minus 1 (n_o<=n_w-1)'])
end

% Get starting and end points of segments
points_start=1:n_w-n_o:(nSegment-1)*(n_w-n_o)+1;
points_end=n_w:n_w-n_o:(nSegment-1)*(n_w-n_o)+n_w;
points=[points_start',points_end'];

end