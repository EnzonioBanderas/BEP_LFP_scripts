function vector_adjusted = buffer_inv(varargin)
% Convert vector after using buffer to overlap segments into a adjusted
% vector using a specified method to combine overlapping segments.
% Input: [buffer matrix of vector, overlap length in points, method]
% Output: [adjusted vector]
% Example: adjusted_vector=buffer_inv(rand(100,20),10,'mean');
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% get number of input arguments
nargin=length(varargin);

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

% assing default method
if nargin<3
    method='sum';
else
    method=varargin{3};
end

% buffered constants
n_w=size(buffered,1); % number of points for segment

% Get starting and end points of segments
points=buffer_points(buffered,n_o);

% number of points of vector before buffering
n=points(end,2);

% Assign ID to each points in segments depending on position in vector before buffering
vector_ID=1:n;
vector_ID=buffer(vector_ID,n_w,n_o,'nodelay');

% List both arrays in column order for accumarray
vector_ID=vector_ID(:);
buffered=buffered(:);

% Use sum function to combine values in buffered according to their index 
% and assign them to the same index in vector_adjusted.
if any(strcmp(method,{'sum','mean'}))
    vector_adjusted=accumarray(vector_ID,buffered,[],@sum);
elseif strcmp(method,'prod')
    vector_adjusted=accumarray(vector_ID,buffered,[],@prod);
end

% Average with number of counted elements of specific ID if the method string is called 'mean'
if strcmp(method,'mean')
    vector_adjusted_n=histc(vector_ID,unique(vector_ID)); % Calculate number of time a vector ID was present in vector ID
    vector_adjusted=vector_adjusted./vector_adjusted_n; % Divide to calculate mean
end
    
end

