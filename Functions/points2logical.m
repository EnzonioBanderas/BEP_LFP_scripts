function logical = points2logical(varargin)
% Take start and end point matrix and convert it to a logical possibly padded
% with false values if second input is given.
% Inputs: point matrix, length in points of to be created logical
% Output: logical
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% Assign variables from input
narginchk(1,2)
points=varargin{1};
if nargin>1
    logical_length=varargin{2};
else
    logical_length=1;
end

if isempty(points)
    logical=false(logical_length,1);
else
    
if nargin<2
    logical_length=points(end,2);
else
    if logical_length<points(end,2)
        error('Defined length shorter than last point')
    end
end

logical=false(logical_length,1);
for j=1:size(points,1)
logical(points(j,1):points(j,2))=true;
end

end

end

