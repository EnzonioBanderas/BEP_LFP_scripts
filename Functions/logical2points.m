function points = logical2points(logical)
% Convert a logical vector to a matrix that contains the start and end points
% in the first and second columns respectively for each consecutive sequence of true values.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% if the logical vector is 1xn instead of nx1
if size(logical,2)>1
    logical=logical';
end

% compute points by differencing, finding indices and shifting end points
logical=[false;logical;false];
logical_diff=diff(logical);
points_start=find(logical_diff==1);
points_end=find(logical_diff==-1)-1;
points=[points_start,points_end];

end

