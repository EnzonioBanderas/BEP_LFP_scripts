function [RESHAPED,TRUNCATED,POINTS] = reshape2(SHAPED,DIMENSION)
% Obsolete function (buffered works better but requires signal processing toolbox)
%
% Reshape array according to defined dimensions and truncate if necessary.
% Save truncated array.
% Save index points of reshaped and truncated sequences.
% (First output is the same as the output of the "reshaped" function when there is no truncation)
% Input: SHAPED sequence, DIMENSIONs of RESHAPED sequence ([DIMENSION(1),DIMENSION(2),...])
% Output: RESHAPED sequence, TRUNCATED sequence, start and end POINTS.

% Get number of elements of shaped, reshaped and truncated sequences
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
nShaped=numel(SHAPED);
nReshaped=prod(DIMENSION);
nTruncated=nShaped-nReshaped;

if nTruncated<0
    error('Number of elements of reshaped sequence higher than number of elements of shaped sequence.')
end

% Save truncated sequence and truncate shaped
SHAPED=SHAPED(:);
TRUNCATED=SHAPED(end-nTruncated+1:end);
SHAPED=SHAPED(1:end-nTruncated);

% Reshape shaped sequence 
RESHAPED=reshape(SHAPED,DIMENSION);

% Get start and end index points of reshaped (1:end-1,[1,2]) and truncated (end,[1,2]) sequences
points_start=([1:DIMENSION(1):nReshaped,nReshaped+1])';
points_end=([DIMENSION(1):DIMENSION(1):nReshaped,nShaped])';
POINTS=[points_start,points_end];

% If there is no truncation remove points representing truncated sequence
if nTruncated==0
    POINTS(end,:)=[];
end

end