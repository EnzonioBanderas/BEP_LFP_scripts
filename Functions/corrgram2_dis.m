function [c,t,l] = corrgram2_dis(S1,S2,points,nMaxLag,nWindowLength,nOverlap,fs)
% corrgram2 applied to discontinuously selected data where the segments are described
% by the "points" for the "S1" and "S2" inputs.
% The points inputs should be aligned in time.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
nSegment=size(points,1);
c=cell(nSegment,1);
for i=1:nSegment
    index=points(i,1):points(i,2);
    c{i}=corrgram2(S1(index),S2(index),nMaxLag,nWindowLength,nOverlap);
end
c=cat(2,c{:});

if nargout>1
t=cell(nSegment,1);
nWindow=floor((points(:,2)-points(:,1)+1-nOverlap)/(nWindowLength-nOverlap));
for i=1:nSegment
    t{i}=(floor(nWindowLength/2):nWindowLength-nOverlap:floor(nWindowLength/2)+(nWindow(i)-1)*(nWindowLength-nOverlap))';
    t{i}=t{i}+points(i,1);
    t{i}=t{i}/fs;
end
t=cat(1,t{:});
end

if nargout>2
    l=(-nMaxLag:nMaxLag)';
    l=l/fs;
end

end

