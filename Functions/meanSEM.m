function [MEAN,SEM] = meanSEM(PSDs)
% Calculate mean and SEM of a matrix
MEAN=mean(PSDs,2);
SEM=std(PSDs,0,2)/sqrt(size(PSDs,2));
end

