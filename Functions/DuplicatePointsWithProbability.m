function data_points_adjusted = DuplicatePointsWithProbability(data_points,probability_vector,nMaxPoint)
% DUPLICATEPOINTSWITHPROBABILITY duplicates points according to their
% probability. For example: A points with a probability equal to 1 for 
% nMaxPoints equal to 100 will be duplicated 100 times, whereas a points 
% with a probability of .32 will be duplicated 32 times.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
nDataPoint=size(data_points,1);
data_points_cell=cell(nDataPoint,1);
nDuplicate=round(probability_vector*nMaxPoint);
for i=1:nDataPoint
    data_points_cell{i}=repmat(data_points(i,:),[nDuplicate(i),1]);
end
data_points_adjusted=cat(1,data_points_cell{:});
end

