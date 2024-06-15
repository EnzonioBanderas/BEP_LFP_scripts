function current_time_chars = GetTime
% Get the current time in a "2018-05-12_18-55-56" format
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
current_time=fix(clock);
current_time_cell=cell(6,1);
for i=1:6
    current_time_cell{i}=num2str(current_time(i));
end
% current_time_cell{1}=current_time_cell{1}(end-1:end);

for i=2:6
if length(current_time_cell{i})<2
    current_time_cell{i}=['0',current_time_cell{i}];
end
end

for i=[2,3,5,6]
    current_time_cell{i}=['-',current_time_cell{i}];
end
current_time_chars=[[current_time_cell{1:3}],'_',[current_time_cell{4:6}]];

end

