function combine_datatables
% Selects *Table*.mat files in the same folder and combines them by
% concatenation.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com

% Select to be analyzed *Table*.mat files
[FileName,PathName] = uigetfile('*Table*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(FileName) %if-statement for the case that only one file is selected
    error('Not enough DataTables selected')
end
cd(PathName)

% Start looping through selected DataTables
nTable=length(FileName);
temp_cell=cell(nTable,1);
for i=1:nTable

    % load DataTable, assume that only one DataTable is present per file
    load(FileName{i},'*Table*')
    
    % list workspace variables, assume that tables within .mat files have si
    TableName=who('*Table*');
    eval(['temp_cell{i}=',TableName{1},';'])

end

% concatenate tables and assign to variable with the same name as loaded
% files
eval([TableName{1},'=cat(1,(temp_cell{:}));'])

% Save combined table
save([TableName{1},'_',GetTime],TableName{1})

end

