% Select Experiment Folder Folder Folder (EFFF) containing Experiment
% Folder Folders (EFFs) containing Experiment Folders (EFs). Can be run
% after each EFF is compiled and has a DataTable.mat file in it. Settings
% defined in LFP_Analyzer_EFF should be the same across EFFs.
clearvars, close all
% Get names of EFFF and EFFs
EFFF=uigetdir('','Select folder containing folders containing experiment folders');
EFF=dir([EFFF,filesep,'*_*_*']);
nEFF=length(EFF);

% Assume Settings are the same over EFFs and define column names
load([EFFF,filesep,EFF(1).name,filesep,'Settings.mat'])
frequencyBins=0:1/Settings.window:Settings.fs/2;
columnNames=[{'Genotype','Method','Channel','Name','Day','Targeted','Time_Length (s)'},num2cell(frequencyBins)];

% Convert DataTable.mat to cell array
DataSheetCell=cell(nEFF,1);
for i=1:nEFF
    load([EFFF,filesep,EFF(i).name,filesep,'DataTable.mat'])
    DataSheetCell{i}=[{DataTable.Genotype{:}}',{DataTable.Method{:}}',...
                      {DataTable.Channel{:}}',{DataTable.Name{:}}',...
                      {DataTable.Day{:}}',{DataTable.Targeted{:}}',...
                      {DataTable.Time_Length{:}}',num2cell([DataTable.Power{:}]')];
end
DataCellArray=[columnNames;cat(1,DataSheetCell{:})];

% Convert cell array to excel file
xlswrite([EFFF,filesep,'DataTableExcel.xlsx'],DataCellArray)