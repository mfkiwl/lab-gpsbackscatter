function [FileTable, Nfile] = GetAllFile(DirName,type)
if strcmp(type,'DAT') 
    File = dir(fullfile(DirName,'*.DAT')); 
    FileTable = {File.name}';               
    Nfile = length(FileTable);
    
elseif strcmp(type,'txt') 
    File = dir(fullfile(DirName,'*.txt')); 
    FileTable = {File.name}';               
    Nfile = length(FileTable);
end
end