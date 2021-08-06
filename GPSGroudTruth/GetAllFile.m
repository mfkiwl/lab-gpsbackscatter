function [FileTable, Nfile] = GetAllFile(DirName,type)
if strcmp(type,'DAT') 
    File = dir(fullfile(DirName,'*.DAT')); 
    FileTable = {File.name}';               
    Nfile = length(FileTable);
end
end