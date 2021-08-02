function [FileTable] = GetMATFile(SourceDirName,Pick)
%GETFILE 此处显示有关此函数的摘要
%   此处显示详细说明
% W203
W203_4Dir = './W203_4';
W203_5Dir = './W203_5';
% C12 TAG_1 10 序列
C12_1_4Dir = './C12_1_4';
C12_1_5Dir = './C12_1_5';
% C12 TAG_2 10 序列
C12_2_4Dir = './C12_2_4';
C12_2_5Dir = './C12_2_5';
% C12 TAG_2 PN 序列
C12_2PN_4Dir = './C12_2PN_4';
C12_2PN_5Dir = './C12_2PN_5';

if Pick == 0
    FileTable = GetAllFile(SourceDirName,'mat'); 
    return;
end

if strcmp(SourceDirName, W203_4Dir)
    FileTable = [
        % '01_01.mat';        % 1
        % '01_02.mat';
        %     '02_01.mat';
        '03_01.mat';
        '04_01.mat';        % 5
        '05_01.mat';
        % '07_01.mat';
        '08_01.mat';
        '09_01.mat';
        %     '10_01.mat';        % 10
        %     '12_01.mat';
        % '13_01.mat';
        % '14_01.mat';
        '15_01.mat';
        % '17_01.mat';        % 15
        % '18_01.mat';
          '19_01.mat';
        %     '20_01.mat';
        % '21_01.mat';
        % '22_01.mat';        % 20
        '23_01.mat';
        % '24_01.mat';
        % '25_01.mat';
        % '26_01.mat';
        % '27_01.mat';        % 25
        % '28_01.mat';
        % '28_02.mat';
        % '29_01.mat';
        '30_01.mat';
        % '32_01.mat';        % 30
        % '33_01.mat';
        '34_01.mat';
        '35_01.mat';
    ];

elseif strcmp(SourceDirName, W203_5Dir)
    FileTable = [
        % '01_01.mat';        % 1
        % '01_02.mat';
        %     '02_01.mat';
        '03_01.mat';
        '04_01.mat';        % 5
        '05_01.mat';
        % '07_01.mat';
        '08_01.mat';
        '09_01.mat';
        %     '10_01.mat';        % 10
        %     '12_01.mat';
        % '13_01.mat';
        % '14_01.mat';
        '15_01.mat';
        % '17_01.mat';        % 15
        % '18_01.mat';
          '19_01.mat';
        %     '20_01.mat';
        % '21_01.mat';
        % '22_01.mat';        % 20
        '23_01.mat';
        % '24_01.mat';
        % '25_01.mat';
        % '26_01.mat';
        % '27_01.mat';        % 25
        % '28_01.mat';
        % '28_02.mat';
        % '29_01.mat';
        '30_01.mat';
        % '32_01.mat';        % 30
        % '33_01.mat';
        '34_01.mat';
        '35_01.mat';
    ];   
    
elseif strcmp(SourceDirName, C12_1_4Dir)
    FileTable = [
        '01_01.mat';
        '02_01.mat';
        '03_01.mat';
        '04_01.mat';
        '05_01.mat';
        '06_01.mat';
        '06_02.mat';
        '07_01.mat';
        '07_02.mat';
        '08_01.mat';
        '08_02.mat';
        '09_01.mat';
        '09_02.mat';
    ];    

elseif strcmp(SourceDirName, C12_1_5Dir)
    FileTable = [
        '01_01.mat';
        '02_01.mat';
        '03_01.mat';
        '04_01.mat';
        '05_01.mat';
        '05_02.mat';
        '06_01.mat';
        '06_02.mat';
        '07_01.mat';
        '07_02.mat';
        '08_01.mat';
        '08_02.mat';
        '09_01.mat';
        '09_02.mat';
    ];  

elseif strcmp(SourceDirName, C12_2_4Dir)
    FileTable = [
        '01_01.mat';
        '01_02.mat';
        '02_01.mat';
        '03_01.mat';
        '04_01.mat';
        '05_01.mat';
        '06_01.mat';
        '07_01.mat';
        '08_01.mat';
        '09_01.mat';
    ];    

elseif strcmp(SourceDirName, C12_2_5Dir)
    FileTable = [
        '01_01.mat';
        '01_02.mat';
        '02_01.mat';
        '03_01.mat';
        '04_01.mat';
        '05_01.mat';
        '06_01.mat';
        '07_01.mat';
        '08_01.mat';
        '09_01.mat';
    ];  



else
    error('Directory has not been record !');
end
end

function [FileTable] = GetAllFile(DirName,type)
if strcmp(type,'mat')
   a = what(DirName).mat; 
   Nfile = length(a);
   for i = 1:Nfile
        FileTable(i,:) = a{i};
   end
    
end

end