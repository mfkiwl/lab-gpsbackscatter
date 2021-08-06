function llaData = ImportGTData(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  LLADATA = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  返回数值数据。
%
%  LLADATA = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  llaData = importfile("D:\Project\lab-gpsbackscatter\lab-gpsbackscatter\GPSGroudTruth\C12Circle\ReceivedTofile-COM21-2021_8_3_20-11-30.DAT", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2021-08-06 16:36:46 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 15);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = ["Var1", "Var2", "VarName3", "Var4", "VarName5", "Var6", "VarName7", "Var8", "Var9", "VarName10", "Var11", "Var12", "Var13", "Var14", "Var15"];
opts.SelectedVariableNames = ["VarName3", "VarName5", "VarName7", "VarName10"];
opts.VariableTypes = ["string", "string", "double", "string", "double", "string", "double", "string", "string", "double", "string", "string", "string", "string", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var6", "Var8", "Var9", "Var11", "Var12", "Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var6", "Var8", "Var9", "Var11", "Var12", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");

% 导入数据
llaData = readtable(filename, opts);

%% 转换为输出类型
llaData = table2array(llaData);
end