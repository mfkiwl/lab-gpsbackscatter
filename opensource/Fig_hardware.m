%% 导入电子表格中的数据
% 用于从以下电子表格导入数据的脚本:
%
%    工作簿: D:\Project\lab-gpsbackscatter\dhx_paper\Hardware\Hardware\gpsbackscatter.xlsx
%    工作表: Sheet1
%
% 由 MATLAB 于 2021-03-10 16:28:09 自动生成

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 20);

% 指定工作表和范围
opts.Sheet = "Sheet1";
opts.DataRange = "A2:T32";

% 指定列名称和类型
opts.VariableNames = ["f", "a", "VarName3", "VmV", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 导入数据
gpsbackscatter1 = readtable("D:\Project\lab-gpsbackscatter\dhx_paper\Hardware\Hardware\gpsbackscatter.xlsx", opts, "UseExcel", false);

%% 转换为输出类型
gpsbackscatter1 = table2array(gpsbackscatter1);

%% 清除临时变量
clear opts

%% 分离数据
jF=1;jVol=(4:2:20);jPower=[2 (5:2:19)];
NumF = size(gpsbackscatter,1);
NumVol = size(jVol,2);

f = gpsbackscatter(:,jF);
Power = gpsbackscatter(:,jPower);
Gain = Power - (-133);
Vol = gpsbackscatter(:,jVol);

figure
plot( f * ones(1, NumVol),Gain);
legend('120','125','130','135','140','145','150','155','160')







