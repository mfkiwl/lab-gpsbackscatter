close all;
clear 
% 需要设定好 
% 1.读取原始 TXT 文件的路径
% 2.创建输出文件（距离误差，hdop）的保存路径，并设定好
% 3.GroundTruth 文件放入到输出文件夹下，设定好总的点数
% 4.设定参与计算的卫星数目

W203Dir = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
% C12 TAG_1 10 序列
C12Dir_1 = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc1';
% C12 TAG_2 10 序列
C12Dir_2 = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc2';
% C12 TAG_2 PN 序列
C12Dir_2PN = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\C12ES\0728测试\TagLoc2_PN';

% 获取文件和目录
SourceDirName = C12Dir_2;
[SaveDirName,FileTable] = GetTXTFile(SourceDirName);
Nfile = size(FileTable,1);

% 参与计算的卫星数目
SatNum = 5;

% 删除旧的数据
SaveDirName = [SaveDirName '_' num2str(SatNum)];
delete([SaveDirName '/*.mat'])
%---------------------------------------------------------
% Read GroundTruth from file
if strcmp(SourceDirName,W203Dir)
    NumGroundTruth = 35;
elseif strcmp(SourceDirName,C12Dir_1) || strcmp(SourceDirName,C12Dir_2) 
    NumGroundTruth = 9;
end
fileID = fopen([SaveDirName '/' 'groundTruth.txt'],'r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 NumGroundTruth])';
GroundTruthLLA = FileGroundTruthLLA(:,2:4);

%---------------------------------------------------------
% 循环读取
ilastPoint = 0;
for iFile = 1:Nfile
    prFileName = FileTable(iFile, :);
    % 解析文件名
    a = strsplit(prFileName,'_');
    iPoint = str2double(a{1}(isstrprop(a{1},'digit')));

    % GroundTruth
    param.llaTrueDegDegM = GroundTruthLLA(iPoint,:);
        
    % Filter
    dataFilter = SetDataFilter;
    [gnssRaw,gnssAnalysis] = ReadGnssLogger(SourceDirName,prFileName,dataFilter);
    if isempty(gnssRaw), continue, end
    % ephemeris downloader
    fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
    utcTime = Gps2Utc([],fctSeconds);
    allGpsEph = GetNasaHourlyEphemeris(utcTime,SourceDirName);
    if isempty(allGpsEph), continue, end

    [gnssMeas] = ProcessGnssMeas(gnssRaw);
    % Satlite filter
    M = length(gnssMeas.Svid);
    if M < SatNum
       warning(['Not enough Satlite: ' num2str(M) '<' num2str(SatNum)])
       continue;
    else
       gnssMeas = SatFilter(gnssMeas,SatNum,0);
    end  

    % Data Seperate
    [gnssMeas_BKS, gnssMeas_NBKS] = Seprate(gnssMeas,prFileName,0);

 
    % Original
    gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph); 
    % Distance
    distance.org = distanceM(gpsPvt.allLlaDegDegM,param.llaTrueDegDegM);
    hdop.org = gpsPvt.hdop;
    if isempty(distance.org),  continue, end    
    
% if 0       
    % BKS
    gpsPvt_BKS = GpsWlsPvt(gnssMeas_BKS,allGpsEph); 
    % Distance
    distance.bks = distanceM(gpsPvt_BKS.allLlaDegDegM,param.llaTrueDegDegM);
    hdop.bks = gpsPvt_BKS.hdop;
    
    % NBKS
    gpsPvt_NBKS = GpsWlsPvt(gnssMeas_NBKS,allGpsEph);
    % Distance
    distance.nbks = distanceM(gpsPvt_NBKS.allLlaDegDegM,param.llaTrueDegDegM);
    hdop.nbks = gpsPvt_NBKS.hdop;
% end
    % Hbrid
    gpsPvt_H = GPSWlsPvtBks(gnssMeas,gnssMeas_BKS,gnssMeas_NBKS,allGpsEph,param.llaTrueDegDegM);
    % Distance
    distance.h = distanceM(gpsPvt_H.allLlaDegDegM,param.llaTrueDegDegM);
    hdop.h = gpsPvt_H.hdop;
    % 位置打点
%     h5 = figure;
%     ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
%     PlotPvtBackscatter(gpsPvt_H,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
%     continue    
    
%---------------------------------------------------------
% save file
    if ilastPoint == iPoint
        id = id + 1;
    else
        id = 1;
    end
    ilastPoint = iPoint;
    fname = sprintf('%02d_%02d',iPoint,id);
    save([SaveDirName '/' fname],'distance','hdop');
end