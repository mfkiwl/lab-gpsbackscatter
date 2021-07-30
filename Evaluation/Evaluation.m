close all;
clear 
SourceDirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';

[FileTable] = GetFile(SourceDirName);

Nfile = size(FileTable,1);

% 删除旧的数据
SaveDirName = './W203_5';
delete([SaveDirName '/*'])
%---------------------------------------------------------
% Read GroundTruth from file
fileID = fopen('groundTruth.txt','r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 35])';
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
    SatNum = 5;
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