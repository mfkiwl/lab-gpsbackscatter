close all;
clear 
SourceDirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';

FileTable = [
'P01_150mV_100mV_Tag2_gnss_log_2021_06_12_17_13_53.txt';
'P01_150mV_100mV_Tag2_gnss_log_2021_06_12_17_32_33.txt';
'P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt';
'P03_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_51.txt';
'P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51.txt';
'P05_150mV_100mV_Tag2_gnss_log_2021_06_12_18_46_52.txt';
'P06_150mV_100mV_Tag2_gnss_log_2021_06_12_17_50_41.txt';
% 'P06_150mv_100mV_Tag2_gnss_log_2021_06_12_18_59_13.txt';
'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_46.txt';
'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_19_10_10.txt';
'P08_150mV_100mV_Tag2_gnss_log_2021_06_12_19_33_19.txt';
'P08_150mV_100mV_Tag2_gnss_log_2021_06_12_19_58_04.txt';
'P09_150mV_100mV_Tag2_gnss_log_2021_06_12_20_09_05.txt';
'P09_150mV_100mV_Tag2_gnss_log_2021_06_12_20_54_12.txt';
'P10_152mV_100mV_Tag2_gnss_log_2021_06_12_21_04_39.txt';
'P12_150mV_100mV_Tag2_gnss_log_2021_06_12_21_15_30.txt';
'P13_150mV_100mV_Tag2_gnss_log_2021_06_12_21_25_55.txt';
'P14_150mV_100mV_Tag2_gnss_log_2021_06_12_21_36_49.txt';
'P15_150mV_100mV_Tag2_gnss_log_2021_06_12_21_48_25.txt';
'P16_150mV_100mV_Tag2_gnss_log_2021_06_12_22_10_28.txt';
'P16_150mV_100mV_Tag2_gnss_log_2021_06_12_22_23_55.txt';
'P17_150mV_100mV_Tag2_gnss_log_2021_06_12_22_34_52.txt';
'P18_150mV_100mV_Tag2_gnss_log_2021_06_12_22_45_33.txt';
'P19_150mV_100mV_Tag2_gnss_log_2021_06_12_22_56_08.txt';
'P20_150mV_100mV_Tag2_gnss_log_2021_06_12_23_06_39.txt';
'P21_150mV_100mV_Tag2_gnss_log_2021_06_13_20_17_14.txt';
'P21_150mV_100mV_Tag2_gnss_log_2021_06_13_20_28_42.txt';
'P22_150mV_100mV_Tag2_gnss_log_2021_06_13_20_39_21.txt';
'P23_150mV_100mV_Tag2_gnss_log_2021_06_13_20_50_35.txt';
'P24_150mV_100mV_Tag2_gnss_log_2021_06_13_21_01_53.txt';
'P25_150mV_100mV_Tag2_gnss_log_2021_06_13_21_13_45.txt';
'P26_150mV_100mV_Tag2_gnss_log_2021_06_13_21_25_00.txt';
'P27_150mV_100mV_Tag2_gnss_log_2021_06_13_21_37_08.txt';
'P28_150mV_100mV_Tag2_gnss_log_2021_06_13_21_49_03.txt';
'P28_150mV_100mV_Tag2_gnss_log_2021_06_14_17_26_41.txt';
'P29_150mV_100mV_Tag2_gnss_log_2021_06_13_22_09_41.txt';
'P29_150mV_100mV_Tag2_gnss_log_2021_06_14_17_15_42.txt';
'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_15_13_24.txt';
'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_15_24_38.txt';
'P30_150mV_100mV_Tag2_gnss_log_2021_06_14_17_04_51.txt';
% 'P31_150mV_100mV_Tag2_gnss_log_2021_06_14_15_42_07.txt'; 
'P32_150mV_100mV_Tag2_gnss_log_2021_06_14_15_53_50.txt';
'P33_150mV_100mV_Tag2_gnss_log_2021_06_14_16_31_26.txt';
'P34_150mV_100mV_Tag2_gnss_log_2021_06_14_16_42_18.txt';
'P35_150mV_100mV_Tag2_gnss_log_2021_06_14_16_54_20.txt';  
];
Nfile = size(FileTable,1);

% 删除旧的数据
SaveDirName = './W203';
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