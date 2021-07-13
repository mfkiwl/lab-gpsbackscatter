close all;
clear all;
dirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
% dirName = 'D:\Filerec\dingding\泛悦城';

prFileName = 'P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt'; 
% prFileName = 'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_19_10_10.txt'; 
% prFileName = 'gnss_log_2021_07_13_14_06_17.txt'; 
%% Read GroundTruth from file
fileID = fopen('groundTruth.txt','r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 35])';
GroundTruthLLA = FileGroundTruthLLA(:,2:4);
% param.llaTrueDegDegM = [30.511739 114.406770 50]; % 设置GroundTruth
param.llaTrueDegDegM = GroundTruthLLA(2,:);
% param.llaTrueDegDegM = [];
%% Filter
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end
%% ephemeris downloader
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end
%% 
[gnssMeas] = ProcessGnssMeas(gnssRaw);
% %% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);
%% Data Seperate
[gnssMeas_BKS, gnssMeas_NBKS] = Seprate(gnssMeas,prFileName);

% temp = gnssMeas_BKS;
% gnssMeas_BKS = gnssMeas_NBKS;
% gnssMeas_NBKS = temp;

if 0
%% plot Pvt results
% Original
gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph); 
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
% h5 = figure;
% PlotPvtStates(gpsPvt_BKS,prFileName);
% return

% BKS
gpsPvt_BKS = GpsWlsPvt(gnssMeas_BKS,allGpsEph); 
h4 = figure;
ts = 'BKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt_BKS,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
% h5 = figure;
% PlotPvtStates(gpsPvt_BKS,prFileName);

% NBKS
gpsPvt_NBKS = GpsWlsPvt(gnssMeas_NBKS,allGpsEph); 
h4 = figure;
ts = 'NBKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt_NBKS,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
% h5 = figure;
% PlotPvtStates(gpsPvt_NBKS,prFileName);
end

%%
if 1
% 这一部分程序根据当前时刻的卫星位置判断出未来卫星位置
N1 = length(gnssMeas_BKS.FctSeconds);
N2 = length(gnssMeas_NBKS.FctSeconds);
% 选取更短的N
N = min([N1 N2]);

weekNum     = floor(gnssMeas_BKS.FctSeconds/GpsConstants.WEEKSEC);
for i= 1:N
    %找弹的数据的第一组
    iValid = find(isfinite(gnssMeas_BKS.PrM(i,:))); %index into valid svid %这一步没看懂
    svid    = gnssMeas_BKS.Svid(iValid)';
    [gpsEph_BKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_BKS.FctSeconds(i)); %从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid) %number of satellites this epoch
        
    prM     = gnssMeas_BKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_BKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_BKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_BKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs,1)*weekNum(i),gnssMeas_BKS.tRxSeconds(i,iValid(iSv))'];    
    prs_BKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    %找没弹过的数据的第一组
    iValid = find(isfinite(gnssMeas_NBKS.PrM(i,:))); %index into valid svid % 这一步没看懂
    svid    = gnssMeas_NBKS.Svid(iValid)';
    [gpsEph_NBKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_NBKS.FctSeconds(i)); % 从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid) %number of satellites this epoch
    
    prM     = gnssMeas_NBKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_NBKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_NBKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_NBKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs,1)*weekNum(i),gnssMeas_NBKS.tRxSeconds(i,iValid(iSv))'];    
    prs_NBKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    % WLS
    xo =zeros(8,1);
    xo(5:7) = zeros(3,1); %initialize speed to zero
    xo(1:3)= Lla2Xyz( GroundTruthLLA(6,:))';
    % [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo);%compute WLS solution
    [xHat,~,~,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph_BKS,gpsEph_NBKS,xo);
    xo = xo + xHat;
    
    %extract position states
    llaDegDegM = Xyz2Lla(xo(1:3)');
    gpsPvt_H.allLlaDegDegM(i,:) = llaDegDegM;
end

% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) < 60),:);
% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) > 30),:);

h5 = figure;
ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt_H,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
end