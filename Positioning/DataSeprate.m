close all;
clear all;
dirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';

prFileName = 'P2_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt'; 

%% 
fileID = fopen('groundTruth.txt','r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 35])';

% param.llaTrueDegDegM = [30.511739 114.406770 50]; % 设置GroundTruth
param.llaTrueDegDegM = FileGroundTruthLLA(:,2:4);
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
%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);
%%
[gnssMeas_BKS, gnssMeas_NBKS] = Seprate(gnssRaw, gnssMeas,prFileName)
%% plot Pvt results
gpsPvt_BKS = GpsWlsPvt(gnssMeas_BKS,allGpsEph); 
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt_BKS,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
h5 = figure;
PlotPvtStates(gpsPvt_BKS,prFileName);

gpsPvt_NBKS = GpsWlsPvt(gnssMeas_NBKS,allGpsEph); 
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt_NBKS,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
h5 = figure;
PlotPvtStates(gpsPvt_NBKS,prFileName);
%% 卫星坐标映射
% 这一部分程序根据当前时刻的卫星位置判断出未来卫星位置
N1 = length(gnssMeas_BKS.FctSeconds);
N2 = length(gnssMeas_NBKS.FctSeconds);
% 选取更短的N
if N1>N2
    N=N2;
else
    N=N1;
end

for i= 1:N
    %找弹的数据的第一组
iValid = find(isfinite(gnssMeas_BKS.PrM(i,:))); %index into valid svid %这一步没看懂
svid    = gnssMeas_BKS.Svid(iValid)';
[gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_BKS.FctSeconds(i)); %从星历中挑选对应的卫星
weekNum     = floor(gnssMeas_BKS.FctSeconds/GpsConstants.WEEKSEC);
numSvs = length(svid); %number of satellites this epoch
prM     = gnssMeas_BKS.PrM(i,iValid(iSv))';
prSigmaM= gnssMeas_BKS.PrSigmaM(i,iValid(iSv))';    
prrMps  = gnssMeas_BKS.PrrMps(i,iValid(iSv))';
prrSigmaMps = gnssMeas_BKS.PrrSigmaMps(i,iValid(iSv))';    
tRx = [ones(numSvs,1)*weekNum(i),gnssMeas_BKS.tRxSeconds(i,iValid(iSv))'];    
prs_BKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

%找没弹过的数据的第一组
iValid = find(isfinite(gnssMeas_NBKS.PrM(i,:))); %index into valid svid % 这一步没看懂
svid    = gnssMeas_NBKS.Svid(iValid)';
[gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_NBKS.FctSeconds(i)); % 从星历中挑选对应的卫星
weekNum     = floor(gnssMeas_NBKS.FctSeconds/GpsConstants.WEEKSEC);
numSvs = length(svid); %number of satellites this epoch
prM     = gnssMeas_NBKS.PrM(i,iValid(iSv))';
prSigmaM= gnssMeas_NBKS.PrSigmaM(i,iValid(iSv))';    
prrMps  = gnssMeas_NBKS.PrrMps(i,iValid(iSv))';
prrSigmaMps = gnssMeas_NBKS.PrrSigmaMps(i,iValid(iSv))';    
tRx = [ones(numSvs,1)*weekNum(i),gnssMeas_NBKS.tRxSeconds(i,iValid(iSv))'];    
prs_NBKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

xo =zeros(8,1);
xo(5:7) = zeros(3,1); %initialize speed to zero
% ground truth 30.969,118.7409,100
% groundTruth= [30.96827 118.74069 150];
xo(1:3)= Lla2Xyz(param.llaTrueDegDegM)';
% [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo);%compute WLS solution
[xHat,~,~,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph,xo);
xo = xo + xHat;
end