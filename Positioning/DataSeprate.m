close all;
clear all;
dirName = 'E:\Users\ASUS\Documents\SynologyDrive\SynologyDrive\GPSBackscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';

% prFileName = 'P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt'; 
% P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51
% prFileName = 'P03_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_51.txt'; 
prFileName = 'P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51.txt'; 
% prFileName = 'P07_150mV_100mV_Tag2_gnss_log_2021_06_12_19_10_10.txt'; 
%% Read GroundTruth from file
fileID = fopen('groundTruth.txt','r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 35])';
GroundTruthLLA = FileGroundTruthLLA(:,2:4);
% param.llaTrueDegDegM = [30.511739 114.406770 50]; % 设置GroundTruth
param.llaTrueDegDegM = GroundTruthLLA(7,:);
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

if 0
%% plot Pvt results
% Original
gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph); 
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
% h5 = figure;
% PlotPvtStates(gpsPvt_BKS,prFileName);

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

%% 初始化赋值
gpsPvt.FctSeconds      = gnssMeas.FctSeconds;
gpsPvt.allLlaDegDegM   = zeros(N,3)+NaN; 
gpsPvt.sigmaLLaM       = zeros(N,3)+NaN;
gpsPvt.allBcMeters     = zeros(N,1)+NaN;
gpsPvt.allVelMps       = zeros(N,3)+NaN;
gpsPvt.sigmaVelMps     = zeros(N,3)+NaN;
gpsPvt.allBcDotMps     = zeros(N,1)+NaN;
gpsPvt.numSvs          = zeros(N,1);
gpsPvt.hdop            = zeros(N,1)+inf;
gpsPvt.pdop            = zeros(N,1)+inf;
gpsPvt.tdop            = zeros(N,1)+inf;
gpsPvt.gdop            = zeros(N,1)+inf;

for i= 1:N
    %找弹的数据的第一组
    iValid = find(isfinite(gnssMeas_BKS.PrM(i,:))); %index into valid svid %这一步没看懂
    svid    = gnssMeas_BKS.Svid(iValid)';
    [gpsEph_BKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_BKS.FctSeconds(i)); %从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
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
    [gpsEph_NBKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_NBKS.FctSeconds(i)); % 从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid); %number of satellites this epoch
    
    prM     = gnssMeas_NBKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_NBKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_NBKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_NBKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs,1)*weekNum(i),gnssMeas_NBKS.tRxSeconds(i,iValid(iSv))'];    
    prs_NBKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    % WLS
    xo =zeros(10,1);
    xo(6:8) = zeros(3,1); %initialize speed to zero
    xo(1:3)= Lla2Xyz( GroundTruthLLA(6,:))';
    % [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo);%compute WLS solution
    [xHat,~,~,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph_BKS,gpsEph_NBKS,xo);
    % xHat是10*1的矩阵，其中123位对应定位位置 4，钟差弹过的，5钟差没弹的卫星；678对应速度，9钟差弹过的，10钟差没弹的
    xo = xo + xHat; %包含速度，钟差等参数
%     xo(1:3) = xo(1:3) + xHat(1:3);
%     30.5118,114.4068,43.758
    %extract position states
    llaDegDegM = Xyz2Lla(xo(1:3)');
%     gpsPvt_H.allLlaDegDegM(i,:) = llaDegDegM;
    gpsPvt.allLlaDegDegM(i,:) = llaDegDegM;
    %%
    RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
    %NOTE: in real-time code compute RE2N once until position changes
    vNed = RE2N*xo(5:7); %velocity in NED
    gpsPvt.allVelMps(i,:) = vNed;
    gpsPvt.allBcDotMps(i) = xo(8);
    
    %compute HDOP
    % 
    numSvs1=length(prs_BKS(:,1));
    numSvs2=length(prs_NBKS(:,1));
    H = [H(:,1:3)*RE2N', ones(numSvs1+numSvs2,1)]; %observation matrix in NED
    P = inv(H'*H);%unweighted covariance
    gpsPvt.hdop(i) = sqrt(P(1,1)+P(2,2));
    gpsPvt.pdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3));
    gpsPvt.tdop(i)  = sqrt(P(4,4));
    gpsPvt.gdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3)+P(4,4));
    disp(['GPS HDOP：',num2str(gpsPvt.hdop(i))])
    disp(['GPS PDOP：',num2str(gpsPvt.pdop(i))])
    disp(['GPS TDOP：',num2str(gpsPvt.tdop(i))])
    disp(['GPS GDOP：',num2str(gpsPvt.gdop(i))])
    
    %compute variance of llaDegDegM
    %inside LsPvt the weights are used like this: 
    %  z = Hx, premultiply by W: Wz = WHx, and solve for x:
    %  x = pinv(Wpr*H)*Wpr*zPr;
    %  the point of the weights is to make sigma(Wz) = 1
    %  therefore, the variances of x come from  diag(inv(H'Wpr'WprH))
    P = inv(H'*(Wpr'*Wpr)*H); %weighted covariance
    gpsPvt.sigmaLLaM(i,:) = sqrt(diag(P(1:3,1:3)));
    
    %similarly, compute variance of velocity
    P = inv(H'*(Wrr'*Wrr)*H); %weighted covariance
    gpsPvt.sigmaVelMps(i,:) = sqrt(diag(P(1:3,1:3)));
end

% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) < 60),:);
% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) > 30),:);
gpsPvt.FctSeconds      = 1:length(gpsPvt.numSvs);%FctSeconds缩短了
h5 = figure;
ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
end