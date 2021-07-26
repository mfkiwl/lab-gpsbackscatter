close all;
clear 

dirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
% dirName = 'E:\Users\ASUS\Documents\SynologyDrive\SynologyDrive\GPSBackscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
% dirName = 'E:\Users\ASUS\Documents\SynologyDrive\SynologyDrive\GPSBackscatter\Data\0719测试'
prFileName = 'P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51.txt'; 
% prFileName = 'P04_150mV_100mV_Tag2_gnss_log_2021_06_12_18_35_51.txt'; 
% P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16
% P03_150mV_100mV_Tag2_gnss_log_2021_06_12_18_01_51
% P05_150mV_100mV_Tag2_gnss_log_2021_06_12_18_46_52
%% Read GroundTruth from file
fileID = fopen('groundTruth.txt','r');
formatSpec = '%d %f %f %f';
FileGroundTruthLLA = fscanf(fileID,formatSpec,[4 35])';
GroundTruthLLA = FileGroundTruthLLA(:,2:4);

% 解析文件名
a = strsplit(prFileName,'_');
iPoint = str2double(a{1}(isstrprop(a{1},'digit')));

% 设置GroundTruth
param.llaTrueDegDegM = GroundTruthLLA(iPoint,:);

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
% Satlite filter
% gnssMeas = SatFilter(gnssMeas,4);
% %% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);
%% Data Seperate
% [gnssMeas_BKS, gnssMeas_NBKS] = Seprate(gnssMeas,prFileName);
[gnssMeas_NBKS, gnssMeas_BKS] = Seprate(gnssMeas,prFileName,0);%位置反转

if 1
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
% 根据当前时刻的卫星位置判断出未来卫星位置
N1 = length(gnssMeas_BKS.FctSeconds);
N2 = length(gnssMeas_NBKS.FctSeconds);
% 选取更短的N
N = min([N1 N2]);
%% gpsPvt 初始化赋值
gpsPvt.FctSeconds      = gnssMeas.FctSeconds;
gpsPvt.allXyzM   = zeros(N,3)+NaN;
gpsPvt.allLlaDegDegM   = zeros(N,3)+NaN; 
gpsPvt.sigmaLLaM       = zeros(N,3)+NaN;
% gpsPvt.allBcMeters     = zeros(N,1)+NaN;
gpsPvt.allBcMeters_BKS     = zeros(N,1)+NaN;
gpsPvt.allBcMeters_NBKS     = zeros(N,1)+NaN;

gpsPvt.allVelMps       = zeros(N,3)+NaN;
gpsPvt.sigmaVelMps     = zeros(N,3)+NaN;
% gpsPvt.allBcDotMps     = zeros(N,1)+NaN;
gpsPvt.allBcDotMps_BKS     = zeros(N,1)+NaN;
gpsPvt.allBcDotMps_NBKS     = zeros(N,1)+NaN;

gpsPvt.numSvs          = zeros(N,1);
gpsPvt.hdop            = zeros(N,1)+inf;
gpsPvt.pdop            = zeros(N,1)+inf;
gpsPvt.tdop            = zeros(N,1)+inf;
gpsPvt.gdop            = zeros(N,1)+inf;

weekNum = floor(gnssMeas_BKS.FctSeconds/GpsConstants.WEEKSEC);
for i= 1:N
    %找弹的数据的第一组
    iValid = find(isfinite(gnssMeas_BKS.PrM(i,:))); %index into valid svid %这一步没看懂
    svid    = gnssMeas_BKS.Svid(iValid)';
    [gpsEph_BKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_BKS.FctSeconds(i)); %从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs_BKS = length(svid); %number of satellites this epoch
        
    prM     = gnssMeas_BKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_BKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_BKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_BKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs_BKS,1)*weekNum(i),gnssMeas_BKS.tRxSeconds(i,iValid(iSv))'];    
    prs_BKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    %找没弹过的数据的第一组
    iValid = find(isfinite(gnssMeas_NBKS.PrM(i,:))); %index into valid svid % 这一步没看懂
    svid    = gnssMeas_NBKS.Svid(iValid)';
    [gpsEph_NBKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_NBKS.FctSeconds(i)); % 从星历中挑选对应的卫星
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs_NBKS = length(svid); %number of satellites this epoch
    
    prM     = gnssMeas_NBKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_NBKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_NBKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_NBKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs_NBKS,1)*weekNum(i),gnssMeas_NBKS.tRxSeconds(i,iValid(iSv))'];    
    prs_NBKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    % WLS
    xo =zeros(10,1);
    xo(6:8) = zeros(3,1); %initialize speed to zero
    xo(1:3)= Lla2Xyz( GroundTruthLLA(6,:))';
    % [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo);%compute WLS solution
    [xHat,~,~,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph_BKS,gpsEph_NBKS,xo);
%     [xHat,~,~,H,Wpr,Wrr] = WlsPvtCoinGPS(prs_BKS,prs_NBKS,gpsEph_BKS,gpsEph_NBKS,xo);
    % xHat是10*1的矩阵，其中123位对应定位位置, 4 弹过的钟差, 5 没弹的钟差, 678 速度, 9 钟差弹过的,10 钟差没弹的
    xo = xo + xHat; %包含速度，钟差等参数

    %extract position states
    llaDegDegM = Xyz2Lla(xo(1:3)');
%     gpsPvt_H.allLlaDegDegM(i,:) = llaDegDegM;
    gpsPvt.allLlaDegDegM(i,:) = llaDegDegM;
    gpsPvt.allXyzM(i,:)   = xo(1:3)';
    gpsPvt.allBcMeters_BKS(i)     = xo(4);
    gpsPvt.allBcMeters_NBKS(i)    = xo(5);
    
    %extract velocity states
    RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
    %NOTE: in real-time code compute RE2N once until position changes
    vNed = RE2N*xo(6:8); %velocity in NED 速度分量改变了
    gpsPvt.allVelMps(i,:) = vNed;
%     gpsPvt.allBcDotMps(i) = xo(9);
    gpsPvt.allBcDotMps_BKS(i) = xo(9);
    gpsPvt.allBcDotMps_NBKS(i)= xo(10);
    
    %compute HDOP
    H = [H(:,1:3)*RE2N', ones(numSvs_BKS + numSvs_NBKS,1)]; %observation matrix in NED
    P = inv(H'*H);%unweighted covariance
%     if P(1,1)+P(2,2)<1.01e-5
% %         A=1;
%         disp(['病态方程',num2str(i)])
%     end
    gpsPvt.hdop(i) = sqrt(P(1,1)+P(2,2));
    gpsPvt.pdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3));
    gpsPvt.tdop(i)  = sqrt(P(4,4));
    gpsPvt.gdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3)+P(4,4));
    
%     disp(['GPS HDOP：',num2str(gpsPvt.hdop(i))])
%     disp(['GPS PDOP：',num2str(gpsPvt.pdop(i))])
%     disp(['GPS TDOP：',num2str(gpsPvt.tdop(i))])
%     disp(['GPS GDOP：',num2str(gpsPvt.gdop(i))])
    
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

%% 位置打点
% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) < 60),:);
% gpsPvt_H.allLlaDegDegM = gpsPvt_H.allLlaDegDegM(find(gpsPvt_H.allLlaDegDegM(:,3) > 30),:);
h5 = figure;
ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvtBackscatter(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
%% bc1-bc2
% allBcMeters   Nx1 common bias computed with llaDegDegM
% allBcDotMps   Nx1 common freq bias computed with velocity
Bc_delta= gpsPvt.allBcMeters_BKS-gpsPvt.allBcMeters_NBKS;
BcV_delta=gpsPvt.allBcDotMps_BKS-gpsPvt.allBcDotMps_NBKS;

% 
gpsPvt.FctSeconds      = 1:length(gpsPvt.numSvs);%FctSeconds缩短了
bc_Rate = polyfit(gpsPvt.FctSeconds,Bc_delta,1)./2;
ToF_deltaMeters=Bc_delta./(-2.*bc_Rate(:,1));
ToF_deltaSeconds=ToF_deltaMeters./GpsConstants.LIGHTSPEED;

% 3D distance
groundTruth=Lla2Xyz(param.llaTrueDegDegM);
distanceVector=gpsPvt.allXyzM-ones(N,1)*groundTruth;
distance=ones(N,1);
for i=1:N
    distance(i)=norm(distanceVector(i,:));
end

figure;
hold on
plot(Bc_delta);
plot(BcV_delta);
title('Clock Bias delta * LightSpeed in Meters'),ylabel('(Meters)')
xs = sprintf('time (seconds)\n%s',prFileName);
xlabel(xs,'Interpreter','none')

figure;
hold on
plot(gpsPvt.allBcMeters_BKS/GpsConstants.LIGHTSPEED);
plot(gpsPvt.allBcMeters_NBKS/GpsConstants.LIGHTSPEED);
title('Clock Bias '),ylabel('(Meters)')
xs = sprintf('time (seconds)\n%s',prFileName);
xlabel(xs,'Interpreter','none')
% title(prFileName);
hold off
%% 分析钟差
gpsPvt.allBcMeters_BKS = filloutliers(gpsPvt.allBcMeters_BKS,'linear');
gpsPvt.allBcMeters_NBKS = filloutliers(gpsPvt.allBcMeters_NBKS,'linear');
gpsPvt.allBcMeters_BKS =smooth(gpsPvt.allBcMeters_BKS);
gpsPvt.allBcMeters_NBKS=smooth(gpsPvt.allBcMeters_NBKS);

gpsPvt.allBcMeters_BKS = filloutliers(gpsPvt.allBcMeters_BKS,'linear');
gpsPvt.allBcMeters_NBKS = filloutliers(gpsPvt.allBcMeters_NBKS,'linear');
gpsPvt.allBcMeters_BKS =smooth(gpsPvt.allBcMeters_BKS);
gpsPvt.allBcMeters_NBKS=smooth(gpsPvt.allBcMeters_NBKS); 

figure;
hold on
plot(gpsPvt.allBcMeters_BKS/GpsConstants.LIGHTSPEED);
plot(gpsPvt.allBcMeters_NBKS/GpsConstants.LIGHTSPEED);
title('Clock Bias (Smoothed)'),ylabel('(Meters)')
xs = sprintf('time (seconds)\n%s',prFileName);
xlabel(xs,'Interpreter','none')
% title(prFileName);
hold off

 R_BKS=polyfit(gpsPvt.FctSeconds,gpsPvt.allBcMeters_BKS/GpsConstants.LIGHTSPEED,1);
 K_BKS=R_BKS(1);
 B_BKS=R_BKS(2);
 P_BKS=polyval(R_BKS,gpsPvt.FctSeconds);
 
 R_NBKS=polyfit(gpsPvt.FctSeconds,gpsPvt.allBcMeters_NBKS/GpsConstants.LIGHTSPEED,1);
 K_NBKS=R_NBKS(1);
 B_NBKS=R_NBKS(2);
 P_NBKS=polyval(R_NBKS,gpsPvt.FctSeconds);
 
 figure;
 hold on;
 plot(gpsPvt.FctSeconds,P_BKS);
 plot(gpsPvt.FctSeconds,P_NBKS);
 hold off
 
 distance=GpsConstants.LIGHTSPEED.*(K_NBKS-K_BKS);
 latitudeDistribution=distance*16/18.86;
 longitudeDistribution=distance*10/18.86;
 
 %% 位置校正

if latitudeDistribution < 20
% % 
caliPosition=ones(size(gpsPvt.allLlaDegDegM)).*[longitudeDistribution.*10e-7 latitudeDistribution.*10e-6 0]
gpsPvt.allLlaDegDegM= gpsPvt.allLlaDegDegM + caliPosition;
else
% % naive方式
caliPosition=ones(size(gpsPvt.allLlaDegDegM)).*[0 0.0002 0];  
gpsPvt.allLlaDegDegM= gpsPvt.allLlaDegDegM + caliPosition;  
end
h8 = figure;
ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
PlotPvtBackscatter(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图


%  CorrectionVector=P(1:3).*distance;
%  CorrectionVectorLla=Xyz2Lla(CorrectionVector);

% figure;
% plot(distance);
% title('3D Positioning error in Meters'),ylabel('(Meters)')
% xs = sprintf('time (seconds)\n%s',prFileName);
% xlabel(xs,'Interpreter','none')

end