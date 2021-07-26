%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you: 
% prFileName = 'gnss_log_2020_08_07_15_09_27.txt'; %with duty cycling, no carrier phase
close all;
clear all;
%% data
% % prFileName = 'gnss_log_2019_08_15_10_46_22.txt';
% prFileName = 'P2_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt';
%  
% % Raw and GPS data Path
% % dirName = 'D:\browser\GnssAnalysisWindowsV3.0.3.0\demofiles\adrDemo';
% dirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0612测试集_3Tag_35Point\Tag2_Loc11';
dirName = 'E:\Users\ASUS\Documents\SynologyDrive\SynologyDrive\GPSBackscatter\Data\0720测试\TagLoc1';
% dirName = 'E:\Users\ASUS\Documents\SynologyDrive\SynologyDrive\GPSBackscatter\Data\0719测试\TagLoc1'
% prFileName = 'P02_150mV_100mV_Tag2_gnss_log_2021_06_12_17_44_16.txt'; 
prFileName = 'P12_146mV_100mV_gnss_log_2021_07_20_18_41_26.txt'; 
%% parameters
param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
% param.llaTrueDegDegM = [22.5760482040084	113.936366904008 0];
% param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
% Get the last second, in millisecond
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
% Convert GPS time (week & seconds), or Full Cycle Time (seconds) to UTC
% e.g: 1280778419432
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
% Analyse the raw data (Only changed the name of func)
% [gnssMeas] = ProcessGnssMeas(gnssRaw);
%  [gnssMeas]= ProcessGnssMeasForBackscatter(gnssRaw); % Process raw measurements read from ReadGnssLogger
 
[gnssMeas]= ProcessGnssMeasForBackscatter(gnssRaw);
% if ~exist([dirName '\' 'gnssMean.mat'], 'file')
%     [gnssMeas]= ProcessGnssMeasForBackscatter(gnssRaw);
%     [gnssMeas_BKS, gnssMeas_NBKS]=Seprate(gnssRaw,gnssMeas,prFileName);
%     save( [dirName '\' 'gnssMean.mat'],'gnssMeas', 'gnssMeas_BKS', 'gnssMeas_NBKS');
% else
%     load([dirName '\' 'gnssMean.mat'],'gnssMeas', 'gnssMeas_BKS', 'gnssMeas_NBKS');
% end
 
% [gnssMeasBackscattered]=gnssMeasSeprate(gnssMeas,prFileName);
%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);
%% seperate BKS
[gnssMeas_BKS, gnssMeas_NBKS]=Seprate(gnssMeas,prFileName);

%% Calculate the position of satellites
%% compute WLS position and velocity
% gpsPvt = GpsWlsPvt(gnssMeas_BKS,allGpsEph);
% [gnssMeas] = PrMAdrSmoother(gnssMeas);
% gnssMeas.PrM = gnssMeas.PrMSmooth;
gpsPvt= GpsWlsPvtBackscatter(gnssMeas,allGpsEph);%

%% plot Pvt results
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;

h5 = figure;
PlotPvtStates(gpsPvt,prFileName);

%% Plot Accumulated Delta Range 
if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
    [gnssMeas]= ProcessAdr(gnssMeas);
    h6 = figure;
    PlotAdr(gnssMeas,prFileName,colors);
    [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
    h7 = figure;
    PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
end

% [gnssMeas] = PrMAdrSmoother(gnssMeas);
%% Calculate err based on GroundTruth
% GroundTruth = [22.576062252004220   113.9363809520042;];
% % distance=distanceCountViaCoordinate(gpsPvt.allLlaDegDegM(:,1:2),GroundTruth)
% distance = zeros(length(gpsPvt.allLlaDegDegM(:,1)),1);
% for i = 1:length(gpsPvt.allLlaDegDegM(:,1))
%     distance(i) = distanceCountViaCoordinate(gpsPvt.allLlaDegDegM(i,1:2),GroundTruth);
% end
% figure
% plot(distance)








