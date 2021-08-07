DirName = 'D:\file\Lab-Drive\Project\GPS_Backscatter\Data\0804测试 (2)';
[FileTable,Nfile] = GetAllFile(DirName,'txt')


% 参与计算的卫星数目
SatNum = 5;

SourceDirName = DirName;
llaAll = [nan nan nan];
for iFile = 1:Nfile
    prFileName = FileTable{iFile};
    % 解析文件名
    a = strsplit(prFileName,'_');
    iPoint = str2double(a{2}(isstrprop(a{2},'digit')));

    % GroundTruth
%     param.llaTrueDegDegM = GroundTruthLLA(iPoint,:);
    param.llaTrueDegDegM = [30.5113070328088	114.401749572179	37.8595181058496];
        
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
    
    % Hbrid
    gpsPvt_H = GPSWlsPvtBks(gnssMeas,gnssMeas_BKS,gnssMeas_NBKS,allGpsEph,param.llaTrueDegDegM);

    % 位置打点
%     hh = figure;
%     ts = 'HBKS_Raw Pseudoranges, Weighted Least Squares solution';
%     PlotPvtBackscatter(gpsPvt_H,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
    
    llaAll = [llaAll ; gpsPvt_H.allLlaDegDegM];
end

% h = figure;
% ts = 'HBKS Raw Pseudoranges, Weighted Least Squares solution';
% PlotPvtBackscatter(gpsPvt_H,prFileName,param.llaTrueDegDegM,ts); drawnow;%绘制位置图
llaAll = llaAll(~ismissing(llaAll(:,1)),:)

figure
plot(llaAll(:,2),llaAll(:,1),'+','MarkerSize',18);
grid on
grid minor

llaAver = mean(llaAll,1);%median position
nedM = Lla2Ned(llaAll,llaAver);
figure
h = plot(nedM(:,2),nedM(:,1),'+','MarkerSize',18);
grid on
grid minor