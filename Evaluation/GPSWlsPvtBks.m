function [gpsPvt]=GPSWlsPvtBks(gnssMeas,gnssMeas_NBKS, gnssMeas_BKS,allGpsEph,TagllaDegDegM)

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

%%
weekNum = floor(gnssMeas_BKS.FctSeconds/GpsConstants.WEEKSEC);
for i=1:N
    %找弹的数据的第一组
    iValid = find(isfinite(gnssMeas_BKS.PrM(i,:))); %index into valid svid 
    svid    = gnssMeas_BKS.Svid(iValid)';
    [gpsEph_BKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_BKS.FctSeconds(i)); 
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs_BKS = length(svid); %number of satellites this epoch
        
    prM     = gnssMeas_BKS.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas_BKS.PrSigmaM(i,iValid(iSv))';    
    prrMps  = gnssMeas_BKS.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas_BKS.PrrSigmaMps(i,iValid(iSv))';    
    tRx = [ones(numSvs_BKS,1)*weekNum(i),gnssMeas_BKS.tRxSeconds(i,iValid(iSv))'];    
    prs_BKS = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

    %找没弹过的数据的第一组
    iValid = find(isfinite(gnssMeas_NBKS.PrM(i,:))); %index into valid svid 
    svid    = gnssMeas_NBKS.Svid(iValid)';
    [gpsEph_NBKS,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas_NBKS.FctSeconds(i)); 
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
    xo(1:3)= Lla2Xyz(TagllaDegDegM)';
    % [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo);%compute WLS solution
    [xHat,~,~,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph_BKS,gpsEph_NBKS,xo);
    % xHat是10*1的矩阵，其中123位对应定位位置, 4 弹过的钟差, 5 没弹的钟差, 678 速度, 9 钟差弹过的,10 钟差没弹的
    xo = xo + xHat; %包含速度，钟差等参数

    %extract position states
    llaDegDegM = Xyz2Lla(xo(1:3)');
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
    gpsPvt.hdop(i) = sqrt(P(1,1)+P(2,2));
    gpsPvt.pdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3));
    gpsPvt.tdop(i)  = sqrt(P(4,4));
    gpsPvt.gdop(i)  = sqrt(P(1,1)+P(2,2)+P(3,3)+P(4,4));
    
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
end