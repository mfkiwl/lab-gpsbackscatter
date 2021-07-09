function [xHat,z,svPos,H,Wpr,Wrr] = WlsPvtBackscatter(prs_BKS,prs_NBKS,gpsEph_BKS, gpsEph_NBKS,xo)
% [xHat,z,svPos,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo)
% 计算WLS PVT 解, xHat
% given pseudoranges, pr rates, and initial state
%
% Inputs: 
%  prs: matrix of raw pseudoranges, and pr rates, each row of the form:
%  [trxWeek,trxSeconds,sv,prMeters,prSigmaMeters,prrMps,prrSigmaMps] 
%   trxWeek, trxSeconds: Rx time of measurement 
%      where trxSeconds = seconds in the current GPS week
%   sv: satellite id number
%   prMeters, prSigmaMeters: pseudorange and standard deviation (meters)
%   prrMps, prrSigmaMps: pseudorange rate and standard deviation (m/s)
%   gpsEph: matching vector of GPS ephemeris struct, defined in ReadRinexNav
%   xo: initial (previous) state, [x,y,z,bc,xDot,yDot,xDot,bcDot]'
%       in ECEF coordinates(meters and m/s)
%
% Outputs: xHat: estimate of state update
%          z = [zPr; zPrr] a-posteriori residuals (measured-calculated)
%          svPos: matrix of calculated sv positions and sv clock error: 
%                 [sv prn, x,y,z (ecef m), dtsv (s),xDot,yDot,zDot, dtsvDot]
%          H: H observation matrix corresponding to svs in svPos(:,1)
%          Wpr,Wrr Weights used in WlsPvt = 1/diag(sigma measurements)
%                  use these matrices to compute variances of xHat
%
% The PVT solution = xo + xHat, in ECEF coordinates
% For unweighted solution, set all sigmas = 1

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

% Check Input
[bOk,numVal] = checkInputs(prs_BKS, gpsEph_BKS, xo);
if ~bOk
    error('prs_BKS输入不对，inputs not right size, or not properly aligned with each other')
end
[bOk,numVal] = checkInputs(prs_NBKS, gpsEph_NBKS, xo);
if ~bOk
    error('prs_NBKS输入不对，inputs not right size, or not properly aligned with each other')
end

% TagTrueDegDegM = [30.96827 118.74069 150];%设置中，背向散射节点的位置
xHat=[]; z=[]; H=[]; svPos=[];
xyz0 = xo(1:3); %初始坐标位置
bc = xo(4);% 初始钟差

if numVal<4  %满足卫星数目大于等于4这一条件才可以，进行求解
  return
end
% ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
% %goes negative, and it is handled in GpsEph2Pvt, where we work with fct
% ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock 
ttxWeek_BKS = prs_BKS(:,jWk);
ttxWeek_NBKS = prs_NBKS(:,jWk);

ttxSeconds_BKS =  prs_BKS(:,jSec) - prs_BKS(:,jPr)/GpsConstants.LIGHTSPEED
ttxSeconds_NBKS =  prs_NBKS(:,jSec) - prs_NBKS(:,jPr)/GpsConstants.LIGHTSPEED
% delta=ttxSeconds_BKS-ttxSeconds_NBKS

% this is accurate satellite time of tx, because we use actual pseudo-ranges 
% here, not corrected ranges
% write the equation for pseudorange to see the rx clock error exactly cancel
% to get precise GPS time: we subtract the satellite clock error from sv time, 
% as done next:
%% 双份星历
dtsv_BKS = GpsEph2Dtsv(gpsEph_BKS,ttxSeconds_BKS);
dtsv_BKS = dtsv_BKS(:); %make into a column for compatibility with other time vectors
ttx_BKS = ttxSeconds_BKS - dtsv_BKS; %subtract dtsv from sv time to get true gps time
%calculate satellite position at ttx
[svXyzTtx_BKS,dtsv_BKS,svXyzDot_BKS,dtsvDot_BKS]=GpsEph2Pvt(gpsEph_BKS,[ttxWeek_BKS,ttx_BKS]);
svXyzTrx_BKS = svXyzTtx_BKS; %initialize svXyz at time of reception

dtsv_NBKS = GpsEph2Dtsv(gpsEph_NBKS,ttxSeconds_NBKS);
dtsv_NBKS = dtsv_NBKS(:); %make into a column for compatibility with other time vectors
ttx_NBKS = ttxSeconds_NBKS - dtsv_NBKS; %subtract dtsv from sv time to get true gps time
%calculate satellite position at ttx
[svXyzTtx_NBKS,dtsv_NBKS,svXyzDot_NBKS,dtsvDot_NBKS]=GpsEph2Pvt(gpsEph_NBKS,[ttxWeek_NBKS,ttx_NBKS]);
svXyzTrx_NBKS = svXyzTtx_NBKS; %initialize svXyz at time of reception
% delat_ttx =ttx_BKS-ttx_NBKS
% delat_sv=svXyzTrx_BKS-svXyzTrx_NBKS

%% 卫星对称映射
sv_length=length(svXyzTrx_BKS(:,1));
svXyzTtx_BKS_mirrored=ones(size(svXyzTrx_BKS));
for i=1: sv_length
    
Temp=mirrorTransform(svXyzTrx_BKS(i,:),xo(1:3).');
% disp(['坐标点：',num2str(xo(1:3)),'卫星坐标：',num2str(svXyzTrx_NBKS(i,:)),'映射坐标：',num2str(Temp)]);
% Coor1=xo(1:3).'
% CoorSv1=svXyzTrx_NBKS(i,:)
% CoorSv2=Temp
% Prm1=norm(CoorSv1-Coor1);
% Prm2=norm(CoorSv2-Coor1);
% delta_PrM=Prm1-Prm2
svXyzTtx_BKS_mirrored(i,:)=Temp;
end

%% % % 这里写下伪距差分的公式，上面计算出来的是真实的卫星位置，我们输入groundtruth计算伪距的差分值
% % %    
% % groundTruth= [30.96827 118.74069 150];
% % GroudtruthXyz=Lla2Xyz(groundTruth);
% % prsTrue=zeros(1,numVal);
% % for i=1:numVal
% %     prsTrue(i)=norm(svXyzTtx(i,:)-GroudtruthXyz);
% % end
% % % 计算获得了真实伪距，与实际伪距做差进行分析
% % prsCorr=prs(:,jPr)-prsTrue(i);
% % prs(:,jPr)=prsTrue;
%% 
%重新构造求解参数矩阵
 prs=[prs_BKS;prs_NBKS];
%  svXyzTrx=[svXyzTrx_BKS_mirrored;svXyzTrx_NBKS];

%%
%%Compute weights ---------------------------------------------------
Wpr = diag(1./prs(:,jPrSig));
Wrr = diag(1./prs(:,jPrrSig));

%iterate on this next part tilL change in pos & line of sight vectors converge
xHat=zeros(4,1); %xHat修改
dx=xHat+inf;
whileCount=0; maxWhileCount=100; 
%we expect the while loop to converge in < 10 iterations, even with initial
%position on other side of the Earth (see Stanford course AA272C "Intro to GPS")
while norm(dx) > GnssThresholds.MAXDELPOSFORNAVM
    whileCount=whileCount+1;
    assert(whileCount < maxWhileCount,...
        'while loop did not converge after %d iterations',whileCount);
    %% 这一步用来做飞行时间矫正，判断飞行时间是否回过于长，弹过和不弹的分开判断
     prs=[prs_BKS;prs_NBKS];
    for i=1:length(gpsEph_BKS)
        % calculate tflight from, bc and dtsv
        dtflight = (prs_BKS(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv_BKS(i);
        % Use of bc: bc>0 <=> pr too big <=> tflight too big.
        %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
        % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
        %   i.e ttx = ttxsv - dtsv
           RawSvXyzTrx_BKS_mirrored=svXyzTtx_BKS_mirrored;
        svXyzTrx_BKS_mirrored(i,:) = FlightTimeCorrection(svXyzTtx_BKS_mirrored(i,:), dtflight);
%         disp(['BKS_dtflight= ',num2str(dtflight)])%输出矫正量
    end
    %% 输出校准量
     CalibrationRes=RawSvXyzTrx_BKS_mirrored-svXyzTrx_BKS_mirrored;
    for i=1:length(gpsEph_NBKS)
        % calculate tflight from, bc and dtsv
        dtflight = (prs_NBKS(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv_NBKS(i);
%         svXyzTrx_NBKS = svXyzTtx_NBKS
           RawSvXyzTrx=svXyzTrx_NBKS;
        svXyzTrx_NBKS(i,:) = FlightTimeCorrection(svXyzTtx_NBKS(i,:), dtflight);
%          disp(['NBKS_dtflight= ',num2str(dtflight)]);%输出校准量
    end
    %% 输出校准量取消误差
    CalibrationRes=RawSvXyzTrx-svXyzTrx_NBKS;
  %calculate line of sight vectors and ranges from satellite to xo
  
  svXyzTrx=[svXyzTrx_BKS_mirrored;svXyzTrx_NBKS];
  numVal=length(svXyzTrx);
  v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
  range = sqrt( sum(v.^2) );
  v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo
dtsv=[dtsv_BKS;dtsv_NBKS];
  svPos=[prs(:,3),svXyzTrx,dtsv(:)];

  %calculate the a-priori range residual
  prHat = range(:) + bc -GpsConstants.LIGHTSPEED*dtsv;
  % Use of bc: bc>0 <=> pr too big <=> rangehat too big
  % Use of dtsv: dtsv>0 <=> pr too small
    
  zPr = prs(:,jPr)-prHat; 
  H = [v', ones(numVal,1)]; % H matrix = [unit vector,1]
  
  %z = Hx, premultiply by W: Wz = WHx, and solve for x:
  dx = pinv(Wpr*H)*Wpr*zPr;

  % update xo, xhat and bc
  xHat=xHat+dx;
  xyz0=xyz0(:)+dx(1:3);
  bc=bc+dx(4);

  %Now calculate the a-posteriori range residual
  zPr = zPr-H*dx;
  
  % display info
  string=['whileCount = ',num2str(whileCount),', norm(dx) =',num2str(norm(dx))];
  disp(string)
  lla0=Xyz2Lla(xyz0(:)');
  stringPos=[' Position: ',num2str(lla0(1)),',',num2str(lla0(2)),',',num2str(lla0(3))];
  disp(stringPos)
%   param.llaTrueDegDegM = [30.511739 114.406770 50];
end
%% 速度部分计算
% % Compute velocities ---------------------------------------------------------
% svXyzDot_BKS,dtsvDot_BKS
rrMps = zeros(numVal,1);
svXyzDot=[svXyzDot_BKS;svXyzDot_NBKS];
dtsvDot=[dtsvDot_BKS;dtsvDot_NBKS];
for i=1:numVal
    %range rate = [satellite velocity] dot product [los from xo to sv]
    rrMps(i) = -svXyzDot(i,:)*v(:,i);
end
prrHat = rrMps + xo(8) - GpsConstants.LIGHTSPEED*dtsvDot;
zPrr = prs(:,jPrr)-prrHat;
%z = Hx, premultiply by W: Wz = WHx, and solve for x:
vHat = pinv(Wrr*H)*Wrr*zPrr;
xHat = [xHat;vHat]; 

z = [zPr;zPrr];

end %end of function WlsPvt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bOk,numVal] = checkInputs(prs, gpsEph, xo)
%utility function for WlsPvt
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

bOk=false;
%check inputs
numVal=size(prs,1);  
% if (max(prs(:,jSec))-min(prs(:,jSec)))> eps %判断时间是否正确 (即所有卫星都来自同一时刻)
if (max(prs(:,jSec))-min(prs(:,jSec)))> 1 %判断时间是否正确 (可以来自相邻时刻)
  return
elseif length(gpsEph)~=numVal
    return
% elseif any(prs(:,jSv) ~= [gpsEph.PRN]') %先注释掉，因为加了七十就不是正常卫星了
%     return
elseif  any(size(xo) ~= [8,1]) % 待求解矩阵xo
    return
elseif size(prs,2)~=7
    return
else
    bOk = true;
end

%We insist that gpsEph and prs are aligned first.
%ClosestGpsEph.m does this, and passes back indices for prs - this is the way to
%do it right, so we don't have nested searches for svId

end %end of function checkInputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

