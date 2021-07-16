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
[bOk,numVal_BKS] = checkInputs(prs_BKS, gpsEph_BKS, xo);
if ~bOk
    error('prs_BKS输入不对，inputs not right size, or not properly aligned with each other')
end
[bOk,numVal_NBKS] = checkInputs(prs_NBKS, gpsEph_NBKS, xo);
if ~bOk
    error('prs_NBKS输入不对，inputs not right size, or not properly aligned with each other')
end
%   load ToF_deltaSeconds.mat; %%加入ToF校准值
%% 
xHat=[]; z=[]; H=[]; svPos=[];
TagLocXyz0 = xo(1:3); %初始坐标位置 tag 的位置
xyz0=xo(1:3);%定位初始位置

%% 钟差也要视为两个参数进行考虑
% bc = xo(4);% 初始钟差
bc_BKS=xo(4);
bc_NBKS=xo(5);

% if numVal<4  %满足卫星数目大于等于4这一条件才可以，进行求解
%   return
% end

% ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
% %goes negative, and it is handled in GpsEph2Pvt, where we work with fct
% ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock 
ttxWeek_BKS = prs_BKS(:,jWk);
ttxWeek_NBKS = prs_NBKS(:,jWk);

ttxSeconds_BKS =  prs_BKS(:,jSec) - prs_BKS(:,jPr)/GpsConstants.LIGHTSPEED;
ttxSeconds_NBKS =  prs_NBKS(:,jSec) - prs_NBKS(:,jPr)/GpsConstants.LIGHTSPEED;
%  delta=ttxSeconds_BKS-ttxSeconds_NBKS

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
    Temp=mirrorTransform(svXyzTrx_BKS(i,:),TagLocXyz0');
    % disp(['坐标点：',num2str(xo(1:3)),'卫星坐标：',num2str(svXyzTrx_NBKS(i,:)),'映射坐标：',num2str(Temp)]);
    % Coor1=xo(1:3).'
    % CoorSv1=svXyzTrx_NBKS(i,:)
    % CoorSv2=Temp
    % Prm1=norm(CoorSv1-Coor1);
    % Prm2=norm(CoorSv2-Coor1);
    % delta_PrM=Prm1-Prm2
    svXyzTtx_BKS_mirrored(i,:)=Temp;
end
svXyzTrx_BKS_mirrored = svXyzTtx_BKS_mirrored; %initialize svXyz at time of reception
% svXyzTrx_BKS_mirrored = svXyzTrx_BKS; %initialize svXyz at time of reception
dtsv=[dtsv_BKS;dtsv_NBKS];

%% % % 这里写下伪距差分的公式，上面计算出来的是真实的卫星位置，我们输入groundtruth计算伪距的差分值
% % %      1.伪距分析
% TagLocXyz0'% Tag位置
% % svXyzTtx_NBKS % 卫星位置
% % 针对弹过的信号
% TrueDistance_BKS=ones(sv_length,1);
% for i= 1:sv_length
% TrueDistance_BKS(i)= norm(svXyzTtx_BKS(i,:)- TagLocXyz0');
% end
% prsCorr_BKS = prs_BKS(:,jPr)-TrueDistance_BKS;
% figure;
% plot(prsCorr_BKS);
%
%          2. 伪距/位置差分
% 校正
% prs_BKS(:,jPr)=prs_BKS(:,jPr)-prsCorr_BKS;
%% 

%  svXyzTrx=[svXyzTrx_BKS_mirrored;svXyzTrx_NBKS];

%%
prs=[prs_BKS;prs_NBKS];
%%Compute weights ---------------------------------------------------
Wpr = diag(1./prs(:,jPrSig));
Wrr = diag(1./prs(:,jPrrSig));

%iterate on this next part tilL change in pos & line of sight vectors converge
xHat=zeros(5,1); %修改 xHat 为待求解参数,为 x,y,z,bc1,bc2, 其中有两个 bias
dx=xHat+inf;
whileCount=0; maxWhileCount=100; 
%we expect the while loop to converge in < 10 iterations, even with initial
%position on other side of the Earth (see Stanford course AA272C "Intro to GPS")
while norm(dx) > GnssThresholds.MAXPRRUNCMPS % 10 % MAXDELPOSFORNAVM  % 20 m
    whileCount=whileCount+1;
    assert(whileCount < maxWhileCount,...
        'while loop did not converge after %d iterations',whileCount);
    %% 这一步用来做飞行时间矫正，判断飞行时间是否回过于长，弹过和不弹的分开判断
%      prs=[prs_BKS;prs_NBKS];
    for i=1:length(gpsEph_BKS)
        % calculate tflight from, bc and dtsv
        %dtflight也分开计算
        dtflight = (prs_BKS(i,jPr)-bc_BKS)/GpsConstants.LIGHTSPEED + dtsv_BKS(i);
        % Use of bc: bc>0 <=> pr too big <=> tflight too big.
        %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
        % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
        %   i.e ttx = ttxsv - dtsv
%         RawSvXyzTrx_BKS_mirrored = svXyzTtx_BKS_mirrored;
        %对称映射
        svXyzTrx_BKS_mirrored(i,:) = FlightTimeCorrection(svXyzTtx_BKS_mirrored(i,:), dtflight);
%         disp(['BKS_dtflight= ',num2str(dtflight)])%输出矫正量
    end
    %输出校准量
%      CalibrationRes=RawSvXyzTrx_BKS_mirrored-svXyzTrx_BKS_mirrored;
     
    for i=1:length(gpsEph_NBKS)
        % calculate tflight from, bc and dtsv
        dtflight = (prs_NBKS(i,jPr)-bc_NBKS)/GpsConstants.LIGHTSPEED + dtsv_NBKS(i);
%         svXyzTrx_NBKS = svXyzTtx_NBKS
%            RawSvXyzTrx=svXyzTrx_NBKS;
        svXyzTrx_NBKS(i,:) = FlightTimeCorrection(svXyzTtx_NBKS(i,:), dtflight);
%          disp(['NBKS_dtflight= ',num2str(dtflight)]);%输出校准量
    end
    % 输出校准量取消误差
%     CalibrationRes=RawSvXyzTrx-svXyzTrx_NBKS;

 %% 分别计算距离矩阵和距离矢量
  %calculate line of sight vectors and ranges from satellite to xo
%   svXyzTrx=[svXyzTrx_BKS_mirrored;svXyzTrx_NBKS];
%   v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
    v_BKS=xyz0(:)*ones(1,numVal_BKS,1) - svXyzTrx_BKS_mirrored'; %对称映射
    v_NBKS=xyz0(:)*ones(1,numVal_NBKS,1) - svXyzTrx_NBKS';  

%   range = sqrt( sum(v.^2) );
    range_BKS=sqrt(sum(v_BKS.^2));
    v_BKS = v_BKS./(ones(3,1)*range_BKS); % line of sight unit vectors from sv to xo
    range_NBKS=sqrt(sum(v_NBKS.^2));
    v_NBKS = v_NBKS./(ones(3,1)*range_NBKS); % line of sight unit vectors from sv to xo

  % 发送信号时卫星的标号，位置，及其钟差
%   svPos=[prs(:,jSv),svXyzTrx,dtsv(:)];
%   dtsv=[dtsv_BKS;dtsv_NBKS];
     svPos_BKS=[prs_BKS(:,jSv),svXyzTrx_BKS_mirrored,dtsv_BKS(:)];
     svPos_NBKS=[prs_BKS(:,jSv),svXyzTrx_BKS_mirrored,dtsv_BKS(:)];
     
  %calculate the a-priori range residual
%   prHat = range(:) + bc - GpsConstants.LIGHTSPEED*dtsv;
    prHat_BKS =range_BKS(:)+ bc_BKS - GpsConstants.LIGHTSPEED*dtsv_BKS;
    prHat_NBKS =range_NBKS(:)+ bc_NBKS - GpsConstants.LIGHTSPEED*dtsv_NBKS;
  % Use of bc: bc>0 <=> pr too big <=> rangehat too big
  % Use of dtsv: dtsv>0 <=> pr too small
    
%    zPr = prs(:,jPr)-prHat; 
    zPr_BKS = prs_BKS(:,jPr)-prHat_BKS; 
    zPr_NBKS = prs_NBKS(:,jPr)-prHat_NBKS; 
    zPr=[zPr_BKS;zPr_NBKS];
    
    % H matrix = [unit vector,1] 
    % H矩阵扩充 原先的H矩阵是18 * 4 前三个维度是卫星的方向矢量，第四个维度是时间 
    H_BKS = [v_BKS', ones(numVal_BKS,1),zeros(numVal_BKS,1)]; 
    H_NBKS = [v_NBKS',zeros(numVal_NBKS,1), ones(numVal_NBKS,1)];
    H=[H_BKS;H_NBKS];
   
%   %z = Hx, premultiply by W: Wz = WHx, and solve for x:
    dx = pinv(Wpr*H)*Wpr*zPr;

  % update xo, xhat and bc
    xHat=xHat+dx;
    xyz0=xyz0(:)+dx(1:3);

%% 钟差校准
  bc_BKS=bc_BKS+dx(4);
%   bc_BKS=bc_BKS+dx(4)-ToF_deltaSeconds(i);
  bc_NBKS=bc_NBKS+dx(5);
  
  %Now calculate the a-posteriori range residual
  zPr = zPr-H*dx;
  
  % display info
  string=['whileCount = ',num2str(whileCount),', norm(dx) =',num2str(norm(dx))];
  disp(string)
  lla0=Xyz2Lla(xyz0(:)');
  stringPos=[' Position: ',num2str(lla0(1)),',',num2str(lla0(2)),',',num2str(lla0(3))];
  disp(stringPos)
  
%%  滤波处理，对于Tag覆盖不到的点进行消除
% 暂时行不通
%     distance=norm(TagLocXyz0'-xyz0);
% if norm(dx) < GnssThresholds.MAXDELPOSFORNAVM
%     if distance < 30
%         break;
%     else
% %         xyz0=TagLocXyz0'- dx(1:3);
% %     xyz0=TagLocXyz0';
% %       xHat(1:3)=dx(1:3);
%         continue;
%     end
% end

% 滤波与dx重新赋值，让dx从当前位置往Tag位置位移一个单位长度
    distance=norm(TagLocXyz0'-xyz0);
    if norm(dx) < GnssThresholds.MAXPRRUNCMPS  % MAXDELPOSFORNAVM 
        if distance < 30
            break;
        else
%         xyz0=TagLocXyz0'- dx(1:3);
%     xyz0=TagLocXyz0';
%       xHat(1:3)=dx(1:3);
            correctionVector=(TagLocXyz0-xyz0);%  ./norm(TagLocXyz0'-xyz0);
            dx(1:3)=correctionVector./10;
            xHat=xHat+dx;
            xyz0=xyz0(:)+dx(1:3);

            stringPos=[' Position: ',num2str(lla0(1)),',',num2str(lla0(2)),',',num2str(lla0(3))];
            disp(stringPos)
%         continue;
        end
    end
% correctionVector=(TagLocXyz0'-xyz0)./norm(TagLocXyz0'-xyz0);
end
%% 速度部分计算
% % Compute velocities ---------------------------------------------------------
% svXyzDot_BKS,dtsvDot_BKS
% rrMps = zeros(numVal,1);
rrMps_BKS = zeros(numVal_BKS,1);
rrMps_NBKS = zeros(numVal_NBKS,1);

for i=1:numVal_BKS
    %range rate = [satellite velocity] dot product [los from xo to sv]
    rrMps_BKS(i) = -svXyzDot_BKS(i,:)*v_BKS(:,i);
end

for i=1:numVal_NBKS
    %range rate = [satellite velocity] dot product [los from xo to sv]
    rrMps_NBKS(i) = -svXyzDot_NBKS(i,:)*v_NBKS(:,i);
end

% prrHat = rrMps + xo(8) - GpsConstants.LIGHTSPEED*dtsvDot;
prrHat_BKS = rrMps_BKS + xo(8) - GpsConstants.LIGHTSPEED*dtsvDot_BKS;
prrHat_NBKS = rrMps_NBKS + xo(8) - GpsConstants.LIGHTSPEED*dtsvDot_NBKS;

zPrr_BKS = prs_BKS(:,jPrr)-prrHat_BKS;
zPrr_NBKS = prs_NBKS(:,jPrr)-prrHat_NBKS;
zPrr=[zPrr_BKS;zPrr_NBKS];

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
% elseif  any(size(xo) ~= [8,1]) % 待求解矩阵xo
elseif  any(size(xo) ~= [10,1]) % 待求解矩阵xo
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

