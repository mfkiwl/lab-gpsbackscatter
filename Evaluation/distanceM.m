function [distM] = distanceM(allLlaDegDegM,llaTrueDegDegM)
%DISTANCE 计算各个经纬度点到某点的平面距离
%   lla 坐标系 转 参考点为中心的 ned 坐标系算距离

nedM = Lla2Ned(allLlaDegDegM,llaTrueDegDegM);
iFi = isfinite(allLlaDegDegM(:,1));
distM = sqrt(sum(nedM(iFi,1:2).^2,2));

end

