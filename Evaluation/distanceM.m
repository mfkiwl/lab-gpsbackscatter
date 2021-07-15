function [distM] = distanceM(allLlaDegDegM,llaTrueDegDegM)
%DISTANCE 此处显示有关此函数的摘要
%   此处显示详细说明

nedM = Lla2Ned(allLlaDegDegM,llaTrueDegDegM);
iFi = isfinite(allLlaDegDegM(:,1));
distM = sqrt(sum(nedM(iFi,1:2).^2,2));

end

