function [distM] = distanceM(allLlaDegDegM,llaTrueDegDegM)
%DISTANCE ���������γ�ȵ㵽ĳ���ƽ�����
%   lla ����ϵ ת �ο���Ϊ���ĵ� ned ����ϵ�����

nedM = Lla2Ned(allLlaDegDegM,llaTrueDegDegM);
iFi = isfinite(allLlaDegDegM(:,1));
distM = sqrt(sum(nedM(iFi,1:2).^2,2));

end

