function [gnssMeasOut] = SatFilter(gnssMeas,SatNunmer,draw)
%SATFILTER 此处显示有关此函数的摘要
%   此处显示详细说明

if nargin < 3
   draw =1;    
end

M = length(gnssMeas.Svid);
if M < SatNunmer
   error(['Not enough Satlite: ' num2str(M) '<' num2str(SatNumer)])
end

if draw
    h3 = figure;
    colors = PlotCno(gnssMeas,'test');
end

% 信号持续时间长，强度大，选取整个采样时段内平均值最强的信号
iNan = find(isnan(gnssMeas.Cn0DbHz));
Cn0iDbHz = gnssMeas.Cn0DbHz;
Cn0iDbHz(iNan) = 0;
Ecn0 = mean(Cn0iDbHz);
[Esort,iMax] = sort(Ecn0,'descend');
iSatOut = iMax(1:SatNunmer);

if draw
    timeSeconds = gnssMeas.FctSeconds-gnssMeas.FctSeconds(1);% elapsed time in seconds
    Nmax = SatNunmer;
    figure
    for i=1:Nmax
        Cn0iDbHz = gnssMeas.Cn0DbHz(:,iMax(i));
        iF = find(isfinite(Cn0iDbHz));
        if ~isempty(iF)
            ti = timeSeconds(iF(end));
            h = plot(timeSeconds,Cn0iDbHz);
            hold on
            set(h,'Color',colors(iMax(i),:));
            ts = int2str(gnssMeas.Svid(iMax(i)));
            text(ti,Cn0iDbHz(iF(end)),ts,'Color',colors(iMax(i),:));
        end
    end
end

% Extract Satlites
cellGnssMeas = struct2cell(gnssMeas);
fieldMeas = fieldnames(gnssMeas);
NumCell = length(cellGnssMeas);
cellGnssMeasOut = cellGnssMeas;

for i=1:NumCell
    if size(cellGnssMeas{i,1},2) > 1
       cellGnssMeasOut{i,1} = cellGnssMeasOut{i,1}(:,iSatOut) ;
    end
end
gnssMeasOut = cell2struct(cellGnssMeasOut, fieldMeas, 1);

end

