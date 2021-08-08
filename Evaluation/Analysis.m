function [ AnalysisData ]=Analysis(DirName,FileName)
% AnalysisData: 4 x 10
Ntype  = 4;
Nvalue = 10;

load([DirName '/' FileName],'distance','hdop');
% Point Num Aver Median Std Max Min CI_Upper CI_lower
jIndex = {1 2 3 4 5 6 7 8 9 10 11};
[jPoint, jNum, jAver, jMed, jStd, jMax, jMin, jCI_Upper, jCI_lower, jHdop] = jIndex{:};

AnalysisData = zeros(Ntype,Nvalue);

a = strsplit(FileName,'_');
AnalysisData(:,jPoint) = ones(Ntype,1) * str2double(a{1}(isstrprop(a{1},'digit')));

% distance filter 
distance = disFilter(distance,15); 

AnalysisData(1,jNum:jMin) = ExtractDis(distance.org);
AnalysisData(2,jNum:jMin) = ExtractDis(distance.bks);
AnalysisData(3,jNum:jMin) = ExtractDis(distance.nbks);
AnalysisData(4,jNum:jMin) = ExtractDis(distance.h);

SEM = AnalysisData(:,jStd) ./ sqrt(AnalysisData(:,jNum));            % Standard Error
ts = 1.96;  % 95%
AnalysisData(:,jCI_Upper) = AnalysisData(:,jAver) + ts*SEM;% Confidence Intervals
AnalysisData(:,jCI_lower) = AnalysisData(:,jAver) - ts*SEM;

% hdop process
AnalysisData(1,jHdop) = ExtractHdop(hdop.org);
AnalysisData(2,jHdop) = ExtractHdop(hdop.bks);
AnalysisData(3,jHdop) = ExtractHdop(hdop.nbks);
AnalysisData(4,jHdop) = ExtractHdop(hdop.h);
if ~isreal(AnalysisData(:,jHdop)) %
    FileName
    AnalysisData(:,jHdop)
end
end

function disFiltered = disFilter(distance,threshold) 
disFiltered.org = distance.org(find(distance.org < 300));
% disFiltered.org = distance.org;
disFiltered.bks = distance.bks(find(distance.bks < threshold));
% disFiltered.nbks = distance.nbks(find(distance.nbks < threshold));
disFiltered.nbks = distance.org(find(distance.org < threshold));
disFiltered.h = distance.h(find(distance.h < threshold));
end

function [AnalysisData] = ExtractDis(distance)
% Point Num Aver Median Std Max Min CI_Upper CI_lower
jIndex = {1 2 3 4 5 6 7 8 9 10 11};
[jPoint, jNum, jAver, jMed, jStd, jMax, jMin, jCI_Upper, jCI_lower] = jIndex{:};

if isempty(distance)
    AnalysisData = zeros(1,jMin + 1 - jNum);
    AnalysisData(1,:) = missing;
    return;
end

AnalysisData(jNum) = length(distance);
AnalysisData(jAver) = mean(distance);
AnalysisData(jMed) = median(distance);
AnalysisData(jStd) = std(distance);
AnalysisData(jMax) = max(distance);
AnalysisData(jMin) = min(distance);

AnalysisData = AnalysisData(jNum:jMin);
end

function [AnalysisData] = ExtractHdop(hdop)
% Point Num Aver Median Std Max Min CI_Upper CI_lower
jIndex = {1 2 3 4 5 6 7 8 9 10 11};
[jPoint, jNum, jAver, jMed, jStd, jMax, jMin, jCI_Upper, jCI_lower, jHdop] = jIndex{:};

AnalysisData(jHdop) = mean(hdop(find(hdop < inf)));

AnalysisData = AnalysisData(jHdop);
end