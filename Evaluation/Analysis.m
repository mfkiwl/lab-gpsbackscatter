function [ AnalysisData ]=Analysis(FileName)
% AnalysisData: 3 x 9

load(FileName,'distance');
% Point Num Aver Median Std Max Min CI_Upper CI_lower
jIndex = {1 2 3 4 5 6 7 8 9 10 11};
[jPoint, jNum, jAver, jMed, jStd, jMax, jMin, jCI_Upper, jCI_lower] = jIndex{:};
Ntype  = 3;
Nvalue = 9;

AnalysisData = zeros(Ntype,Nvalue);

a = strsplit(FileName,'_');
AnalysisData(:,jPoint) = ones(3,1) * str2double(a{1}(isstrprop(a{1},'digit')));

% filter
distance = disFilter(distance); 

if isempty(distance.org)
    AnalysisData(:,:) = missing;
    return;
end
AnalysisData(1,jNum:jMin) = ExtractDis(distance.org);
AnalysisData(2,jNum:jMin) = ExtractDis(distance.bks);
AnalysisData(3,jNum:jMin) = ExtractDis(distance.nbks);

SEM = AnalysisData(:,jStd) ./ sqrt(AnalysisData(:,jNum));            % Standard Error
ts = 1.96;  % 95%
AnalysisData(:,jCI_Upper) = AnalysisData(:,jAver) + ts*SEM;% Confidence Intervals
AnalysisData(:,jCI_lower) = AnalysisData(:,jAver) - ts*SEM;
end

function disFiltered = disFilter(distance) 
disFiltered.org = distance.org(find(distance.org < 30));
disFiltered.bks = distance.bks(find(distance.bks < 30));
disFiltered.nbks = distance.nbks(find(distance.nbks < 30));
end

function [AnalysisData] = ExtractDis(distance)
% Point Num Aver Median Std Max Min CI_Upper CI_lower
jIndex = {1 2 3 4 5 6 7 8 9 10 11};
[jPoint, jNum, jAver, jMed, jStd, jMax, jMin, jCI_Upper, jCI_lower] = jIndex{:};

AnalysisData(jNum) = length(distance);
AnalysisData(jAver) = mean(distance);
AnalysisData(jMed) = median(distance);
AnalysisData(jStd) = std(distance);
AnalysisData(jMax) = max(distance);
AnalysisData(jMin) = min(distance);

AnalysisData = AnalysisData(jNum:jMin);
end