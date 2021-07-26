function [ distance ]=AnalysisCdf(DirName,FileName)
% AnalysisData: 4 x 10

load([DirName '/' FileName],'distance','hdop');

a = strsplit(FileName,'_');

% distance filter 
distance = disFilter(distance,50); 
length(distance.org)
length(distance.h)
end

function disFiltered = disFilter(distance,R) 
disFiltered.org = distance.org(find(distance.org < R));
disFiltered.h = distance.h(find(distance.h < R));

disFiltered.org = sort(disFiltered.org);
disFiltered.org = disFiltered.org(1:floor(length(disFiltered.org)/2));

disFiltered.h = sort(disFiltered.h);
disFiltered.h = disFiltered.h(1:floor(length(disFiltered.h)/2));
end

