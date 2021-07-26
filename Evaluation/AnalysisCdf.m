function [ distance ]=AnalysisCdf(DirName,FileName)
% AnalysisData: 4 x 10

load([DirName '/' FileName],'distance','hdop');

a = strsplit(FileName,'_');

% distance filter 
distance = disFilter(distance,100); 
length(distance.org)
length(distance.h)
end

function disFiltered = disFilter(distance,R) 
disFiltered.org = distance.org(find(distance.org < R));
% disFiltered.bks = distance.bks(find(distance.bks < R));
% disFiltered.nbks = distance.nbks(find(distance.nbks < R));
disFiltered.h = distance.h(find(distance.h < R));
end

