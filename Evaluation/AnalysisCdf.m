function [ distance hdop]=AnalysisCdf(DirName,FileName)
% AnalysisData: 4 x 10

load([DirName '/' FileName],'distance','hdop');

a = strsplit(FileName,'_');

% distance filter 
distance = disFilter(distance,50); 
length(distance.org)
length(distance.h)

% hdop filter 
hdop = hdopFilter(hdop,30); 
length(hdop.org)
length(hdop.h)
end

function disFiltered = disFilter(distance,R) 
disFiltered.org = distance.org(find(distance.org < R));
disFiltered.h = distance.h(find(distance.h < R));

disFiltered.org = sort(disFiltered.org);
disFiltered.org = disFiltered.org(1:floor(length(disFiltered.org)/2));

disFiltered.h = sort(disFiltered.h);
disFiltered.h = disFiltered.h(1:floor(length(disFiltered.h)/2));
end

function hdopFiltered = hdopFilter(hdop,R) 
hdop.org = hdop.org(isfinite(hdop.org));
hdopFiltered.org = hdop.org(find(hdop.org < R));
hdopFiltered.h = hdop.h(find(real(hdop.h) < R));

hdopFiltered.org = hdopFiltered.org;
hdopFiltered.h = real(hdopFiltered.h);

% hdopFiltered.org = sort(hdopFiltered.org);
% hdopFiltered.org = hdopFiltered.org(1:floor(length(hdopFiltered.org)/2));
% 
% hdopFiltered.h = sort(hdopFiltered.h);
% hdopFiltered.h = hdopFiltered.h(1:floor(length(hdopFiltered.h)/2));
end
