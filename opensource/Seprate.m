function  [gnssMeas_BKS, gnssMeas_NBKS]=Seprate(gnssMeas,prFileName)
% figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% figure;
% PlotCno(gnssMeas,prFileName,colors);

% pick out the strongest sat signal as standard
[strongestSig,SigIndex] = max(gnssMeas.Cn0DbHz);
[MaxCn0DbHz,Index] = max(strongestSig);
StongestSvid = gnssMeas.Svid(Index);

% Index ;
Index
gnssCnoMax = gnssMeas.Cn0DbHz(:,Index);

% set the nan as 0
for i=1:length(gnssCnoMax)
    if isnan(gnssCnoMax(i)) 
       gnssCnoMax(i)=0; 
    end
end

% plot for scale view
figure;
plot(gnssCnoMax);

% window length 4
% split data into two parts
LenGnss = length(gnssCnoMax);
% pickUp -> backscattered sat signal  BKS
% pickDn -> original sat signal       NBKS
pickUp = zeros(1, LenGnss);
pickDn = zeros(1, LenGnss);
% midThre = 0;
cntPick = ones(1,2);
% for i=1:LenGnss
%     if i+9 <= LenGnss
%         if mod(i,10) == 1
%            dataBuf = gnssCnoMax(i:i+9); 
%            midThre = median(dataBuf); 
%         end
%     end
%     
%     if gnssCnoMax(i) >= midThre - 1
%         pickUp(i) = 1;
%         cntPick(1) = cntPick(1) + 1;
%     else
%         pickDn(i) = 1;
%         cntPick(2) = cntPick(2) + 1;
%     end
%     
% end

CnoSmooth = smooth(gnssCnoMax, 7);
for i=1:LenGnss
    if gnssCnoMax(i) > CnoSmooth(i)
        pickUp(i) = 1;
        cntPick(1) = cntPick(1) + 1;
    else
        pickDn(i) = 1;
        cntPick(2) = cntPick(2) + 1;
    end
    
end

% plot to check seperation
gnssMeasSep.Bks = pickUp .* gnssCnoMax';
gnssMeasSep.Org = pickDn .* gnssCnoMax';
gnssMeasSep.smooth = smooth(gnssCnoMax, 7);
figure
plot(gnssCnoMax);
hold on
plot(gnssMeasSep.smooth)
plot(gnssMeasSep.Bks, 'r*');
plot(gnssMeasSep.Org, 'g*');
hold off



% Seperate BKS and NBKS 
iBKS = pickUp==1;
iORG = pickDn==1;
cellGnssMeas = struct2cell(gnssMeas);
fieldMeas = fieldnames(gnssMeas);
NumCell = length(cellGnssMeas);
cellGnssMeas_BKS = cellGnssMeas;
cellGnssMeas_NBKS = cellGnssMeas;
for i=1:NumCell
    if length(cellGnssMeas{i,1}(:,1)) > 1
       cellGnssMeas_BKS{i,1} = cellGnssMeas{i,1}(iBKS,:) ;
       cellGnssMeas_NBKS{i,1} = cellGnssMeas{i,1}(iORG,:) ;
    end
end
gnssMeas_BKS = cell2struct(cellGnssMeas_BKS, fieldMeas, 1);
gnssMeas_NBKS = cell2struct(cellGnssMeas_NBKS, fieldMeas, 1);

% save seperate date
% save( 'gnssMean.mat','gnssMeas', 'gnssMeas_BKS', 'gnssMeas_NBKS');

% figure;
% hold on
% plot(gnssMeas_BKS.Cn0DbHz(:,1));

figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
PlotCno(gnssMeas_BKS,prFileName,colors);

end