% Assumes pm and zabers are already initialized as done in PiezoScan code

pmNread = 100;      % number of power measurments at each point
posPts  = 100;      % number of points

ptsStrt = 18.5;
ptsStop = 23; 

pmPos   = linspace(ptsStrt,ptsStop,posPts);
pmReads = nan(posPts,pmNread);


for jj = 1:posPts
    VFN_Zab_move(Zabs.pmX, pmPos(jj));
    pause(0.3)
    pmReads(jj,:)=VFN_PM100D_read(PM,pmNread);
end
if sum(isnan(pmReads))
    fprintf('nan present')
end
figure()
plot(pmPos, mean(pmReads,2))

% Calculate beam size:
pmVals = mean(pmReads,2);
pmErr = std(transpose(pmReads));
pmDev = gradient(gradient(pmVals));
[valMax, indMax] = max(pmDev(:));
[valMin, indMin] = min(pmDev(:));
fprintf('Beam Size: %f\n',pmPos(indMin)-pmPos(indMax))

figure()
plot(pmPos, pmDev)

flnm = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\200615_COV1\Realigned_BeamDiam635nmLasSamp2';

save(flnm, 'pmPos', 'pmReads')

datFl = strcat(flnm, 'data.txt');
fileID = fopen(datFl, 'w');
for jj = 1:posPts
    fprintf(fileID, '%6.4f %8.6e %8.6e\n', pmPos(jj), pmVals(jj), pmErr(jj));
end
fclose(fileID);
