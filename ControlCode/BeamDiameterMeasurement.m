% Assumes pm and zabers are already initialized as done in PiezoScan code

pmNread = 100;      % number of power measurments at each point
posPts  = 200;      % number of points

ptsStrt = 13;
ptsStop = 17.5; 

pmPos   = linspace(ptsStrt,ptsStop,posPts);
pmReads = nan(posPts,pmNread);


for jj = 1:posPts
    VFN_Zab_move(pmX, pmPos(jj));
    pause(0.5)
    for ii = 1:pmNread
        pmReads(jj,ii)=str2num(query(obj1, 'measure:power?'));
    end    
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

flnm = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\021919_FNM3\BeamDiam';

save(flnm, 'pmPos', 'pmReads')

datFl = strcat(flnm, 'data.txt');
fileID = fopen(datFl, 'w');
for jj = 1:posPts
    fprintf(fileID, '%6.4f %8.6e %8.6e\n', pmPos(jj), pmVals(jj), pmErr(jj));
end
fclose(fileID);
