% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\070119_VPL9';

% Experiment name for figures
expNm = '13BWS_Var640BWScan2';

VarPow = [.462 1.205 2.1218 3.041 3.74]
X = [6 12 18 24 30]
Y = [0.004209 0.01129 0.019956 0.028571 0.036220]

Y = Y ./ VarPow


figure();
scatter(X,Y)
titStr = sprintf('Eta_s Vs Bandwidth (Normalized)');
title(titStr)
xlabel(['Varia Bandwidth about 640nm'])
ylabel('Eta s')
saveas(gcf, [svFld filesep expNm '_BWGraphNorm_CoupMap.png'])