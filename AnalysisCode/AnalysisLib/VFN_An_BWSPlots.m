% Plotting Script for analyzing the data for the monochromatic OSA paper
% NOTE: this script isn't capable of handling multiple wavelength centers
close all
clear all

%% Load data
% Default fonts for saved figures
fontsize1 = 26;
fontsize2 = 18;

%Addpath to matlab figure exporting functions
%addpath('C:\Users\danie\Documents\MATLAB\SupportPackages\altmany-export_fig-d570645')

% add path to intelligent polar transform function 
%addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib')

% add path to analysis library
%addpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\AnalysisCode\AnalysisLib')
addpath('C:\Users\Thomas\Desktop\VFN-Lab-master\VFN-Lab-master\AnalysisCode\AnalysisLib')
addpath('D:\Programs\Matlab\altmany-export_fig-b1a7288')


% Folder with data
%fldnm = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling';
fldnm = 'C:\Users\Thomas\Desktop\vfn_test';

%flSimm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\ProcessedData\';
% Folder for saving results
%savefld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\071019_VAR1test';
savefld = 'C:\Users\Thomas\Desktop\bwsplots_test';

% Scaling factor for converting between femto gain settings
an_params.gnFactor = 9.97;

% Throughput of the final lens (used for the red thorlabs pm normalization)
an_params.lensTHPT = 0.9973;

% Conversion factor between red thorlabs PM in uW to femto volts at 10^6
an_params.uw2fmto = .9;%0.6318;        %P_femto[v]/P_red[uW]
an_params.rdOfmto = 1/an_params.uw2fmto;   % Ratio as defined by nem (red[uW]/fmto[V]

% Conversion from piezo volts to microns
an_params.V2um = 1/8.72; %<-- Value determined via PZT characterization
                  % Spec value: 1/10;

% Other normalization terms
an_params.FR = 0.9654;      % Fresnel reflection term in Nem's equation
an_params.PL = 0.9966;      % Propagation loss term in nem's equation.

% Read all .fits files in the directory
%nmWpath = { [fldnm '\071019_VAR1test\*.fits']};
nmWpath = { [fldnm '\*.fits']};
an_params.STRNMS = VFN_An_getSTRNMS(nmWpath);

% 
filetag = '74BWS_Var780Dat7_fullCube';
% Sections to Run
Sec2Run = [ 
    false;          %Section 1  ** Donut Sensor Matrix
    false;          %Section 2  ** Donut Sensor Matrix (Log scale)
    true;          %Section 3  ** etaS linear plot
    false;          %Section 4  ** etaS ylog plot
    false;          %Section 5  ** etaS xylog plot
    true;          %Section 6  ** etaP linear plot
    false;          %Section 7  ** etaP ylog plot
    false;          %Section 8  ** etaP xylog plot
    false;          %Section END  -- REFERENCE ONLY; DO NOT ENABLE
    ];
    
% Flag to print locations
isPrintLoc = false;

% Flag to mark whether (.mat) data should be saved
isSave = false;

% Flag to mark whether figures should be saved
isSaveFig = false;

if ~isPrintLoc
    warning('Location printing is disabled with "isPrintLoc"')
end

[BWs, Wvls] = VFN_An_getWvlScanPos(filetag, an_params);
an_params.BWs = BWs;
an_params.Wvls = Wvls;

% Conversion factor from microns to L/D in focal plane: [L/D / micron]
an_params.um2LD = 1/((Wvls(1)/1000)*10.96/2.1); % Conversion to lambda/D: 1/(lambda*F#) 
                       % f = 10.96 taken from Thorlabs A397 datasheet

% Manual input for normalization values, leave empty if values are in
    %fits file
an_params.NormVals = [6.46973 13.43771 20.28842 26.57895 32.91610 38.78649 44.59445 49.96976 55.70637 60.81887];
%% Section 1: Donut Sensor Matrix
if Sec2Run(1)
	disp("______Section 1:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(nrmFine(i));

        %- Print eta_s info
        fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine(i))
        if isPrintLoc
            fprintf([' eta_s Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indFine(i,2), indFine(i,1), indFine(i,4), indFine(i,5), indFine(i,6));
        end
        %- Calculate eta_p
        %[eta_pFine(i), indPFine(i)] = VFN_An_getEta_p(nrmFine(i));
        %fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine(i));
        %if isPrintLoc
        %    fprintf([' eta_p Location: ', ...
        %             '\n  X index = %3i', ...
        %             '\n  Y index = %3i', ...
        %             '\n  Foc ind = %3i', ...
        %             '\n  VX ind  = %3i', ...
        %             '\n  VY ind  = %3i\n'], ...
        %             indPFine(i,2), indPFine(i,1), indPFine(i,4), indPFine(i,5), indPFine(i,6));
        %end
    end
    
    for j = 1:length(BWs)
        %- Get axes for image (relative to null point):
        [xaxFine, yaxFine] = VFN_An_getAxes(nmFine, squeeze(indFine(j,:)),an_params);
        
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        %fighAsym = VFN_An_fitsDisp(nrmFine, xaxFine, yaxFine, figh);
        fighAsym = VFN_An_fitsDisp(squeeze(nrmFine(j,:,:)), xaxFine, yaxFine, figh);
        %imagesc(squeeze(nrmFine(j,:,:)))
        title(['Scan at ' num2str(Wvls(1)) 'nm with ' num2str(BWs(j)) 'nm bandwidth'])%, 'Position', [-0.0183 1.1178-0.05 0])
        xlabel('\lambda/D')
        ylabel('\lambda/D')
        set(gca,'fontsize',fontsize2)
        
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_Lin_BW' num2str(BWs(j)) '_CoupMap.png'],'-r300', '-painters');
        end
    end
    
       
    fprintf("______End Section 1\n\n")
end

%% Section 2: Donut Sensor Matrix (Log)
if Sec2Run(2)
    disp('______Section 2:______')
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(nrmFine(i));

        %- Print eta_s info
        fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine(i))
        if isPrintLoc
            fprintf([' eta_s Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indFine(i,2), indFine(i,1), indFine(i,4), indFine(i,5), indFine(i,6));
        end
        %- Calculate eta_p
        %[eta_pFine(i), indPFine(i)] = VFN_An_getEta_p(nrmFine(i));
        %fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine(i));
        %if isPrintLoc
        %    fprintf([' eta_p Location: ', ...
        %             '\n  X index = %3i', ...
        %             '\n  Y index = %3i', ...
        %             '\n  Foc ind = %3i', ...
        %             '\n  VX ind  = %3i', ...
        %             '\n  VY ind  = %3i\n'], ...
        %             indPFine(i,2), indPFine(i,1), indPFine(i,4), indPFine(i,5), indPFine(i,6));
        %end
    end
    
    for j = 1:length(BWs)
        %- Get axes for image (relative to null point):
        [xaxFine, yaxFine] = VFN_An_getAxes(nmFine, squeeze(indFine(j,:)),an_params);
        
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        for k = 1:length(nrmFine(1,1,:))
            for i = 1:length(nrmFine(1,1,:))
                nrmFine(j,k,i) = log10(nrmFine(j,k,i));
            end
        end
        fighAsym = VFN_An_fitsDisp(squeeze(nrmFine(j,:,:)), xaxFine, yaxFine, figh);
        %imagesc(squeeze(nrmFine(j,:,:)))
        title(['Scan at ' num2str(Wvls(1)) 'nm with ' num2str(BWs(j)) 'nm bandwidth [log]'])%, 'Position', [-0.0183 1.1178-0.05 0])
        xlabel('\lambda/D')
        ylabel('\lambda/D')
        set(gca,'fontsize',fontsize2)
        
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_Log_BW' num2str(BWs(j)) '_CoupMap.png'],'-r300', '-painters');
        end
        
    end
    
    fprintf("______End Section 2\n\n")
end

%% Section 3: Eta S linear plot
if Sec2Run(3)
	disp("______Section 3:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(nrmFine(i));

        %- Print eta_s info
        fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine(i))
        if isPrintLoc
            fprintf([' eta_s Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indFine(i,2), indFine(i,1), indFine(i,4), indFine(i,5), indFine(i,6));
        end
    end
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,eta_sFine,25, 'filled')
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'])%, [-0.0183 1.1178-0.05 0])
        xlabel(['Varia Bandwidth about ' num2str(Wvls(j)) 'nm [\delta\lambda/\lambda_0]'])
        ylabel('\eta_{s}')
        set(gca,'fontsize',fontsize2)
        
        P = polyfit(xax,eta_sFine,1);
        yfit = P(1)*xax+P(2);
        hold on;
        plot(xax,yfit,'r-.');
        text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_Lin_etaS_CoupMap.png'],'-r300', '-painters');
        end
    end
       
    fprintf("______End Section 3\n\n")
end

%% Section 4: Eta S log plot
if Sec2Run(4)
    disp('______Section 4:______')
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(nrmFine(i));
        eta_sFine(i) = log10(eta_sFine(i));

        %- Print eta_s info
        fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine(i))
        if isPrintLoc
            fprintf([' eta_s Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indFine(i,2), indFine(i,1), indFine(i,4), indFine(i,5), indFine(i,6));
        end
    end
    %- Show Null frame
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %fighAsym = VFN_An_fitsDisp(nrmFine, xaxFine, yaxFine, figh);
    xax = BWs ./ Wvls(1);
    scatter(xax,eta_sFine,25, 'filled')
    title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm [log]'])%, 'Position', [-0.0183 1.1178-0.05 0])
    xlabel(['Varia Bandwidth about ' num2str(Wvls(j)) 'nm [\delta\lambda/\lambda_0]'])
    ylabel('\eta_{s} [log]')
    set(gca,'fontsize',fontsize2)

    P = polyfit(xax,eta_sFine,1);
    yfit = P(1)*xax+P(2);
    hold on;
    plot(xax,yfit,'r-.');
    text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])

    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmFine);
        export_fig([savefld filesep tmpnm '_An_BWGraph_Log_etaS_CoupMap.png'],'-r300', '-painters');
    end
    
    fprintf("______End Section 4\n\n")
end


%% Section 5: Donut scan with vortex centered for symmetry - MNO1 [10]
%  --- The null/planet coupling are bad but the symmetry is good
if Sec2Run(5)
	disp("______Section 5:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(nrmFine(i));

        %- Print eta_s info
        fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine(i))
        if isPrintLoc
            fprintf([' eta_s Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indFine(i,2), indFine(i,1), indFine(i,4), indFine(i,5), indFine(i,6));
        end
        %- Calculate eta_p
        %[eta_pFine(i), indPFine(i)] = VFN_An_getEta_p(nrmFine(i));
        %fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine(i));
        %if isPrintLoc
        %    fprintf([' eta_p Location: ', ...
        %             '\n  X index = %3i', ...
        %             '\n  Y index = %3i', ...
        %             '\n  Foc ind = %3i', ...
        %             '\n  VX ind  = %3i', ...
        %             '\n  VY ind  = %3i\n'], ...
        %             indPFine(i,2), indPFine(i,1), indPFine(i,4), indPFine(i,5), indPFine(i,6));
        %end
    end
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,eta_sFine,25, 'filled')
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'])%, [-0.0183 1.1178-0.05 0])
        xlabel(['Varia Bandwidth about ' num2str(Wvls(j)) 'nm [\delta\lambda/\lambda_0]'])
        ylabel('\eta_{s}')
        set(gca,'fontsize',fontsize2)
        
        P = polyfit(xax,eta_sFine,1);
        yfit = P(1)*xax+P(2);
        hold on;
        plot(xax,yfit,'r-.');
        text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_Lin_etaS_CoupMap.png'],'-r300', '-painters');
        end
    end
       
    fprintf("______End Section 5\n\n")
end 

%% Section 6: Eta P linear plot
if Sec2Run(6)
	disp("______Section 6:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Calculate eta_p
        [eta_pFine(i), indPFine(i,:)] = VFN_An_getEta_p(nrmFine(i));
        fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine(i));
        if isPrintLoc
            fprintf([' eta_p Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indPFine(i,2), indPFine(i,1), indPFine(i,4), indPFine(i,5), indPFine(i,6));
        end
    end
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,eta_pFine,25, 'filled')
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'])%, [-0.0183 1.1178-0.05 0])
        xlabel(['Varia Bandwidth about ' num2str(Wvls(j)) 'nm [\delta\lambda/\lambda_0]'])
        ylabel('\eta_{p}')
        set(gca,'fontsize',fontsize2)
        
        P = polyfit(xax,eta_pFine,1);
        yfit = P(1)*xax+P(2);
        hold on;
        plot(xax,yfit,'r-.');
        text(xax(1),eta_pFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_Lin_etaP_CoupMap.png'],'-r300', '-painters');
        end
    end
       
    fprintf("______End Section 6\n\n")
end 

%% Section 7: Eta P log plot
%  NOTE: extracts specific frames from the scanned cube
if Sec2Run(7)
	disp("______Section 7:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['vfn_test' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Calculate eta_p
        [eta_pFine(i), indPFine(i,:)] = VFN_An_getEta_p(nrmFine(i));
        eta_pFine(i) = log10(eta_pFine(i));
        fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine(i));
        if isPrintLoc
            fprintf([' eta_p Location: ', ...
                     '\n  X index = %3i', ...
                     '\n  Y index = %3i', ...
                     '\n  Foc ind = %3i', ...
                     '\n  VX ind  = %3i', ...
                     '\n  VY ind  = %3i\n'], ...
                     indPFine(i,2), indPFine(i,1), indPFine(i,4), indPFine(i,5), indPFine(i,6));
        end
    end
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,eta_pFine,25, 'filled')
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'])%, [-0.0183 1.1178-0.05 0])
        xlabel(['Varia Bandwidth about ' num2str(Wvls(j)) 'nm [\delta\lambda/\lambda_0]'])
        ylabel('\eta_{p}')
        set(gca,'fontsize',fontsize2)
        
        P = polyfit(xax,eta_pFine,1);
        yfit = P(1)*xax+P(2);
        hold on;
        plot(xax,yfit,'r-.');
        text(xax(1),eta_pFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_Lin_etaP_CoupMap.png'],'-r300', '-painters');
        end
    end
    
       
    fprintf("______End Section 7\n\n")
end 

%% Section END: NOT USED - Reference only
%---  This is an example of how I used to do the azimuthal averaging
if Sec2Run(end)
    disp('______Section END:______')
    
    %-- Calculate max values around donut
    %find x,y coords for null
    [nlRow, nlCol] = find(fits == nlVal);
    %calculate polar coords for cartesian
    [theta, r] = fitsPol(fits, nlCol, nlRow);
    %preallocate loop values
    j   = 0 : 0.1 : 1;          %radii over which to analyze
    mns = -5*ones(length(j),1); %array to hold means. *-5 to check for errs
    tol = 0.05;     %Tolerance for radii 
    %iterate through radii and calculate average
    for k = 1:length(j)
        i       = j(k);
        %uncomment below to plot circles
        %fitspl  = fits;
        %fitspl((i-tol<=r)&(r<=i+tol)) = 1000;
        %fitsDisp(fitspl);
        mns(k) = mean(fits((i-tol<=r)&(r<=i+tol)));
    end
    
    %-- Display
    %Average for the donut ring with the peak coupling
    mnpkVal = max(mns);
    fprintf('--Radius of peak average: r = %f\n', j(mns==mnpkVal))
    
    %-- Plot average coupling for varying radii
    %Linear Scale
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
    plot(0:imLim/(length(j)-1):imLim,((mns-avgBkgd)/fibtip)*100,'-o', 'Color', [0.8500    0.3250    0.0980])
    title('Average Coupling')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
    set(gca,'fontsize',fontsize2)
    set(findall(gca, 'Type', 'Line'),'LineWidth', 2);
    grid on
    %export_fig([savefld filesep 'ActualCouplingLine.png'],'-r300', '-painters');
    
    fprintf('______End Section END\n\n')
end

%% Save data as needed
if isSave
    savenm = [savefld filesep 'PlotData_' datestr(now,'mmddyy_HH_MM') '.mat'];
    % Save something in the file so that we can append later
    save(savenm, 'gnFactor', 'lensTHPT', 'rdOfmto', 'um2LD', 'V2um', 'FR', 'PL', ...
                 'nrm*', 'xax*', 'yax*');
end