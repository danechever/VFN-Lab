% Plotting Script for analyzing the data for the monochromatic OSA paper
close all
clear all

%% Load data
% Default fonts for saved figures
fontsize1 = 26;
fontsize2 = 18;

%Addpath to matlab figure exporting functions
%addpath('C:\Users\danie\Documents\MATLAB\SupportPackages\altmany-export_fig-d570645')

% add path to intelligent polar transform function 
addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib')

% add path to analysis library
addpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\AnalysisLib')

% Folder with data
fldnm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\PupilVFN\CouplingCurves';
flSimm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\ProcessedData\';
% Folder for saving results
savefld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\ProcessedData';

% Scaling factor for converting between femto gain settings
an_params.gnFactor = 9.97;

% Throughput of the final lens (used for the red thorlabs pm normalization)
an_params.lensTHPT = 0.9973;

% Conversion factor between red thorlabs PM in uW to femto volts at 10^6
an_params.uw2fmto = 0.6318;        %P_femto[v]/P_red[uW]
an_params.rdOfmto = 1/an_params.uw2fmto;   % Ratio as defined by nem (red[uW]/fmto[V]

% Conversion factor from microns to L/D in focal plane: [L/D / micron]
an_params.um2LD = 1/(0.635*10.96/2.1); % Conversion to lambda/D: 1/(lambda*F#) 
                       % f = 10.96 taken from Thorlabs A397 datasheet

% Conversion from piezo volts to microns
an_params.V2um = 1/8.72; %<-- Value determined via PZT characterization
                  % Spec value: 1/10;

% Other normalization terms
an_params.FR = 0.9654;      % Fresnel reflection term in Nem's equation
an_params.PL = 0.9966;      % Propagation loss term in nem's equation.

% Read all .fits files in the directory
nmWpath = { [fldnm '\021619_FNM2\*.fits'];
            [flSimm '\*.fits'];
            [fldnm '\*_FTO3\*.fits']; 
            [fldnm '\*_FTO4\*.fits']};
an_params.STRNMS = VFN_An_getSTRNMS(nmWpath);

% Sections to Run
Sec2Run = [ 
    false;          %Section 1  ** FNM2 [21] Best null
    true;          %Section 2  ** FTO3 [31] Best null (old Fnum - October)
    true;          %Section 3  ** FNM2 [21] Radial average
    false;          %Section 4  -- FTO4 [11] (Ultra Fine, Full FOV, Asymmetric)
    false;          %Section 5  -- MNO1 [10] Symmetric donut
    true;          %Section 6  ** FNM2 [23] Peak Coupling on PSF
    false;          %Section 7  -- MNO1 [5] Circular, Symmetric Donuts
    false;          %Section END  -- REFERENCE ONLY; DO NOT ENABLE
    ];
    
% Flag to print locations
isPrintLoc = true;

% Flag to mark whether (.mat) data should be saved
isSave = false;

% Flag to mark whether figurs should be saved
isSaveFig = false;

if ~isPrintLoc
    warning('Location printing is disabled with "isPrintLoc"')
end

%% Section 1: Best null full map - FNM2 [21]
if Sec2Run(1)
	disp("______Section 1:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['021619_FNM2' filesep '21Don_FineFull2'];
    nrmFine = VFN_An_fitsNormed(nmFine, [1,1,1], an_params);

    %- Extract eta_s and location
    [eta_sFine, indFine] = VFN_An_getEta_s(nrmFine);
    
    %- Print eta_s info
    fprintf('\nEta_s in Fine (full) Frame:    %f\n',eta_sFine)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFine(2), indFine(1), indFine(4), indFine(5), indFine(6));
    end

    %- Calculate eta_p
    [eta_pFine, indPFine] = VFN_An_getEta_p(nrmFine);
    fprintf('\nEta_p in Fine (full) Frame:    %f\n',eta_pFine);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPFine(2), indPFine(1), indPFine(4), indPFine(5), indPFine(6));
    end
    
    %- Get axes for image (relative to null point):
    [xaxFine, yaxFine] = VFN_An_getAxes(nmFine, indFine);
        
    %- Show Null frame
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    fighAsym = VFN_An_fitsDisp(nrmFine, xaxFine, yaxFine, figh);
    title('Fine Frame - Full: [FNM2 #21]', 'Position', [-0.0183 1.1178-0.05 0])
    set(gca,'fontsize',fontsize2)
    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmFine);
        export_fig([savefld filesep tmpnm '_CoupMap.png'],'-r300', '-painters');
    end
       
    fprintf("______End Section 1\n\n")
end

%% Section 2: Best null w/ Full Donut (October - old Fnumber) - FTO3 [31]
if Sec2Run(2)
    disp('______Section 2:______')
    %*** NOTE: CHANGE um2LD temporarily when analyzing data w/ old F-number
    um2LD_OLD = an_params.um2LD;
    an_params.um2LD = 1/(0.635*11/3.59); 
    
    %-- Get the normalized data (full eta map)
    nmFull  = ['101718_FTO3' filesep '31Don_NullCon'];
    nrmFull = VFN_An_fitsNormed(nmFull, [1,1,1], an_params);

    %-- Extract eta_s and location
    [eta_sFull, indFull] = VFN_An_getEta_s(nrmFull);
    
    %-- Print eta_s info
    fprintf('\nEta_s in Full Frame:    %f\n',eta_sFull)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFull(2), indFull(1), indFull(4), indFull(5), indFull(6));
    end
        
    %-- Calculate eta_p
    [eta_pFull, indPFull] = VFN_An_getEta_p(nrmFull);
    fprintf('\nEta_p in Full Frame:    %f\n',eta_pFull);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPFull(2), indPFull(1), indPFull(4), indPFull(5), indPFull(6));
    end

    %-- Get axes for image (relative to null point):
    [xaxFull, yaxFull] = VFN_An_getAxes(nmFull, indFull, an_params);
       
    %-- Show Full frame
    fighFull = VFN_An_fitsDisp(nrmFull, xaxFull, yaxFull);
    title('Full Frame Best Null: [FTO3 #31]')
    
    %*** NOTE: CHANGE um2LD back to match current system
    an_params.um2LD = um2LD_OLD;
    
    fprintf("______End Section 2\n\n")
end

%% Section 3: Best null radial average - FNM2 [21], Simmulations, FTO3 [31]
if Sec2Run(3)
	disp('______Section 3:______')
    
    %################################################################
    %################# Experimental data FNM2 [21]  #################
    %################################################################
    
    %- Get the normalized data (full eta map)
    nmFine  = ['021619_FNM2' filesep '21Don_FineFull2'];
    nrmFine = VFN_An_fitsNormed(nmFine, [1,1,1], an_params);    

    %- Extract eta_s and location
    [eta_sFine, indFine] = VFN_An_getEta_s(nrmFine);
    
    %- Transform from cartesian to polar coordinates, interpolating 
    cent = [indFine(2), indFine(1)];  %[x, y] coords of the null  
    % Calculate maximum radius without clipping
        % max(indFine(1:2)) takes the center value closest to the farther edge
        % this is then subtracted form both edges to see which is closer
        % We also compare to the closer edge and take the min between all of these
    rmax = min(min(size(nrmFine)-max(indFine(1:2)), indFine(1:2))); 
    % use rmax for both radial values to get a 1-1 sampling
    [FinePl,FineR,FineTheta] = polarTransform(nrmFine,cent,rmax,rmax,360,'linear');
    
    %- Calculate radial average
    mns = mean(FinePl,2);
    
    %- Get axes for image (relative to null point):
    [xaxFine, yaxFine] = VFN_An_getAxes(nmFine, indFine, an_params);
    raxFine = xaxFine(indFine(2):indFine(2)+rmax-1);
        
    %-- Display
    %Average for the donut ring with the peak coupling
    mnpkVal = max(mns);
    fprintf('--Actual peak average : Avg = %f\n', mnpkVal)
    fprintf('--Radius of peak average: r = %f\n', raxFine(mns==mnpkVal))
    
    %-- Plot average coupling for varying radii
    %Linear Scale
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
    z = plot(raxFine,mns*100, 'LineWidth',4);
    title('Average Coupling - New F#=5/24', 'Position', [0.5000   16.0859+0.4 0])
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
    set(gca,'fontsize',fontsize2)
    grid on
    
    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmFine);
        export_fig([savefld filesep tmpnm '_RadAvg.png'],'-r300', '-painters');
    end
    
    %################################################################
    %####################### Simulated Data  ########################
    %################################################################
    %- Get the simulated data (full eta map)
    %nmSimm  = [flSimm 'Simulated_CoupMap.fits'];
    %nmSimm  = [flSimm 'Simulated_GaussCoupMap_Df1p0.fits'];
    %nmSimm  = [flSimm 'Simulated_GaussCoupMap_Df0p95.fits'];
    inds = contains(an_params.STRNMS, 'Simulated_');
    nmSimm  = an_params.STRNMS(inds);
    % Remove gaussian simulations
    inds = ~contains(nmSimm, {'GaussCoupMap_Df','Bess'});
    nmSimm = nmSimm(inds);
    
    apRad = 150;
    
    crpVl = 30; %1/2 size of crop window for simulated data
    radPts = 100; % Number of radial points for interpolation
    angPts = 360; % Number of angular points for interpolation
    SimmPl = ones(length(nmSimm),radPts,angPts); %cube of data
    SimmR = ones(length(nmSimm),radPts);    %Radii scaled to lambda/D
    mnsSimm = ones(length(nmSimm),radPts);  %Average around radii
    legLabels = {};     %Legend labels
    
    for i = 1:length(nmSimm)
        % Read file
        tmpSimm = fitsread(nmSimm(i));
        
        % Find null (assuming it's close to origin
        tmpSize  = size(tmpSimm)/2+1;
        sbSz    = 10;    %subwindow size for min search
        tmpSimm2 = tmpSimm(tmpSize(1)-sbSz:tmpSize(1)+sbSz,tmpSize(2)-sbSz:tmpSize+sbSz);
        [mnVal, mnInd] = min(tmpSimm2(:));
        [mnRow, mnCol] = ind2sub(size(tmpSimm2),mnInd);
        mnRow = mnRow+tmpSize(1)-sbSz-1;    %Not sure why you need the 1 but it works
        mnCol = mnCol+tmpSize(2)-sbSz-1;
        
        % Crop to region of interest
        nrmSimm = tmpSimm(mnRow-crpVl:mnRow+crpVl+1,mnCol-crpVl:mnCol+crpVl+1);
        
        %- Transform from cartesian to polar coordinates, interpolating
        cent = [crpVl+1, crpVl];  %[x, y] coords of the null
        % Calculate maximum radius without clipping
        % max(indFine(1:2)) takes the center value closest to the farther edge
        % this is then subtracted form both edges to see which is closer
        % We also compare to the closer edge and take the min between all of these
        rmax = min(min(size(nrmSimm)-max(cent), cent));
        % use rmax for both radial values to get a 1-1 sampling
        [tmpPl,tmpR,~] = polarTransform(nrmSimm,cent,rmax,radPts,angPts,'linear');
        SimmPl(i,:,:) = tmpPl;
        SimmR(i,:) = tmpR;
        %- Calculate radial average
        mnsSimm(i,:) = mean(SimmPl(i,:,:),3);
        %-- Scale coupling (this is to account for miscellaneous losses)
        mnSCVL = 1;%0.705387/0.7806;   %<-- peak coup on PSF in lab / peak coup on PSF in simulator (given known WFE)
        mnsSimm(i,:) = mnsSimm(i,:)*mnSCVL;
        
        %- Get axes for image (from simmulation code):
        lambdaOverD = ((tmpSize(1)-1)*2)/apRad/2;
        SimmR(i,:) = SimmR(i,:)/lambdaOverD;    % Rescale to lambda/D
        
        %- Get just name part of filename
        [~, tmpnm, ~] = fileparts(nmSimm(i));
        tmpnm = erase(tmpnm,[{'Simulated_'},{'CoupMap'},{'_FN5p22'}]);
        for j = 0:9
            tmpnm = strrep(tmpnm, [num2str(j),'p'], [num2str(j),'.']);
        end
        
        %-- Display
        %Average for the donut ring with the peak coupling
        mnpkVal = max(mnsSimm(i,:));
        fprintf('--%s  Values\n', tmpnm)
        fprintf('--    Simulated peak average: Avg = %f\n', mnpkVal)
        tmpR = SimmR(i,:);
        fprintf('--    Radius of peak average  : r = %f\n', tmpR(mnsSimm(i,:)==mnpkVal))
        
        %-- Plot average coupling for varying radii
        if i == 1
            figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
            hold on
        end
        %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
        z = plot(SimmR(i,:),mnsSimm(i,:)*100, '--', 'LineWidth',1.5);
        %axis([0 1.2 0 19]);
        title('Average Coupling - Simmulated', 'Position', [1.25, 20.11+.4 0]);
        xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
        %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
        ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
        set(gca,'fontsize',fontsize2)
        grid on
        
        % Add Legend Label
        legLabels{i} = tmpnm;
    end
    hold off
    %axis([0 2.5 0 25])
    legend(legLabels, 'Interpreter', 'none');
    
    if isSaveFig
        export_fig([savefld filesep 'SimulateCoupling_RadAvg.png'],'-r300', '-painters');
    end
    
	%################################################################
    %################## Old F-number (FTO3 [31]) ####################
    %################################################################
    %*** NOTE: CHANGE um2LD temporarily when analyzing data w/ old F-number
    um2LD_OLD = an_params.um2LD;
    an_params.um2LD = 1/(0.635*11/3.59); 
    
    %-- Get the normalized data (full eta map)
    nmFull  = ['101718_FTO3' filesep '31Don_NullCon'];
    nrmFull = VFN_An_fitsNormed(nmFull, [1,1,1], an_params);

    %-- Extract eta_s and location
    [eta_sFull, indFull] = VFN_An_getEta_s(nrmFull);
    
    %-- Get axes for image (relative to null point):
    [xaxFull, yaxFull] = VFN_An_getAxes(nmFull, indFull, an_params);
    
    %- Transform from cartesian to polar coordinates, interpolating 
    cent = [indFull(2), indFull(1)];  %[x, y] coords of the null  
    % Calculate maximum radius without clipping
        % max(indFine(1:2)) takes the center value closest to the farther edge
        % this is then subtracted form both edges to see which is closer
        % We also compare to the closer edge and take the min between all of these
    rmax = min(min(size(nrmFull)-max(indFull(1:2)), indFull(1:2))); 
    % use rmax for both radial values to get a 1-1 sampling
    [FullPl,FullR,~] = polarTransform(nrmFull,cent,rmax,rmax,360,'linear');
    
    %- Calculate radial average
    mnsFull = mean(FullPl,2);
    
    %- Get axes for image (relative to null point):
    raxFull = xaxFull(indFull(2):indFull(2)+rmax-1);
        
    %-- Display
    %Average for the donut ring with the peak coupling
    mnpkVal = max(mnsFull);
    fprintf('--Actual peak average : Avg = %f\n', mnpkVal)
    fprintf('--Radius of peak average: r = %f\n', raxFull(mnsFull==mnpkVal))
    
    %-- Plot average coupling for varying radii
    %Linear Scale
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
    z = plot(raxFull,mnsFull*100, 'LineWidth',4);
    title('Average Coupling - Old F#=3.06', 'Position', [0.75000   12.04+0.4 0])
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
    set(gca,'fontsize',fontsize2)
    grid on
    
    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmFull);
        export_fig([savefld filesep tmpnm '_RadAvg.png'],'-r300', '-painters');
    end
    
    %*** NOTE: CHANGE um2LD back to match current system
    an_params.um2LD = um2LD_OLD;
    
    %################################################################
    %####################### Plot for All  ##########################
    %################################################################
    
    %-- Plot both average coupling for varying radii
    %Linear Scale
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
    z2 = plot(raxFine,mns*100, 'LineWidth',4);
    hold on
    % Plot old F#=3.06 data
    z3 = plot(raxFull,mnsFull*100, 'LineWidth',4);
    % Iterate through simmulation data
    for i = 1:length(nmSimm)
        z = plot(SimmR(i,:),mnsSimm(i,:)*100, '--', 'LineWidth',1.5);
    end
    hold off
    %axis([0 1.2 0 19]);
    %axis([0 2 0 25])
    legLabels = [{'Lab F#=5.24'} {'Lab F#=3.06'} legLabels];
    legend(legLabels, 'Interpreter', 'none');
    title('Average Coupling - All', 'Position', [1.25, 20.11+0.4 0]);
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
    set(gca,'fontsize',fontsize2)
    grid on
    
    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmSimm);
        export_fig([savefld filesep tmpnm '_RadAvg.png'],'-r300', '-painters');
    end
    
    %{
    %calculate polar coords for cartesian
    [theta, r] = fitsPol(nrmFine, indFine(2), indFine(1));
    %preallocate loop values
    j   = 0 : 0.02 : 1;          %radii over which to analyze
    mns = nan(length(j),1); %array to hold means. nan to check for errors
    tol = 0.01;     %Tolerance for radii 
    %iterate through radii and calculate average
    for k = 1:length(j)
        i       = j(k);
        %uncomment below to plot circles
        fitspl  = nrmFine;
        fitspl((i-tol<=r)&(r<=i+tol)) = 1000;
        fitsDisp(fitspl);
        mns(k) = mean(nrmFine((i-tol<=r)&(r<=i+tol)));
    end
    
    %-- Display
    %Average for the donut ring with the peak coupling
    mnpkVal = max(mns);
    fprintf('--Radius of peak average: r = %f\n', j(mns==mnpkVal))
    
    %-- Plot average coupling for varying radii
    %Linear Scale
    figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
    %plot(0:imLim/(length(j)):imLim,((mns-avgBkgd)/(mnpkVal-avgBkgd)),'-o', 'Color', [0.8500    0.3250    0.0980])
    plot(0:xaxFine(end)/(length(j)-1):xaxFine(end),mns*100,'-o', 'Color', [0.8500    0.3250    0.0980])
    title('Average Coupling')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    %ylabel('Coupling Normalized by Peak', 'FontSize', fontsize1)
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1)
    set(gca,'fontsize',fontsize2)
    set(findall(gca, 'Type', 'Line'),'LineWidth', 2);
    grid on
    %export_fig([savefld filesep 'ActualCouplingLine.png'],'-r300', '-painters');
    %}
    fprintf('______End Section 3\n\n')
    %{
    disp("______Section 2:______")
    
    %################################################################
    %################## FOV1 (null) frame (#5)  #####################
    %################################################################
    %- Get the normalized data (full eta map)
    nmFOV1  = ['102018_FTO5' filesep '5Don_NullConf2'];
    nrmFOV1 = fitsNormed(nmFOV1);

    %- Extract eta_s and location
    [eta_sFOV1, indFOV1] = getEta_s(nrmFOV1);
    
    %- Print eta_s info
    fprintf('\nEta_s in FOV1 (null) Frame:    %f\n',eta_sFOV1)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFOV1(2), indFOV1(1), indFOV1(4), indFOV1(5), indFOV1(6));
    end

    %- Calculate eta_p
    % NOTE: Eta_p nonsensical for a null frame; FOV doesn't include planet
    
    %- Get axes for image (relative to null point):
    [xaxFOV1, yaxFOV1] = getAxes(nmFOV1, indFOV1);
        
    %- Show Null frame
    fighFOV1 = fitsDisp(nrmFOV1, xaxFOV1, yaxFOV1);
    title('FOV1 Frame - Null: [FTO5 #5]')
    
    
    %################################################################
    %################## FOV2 frame (#6)  ############################
    %################################################################
    %- Get the normalized data (full eta map)
    nmFOV2  = ['102018_FTO5' filesep '6Don_NullConf3'];
    nrmFOV2 = fitsNormed(nmFOV2);

    %- Extract eta_s and location
    [eta_sFOV2, indFOV2] = getEta_s(nrmFOV2);
    
    %- Print eta_s info
    fprintf('\nEta_s in FOV2 (medm) Frame:    %f\n',eta_sFOV2)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFOV2(2), indFOV2(1), indFOV2(4), indFOV2(5), indFOV2(6));
    end
    
    %- Calculate eta_p
    % NOTE: Eta_p nonsensical for a null frame; FOV doesn't include planet
    
    %- Get axes for image (relative to null point):
    [xaxFOV2, yaxFOV2] = getAxes(nmFOV2, indFOV2);
        
    %- Show Null frame
    fighFOV2 = fitsDisp(nrmFOV2, xaxFOV2, yaxFOV2);
    title('FOV2 Frame: [FTO5 #6]')
    
    
    %################################################################
    %################## FOV3 (peak) frame (#7)  #####################
    %################################################################
    %- Get the normalized data (full eta map)
    nmFOV3  = ['102018_FTO5' filesep '7Don_NullConf4'];
    nrmFOV3 = fitsNormed(nmFOV3);

    %- Extract eta_s and location
    [eta_sFOV3, indFOV3] = getEta_s(nrmFOV3);
    
    %- Print eta_s info
    fprintf('\nEta_s in FOV3 (peak) Frame:    %f\n',eta_sFOV3)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFOV3(2), indFOV3(1), indFOV3(4), indFOV3(5), indFOV3(6));
    end

    %- Calculate eta_p
    [eta_pFOV3, indPFOV3] = getEta_p(nrmFOV3);
    fprintf('\nEta_p in FOV3 (peak) Frame:    %f\n',eta_pFOV3);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPFOV3(2), indPFOV3(1), indPFOV3(4), indPFOV3(5), indPFOV3(6));
    end

    %- Get axes for image (relative to null point):
    [xaxFOV3, yaxFOV3] = getAxes(nmFOV3, indFOV3);
        
    %- Show Null frame
    fighFOV3 = fitsDisp(nrmFOV3, xaxFOV3, yaxFOV3);
    title('FOV3 Frame - peak: [FTO5 #7]')
    
    
    %################################################################
    %################## FOV4 (peak) frame (#8)  #####################
    %################################################################
    %- Get the normalized data (full eta map)
    nmFOV4  = ['102018_FTO5' filesep '8Don_NullConf5'];
    nrmFOV4 = fitsNormed(nmFOV4);

    %- Extract eta_s and location
    [eta_sFOV4, indFOV4] = getEta_s(nrmFOV4);
    
    %- Print eta_s info
    fprintf('\nEta_s in FOV4 (full) Frame:    %f\n',eta_sFOV4)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFOV4(2), indFOV4(1), indFOV4(4), indFOV4(5), indFOV4(6));
    end

    %- Calculate eta_p
    [eta_pFOV4, indPFOV4] = getEta_p(nrmFOV4);
    fprintf('\nEta_p in FOV4 (full) Frame:    %f\n',eta_pFOV4);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPFOV4(2), indPFOV4(1), indPFOV4(4), indPFOV4(5), indPFOV4(6));
    end

    %- Get axes for image (relative to null point):
    [xaxFOV4, yaxFOV4] = getAxes(nmFOV4, indFOV4);
        
    %- Show Null frame
    fighFOV4 = fitsDisp(nrmFOV4, xaxFOV4, yaxFOV4);
    title('FOV4 Frame - Full: [FTO5 #8]')
        
    fprintf("______End Section 2\n\n")
    %}
end

%% Section 4: Single, Asymmetric, Full FOV scan - FTO4 [11]
%  --- This is a good one for the line scan but only goes to ~1.5 L/D
if Sec2Run(4)
	disp("______Section 4:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['101818_FTO4' filesep '11Don_WideFOVFineScan'];
    nrmFine = VFN_An_fitsNormed(nmFine, [1,1,1], an_params);

    %- Extract eta_s and location
    [eta_sFine, indFine] = VFN_An_getEta_s(nrmFine);
    
    %- Print eta_s info
    fprintf('\nEta_s in Asym (fine) Frame:    %f\n',eta_sFine)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indFine(2), indFine(1), indFine(4), indFine(5), indFine(6));
    end

    %- Calculate eta_p
    [eta_pFine, indPFine] = VFN_An_getEta_p(nrmFine);
    fprintf('\nEta_p in Asym (fine) Frame:    %f\n',eta_pFine);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPFine(2), indPFine(1), indPFine(4), indPFine(5), indPFine(6));
    end
    
    %- Get axes for image (relative to null point):
    [xaxFine, yaxFine] = VFN_An_getAxes(nmFine, indFine, an_params);
        
    %- Show Null frame
    fighAsym = VFN_An_fitsDisp(nrmFine, xaxFine, yaxFine);
    title('Asym Frame - Fine: [FTO4 #11]')
       
    fprintf("______End Section 4\n\n")
end 
    
%% Section 5: Donut scan with vortex centered for symmetry - MNO1 [10]
%  --- The null/planet coupling are bad but the symmetry is good
if Sec2Run(5)
	disp("______Section 5:______")
    
    %- Get the normalized data (full eta map)
    nmSymm  = ['102618_MNO1' filesep '10Don_Symmetry8'];
    nrmSymm = VFN_An_fitsNormed(nmSymm, [1,1,1], an_params);

    %- Extract eta_s and location
    [eta_sSymm, indSymm] = VFN_An_getEta_s(nrmSymm, 4);
    
    %- Print eta_s info
    fprintf('\nEta_s in Symm (symm) Frame:    %f\n',eta_sSymm)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indSymm(2), indSymm(1), indSymm(4), indSymm(5), indSymm(6));
    end

    %- Calculate eta_p
    [eta_pSymm, indPSymm] = VFN_An_getEta_p(nrmSymm);
    fprintf('\nEta_p in Symm (symm) Frame:    %f\n',eta_pSymm);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPSymm(2), indPSymm(1), indPSymm(4), indPSymm(5), indPSymm(6));
    end
    
    %- Get axes for image (relative to null point):
    [xaxSymm, yaxSymm] = VFn_An_getAxes(nmSymm, indSymm, an_params);
        
    %- Show Null frame
    fighSymm = VFN_An_fitsDisp(nrmSymm, xaxSymm, yaxSymm);
    title('Symm Frame - symm: [MNO1 #10]')
       
    fprintf("______End Section 5\n\n")
end 

%% Section 6: Peak coupling on System PSF w/ vort out - FNM2 [23]
if Sec2Run(6)
	disp("______Section 6:______")
    
    %- Get the normalized data (full eta map)
    nmPeak  = ['021619_FNM2' filesep '23PSF_PeakNewFNM2'];
    nrmPeak = VFN_An_fitsNormed(nmPeak, [1,1,1], an_params);
    
    %- Calculate peak Coupling -- SAVED AS eta_p
      % the eta_p function will have the same effect here even though this
      % is actually not a donut scan.
    [eta_pPeak, indPPeak] = VFN_An_getEta_p(nrmPeak);
    fprintf('\nEta_p in Peak (psf) Frame:    %f\n',eta_pPeak);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPPeak(2), indPPeak(1), indPPeak(4), indPPeak(5), indPPeak(6));
    end
    
    %- Get axes for image (relative to peak point):
    [xaxPeak, yaxPeak] = VFN_An_getAxes(nmPeak, indPPeak, an_params);
        
    %- Show Null frame
    fighPeak = VFN_An_fitsDisp(nrmPeak, xaxPeak, yaxPeak);
    title('Peak Frame - psf: [FNM2 #23]')
       
    fprintf("______End Section 6\n\n")
end 

%% Section 7: Symmetric, Circular coupling maps - MNO1 [5]
%  NOTE: extracts specific frames from the scanned cube
if Sec2Run(7)
	disp("______Section 7:______")
    
    
    %################################################################
    %################## Crc1 (circ) frame (1,4,2)  ##################
    %################################################################
    
    %- Get the normalized data (full eta map)
    nmCrc1  = ['102618_MNO1' filesep '5Don_Symmetry3'];
    nrmCrc1 = VFN_An_fitsNormed(nmCrc1,[1, 4, 2], an_params);
    
    %- Extract eta_s and location
    [eta_sCrc1, indCrc1] = VFN_An_getEta_s(nrmCrc1, 4);
    
    %- Print eta_s info
    fprintf('\nEta_s in Crc1 (circ) Frame:    %f\n',eta_sCrc1)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indCrc1(2), indCrc1(1), indCrc1(4), indCrc1(5), indCrc1(6));
    end
    
    %- Calculate peak Coupling
    [eta_pCrc1, indPCrc1] = VFN_An_getEta_p(nrmCrc1);
    fprintf('\nEta_p in Crc1 (circ) Frame:    %f\n',eta_pCrc1);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPCrc1(2), indPCrc1(1), indPCrc1(4), indPCrc1(5), indPCrc1(6));
    end
    
    %- Get axes for image (relative to null point):
    [xaxCrc1, yaxCrc1] = VFN_An_getAxes(nmCrc1, indCrc1, an_params);
        
    %- Show Null frame
    fighCrc1 = VFN_An_fitsDisp(nrmCrc1, xaxCrc1, yaxCrc1);
    title('Crc1 Frame - circ: [MNO1 #15,1]')
    
    %################################################################
    %################## Crc2 (circ) frame (1,4,3)  ##################
    %################################################################
    
    %- Get the normalized data (full eta map)
    nmCrc2  = ['102618_MNO1' filesep '5Don_Symmetry3'];
    nrmCrc2 = VFN_An_fitsNormed(nmCrc2,[1, 4, 3], an_params);
    
    %- Extract eta_s and location
    [eta_sCrc2, indCrc2] = VFN_An_getEta_s(nrmCrc2, 4);
    
    %- Print eta_s info
    fprintf('\nEta_s in Crc2 (circ) Frame:    %f\n',eta_sCrc2)
    if isPrintLoc
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indCrc2(2), indCrc2(1), indCrc2(4), indCrc2(5), indCrc2(6));
    end
    
    %- Calculate peak Coupling
    [eta_pCrc2, indPCrc2] = VFN_An_getEta_p(nrmCrc2);
    fprintf('\nEta_p in Crc2 (circ) Frame:    %f\n',eta_pCrc2);
    if isPrintLoc
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indPCrc2(2), indPCrc2(1), indPCrc2(4), indPCrc2(5), indPCrc2(6));
    end
    
    %- Get axes for image (relative to null point):
    [xaxCrc2, yaxCrc2] = VFN_An_getAxes(nmCrc2, indCrc2, an_params);
        
    %- Show Null frame
    fighCrc2 = VFN_An_fitsDisp(nrmCrc2, xaxCrc2, yaxCrc2);
    title('Crc2 Frame - circ: [MNO1 #15,2]')
       
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