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
addpath('C:\Users\Thomas\Desktop\VFN-Simulations-master\VFN-Simulations-master\VFNlib')


% Folder with data
%fldnm = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling';
fldnm = 'C:\Users\Thomas\Desktop\071019_VAR1';

%flSimm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\ProcessedData\';
% Folder for saving results
%savefld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\071019_VAR1test';
savefld = 'C:\Users\Thomas\Desktop\VFN_Data';

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
    true;          %Section 1  ** Donut Sensor Matrix
    true;          %Section 2  ** Donut Sensor Matrix (Log scale)
    false;          %Section 3  ** etaS linear plot
    false;          %Section 4  ** etaS ylog plot
    false;          %Section 5  ** etaS xylog plot
    false;          %Section 6  ** etaP linear plot
    false;          %Section 7  ** etaP ylog plot
    false;          %Section 8  ** etaP xylog plot
    false;          %Section END  -- REFERENCE ONLY; DO NOT ENABLE
    ];
% Flag for fitting a line to plots
isLinFit = false;

% Flag for using radial average instead of direct min and max values
isRadAvg = true;

% Flag for using a fixed fiber position (fixed null and radAvg position)
% for all bandwidths
isFixed = true;
    
% Flag to print locations
isPrintLoc = false;

% Flag to mark whether (.mat) data should be saved
isSave = false;

% Flag to mark whether figures should be saved
isSaveFig = true;

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
    nmFine  = ['071019_VAR1' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));

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
        %[eta_pFine(i), indPFine(i)] = VFN_An_getEta_p(nrmFine(i,:,:));
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
    nmFine  = ['071019_VAR1' filesep filetag];
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));

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
        %[eta_pFine(i), indPFine(i)] = VFN_An_getEta_p(nrmFine(i,:,:));
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
    nmFine  = ['071019_VAR1' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));

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
        eta_sFine = eta_sFine .* 1000;
        scatter(xax,eta_sFine,50, 'filled')
        xlabel('Bandwidth  [\Delta\lambda/\lambda_0]')
        ylabel('\eta_{s} [\times 10^{-3}]')
        if max(eta_sFine)*1.2>1
            ylim([0 max(eta_sFine)*1.2]);
        else
            ylim([0 1]);
        end
        set(gca,'fontsize',fontsize2)
        grid on;
        range = [xlim;ylim];
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
        if isLinFit
            P = polyfit(xax,eta_sFine,1);
            yfit = P(1)*xax+P(2);
            hold on;
            plot(xax,yfit,'r-.');
            text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        end
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
    nmFine  = ['071019_VAR1' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));
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
    scatter(xax,eta_sFine,50, 'filled')
    xlabel('Bandwidth  [\Delta\lambda/\lambda_0]')
    ylabel('\eta_{s} [log]')
    set(gca,'fontsize',fontsize2)
    range = [xlim;ylim];
    title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
    if isLinFit
        P = polyfit(xax,eta_sFine,1);
        yfit = P(1)*xax+P(2);
        hold on;
        plot(xax,yfit,'r-.');
        text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
    end

    if isSaveFig
        [~, tmpnm, ~] = fileparts(nmFine);
        export_fig([savefld filesep tmpnm '_An_BWGraph_Log_etaS_CoupMap.png'],'-r300', '-painters');
    end
    
    fprintf("______End Section 4\n\n")
end


%% Section 5: Eta S xyLog plot
%  --- The null/planet coupling are bad but the symmetry is good
if Sec2Run(5)
	disp("______Section 5:______")
    
    %- Get the normalized data (full eta map)
    nmFine  = ['071019_VAR1' filesep filetag];
    
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);

        %- Extract eta_s and location
        [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));

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
        xax = log10(xax);
        eta_sFine = log10(eta_sFine);
        scatter(xax,eta_sFine,50, 'filled')
        xlabel(['Bandwidth about ' num2str(Wvls(j)) 'nm [log(\Delta\lambda/\lambda_0)]'])
        ylabel('\eta_{s} [log]')
        set(gca,'fontsize',fontsize2)
        range = [xlim;ylim];
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
        if isLinFit
            P = polyfit(xax,eta_sFine,1);
            yfit = P(1)*xax+P(2);
            hold on;
            plot(xax,yfit,'r-.');
            text(xax(1),eta_sFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        end
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_xyLog_etaS_CoupMap.png'],'-r300', '-painters');
        end
    end
       
    fprintf("______End Section 5\n\n")
end 
%% Setup for Eta P plots
if Sec2Run(6) || Sec2Run(7) || Sec2Run(7)
    %- Get the normalized data (full eta map)
    nmFine  = ['071019_VAR1' filesep filetag];
    
    eta_pFine = nan(length(BWs),1);
    
    for i = 1:length(BWs)
        nrmFine(i,:,:) = VFN_An_fitsNormed(nmFine,[1,i,1],an_params);
        if isRadAvg && ~isFixed
            %- Extract eta_s and location
            [eta_sFine(i), indFine(i,:)] = VFN_An_getEta_s(squeeze(nrmFine(i,:,:)));

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
            %-- Get axes for image (relative to null point):
            [xax2a, yax2a]  = VFN_An_getAxes(nmFine, squeeze(indFine(i,:)), an_params);
            %-- Calculate average radial profile
            [radvg2a,rvec2a] = VFN_An_radAverage(squeeze(nrmFine(i,:,:)), [indFine(i,2),indFine(i,1)], length(nrmFine(1,1,:)));
            
            %-- Get radial axis 
            % Use xax to get L/D per 'pixel' value then rescale rvec by this value
            rax2a = (xax2a(2)-xax2a(1))*rvec2a;

            %-- Print Peak
            %Average for the donut ring with the peak coupling
            eta_pFine(i) = max(radvg2a);
            %figure; scatter(rax2a, radvg2a); ylim([0, 0.062]);
        elseif ~isRadAvg
            %- Calculate eta_p
            [eta_pFine(i), indPFine(i,:)] = VFN_An_getEta_p(nrmFine(i,:,:));
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
    end
    
    if isRadAvg && isFixed
        % Average 
        avgNrmFine = squeeze(mean(nrmFine,1));
        %- Extract eta_s and location
        [eta_sFix, ind_SFix] = VFN_An_getEta_s(avgNrmFine);
        
        %- Extract position of max coupling
        
        %-- Get axes for image (relative to null point):
        [xax2a, yax2a]  = VFN_An_getAxes(nmFine, ind_SFix, an_params);
        %-- Calculate average radial profile
        [radvg2a,rvec2a] = VFN_An_radAverage(avgNrmFine, [ind_SFix(2),ind_SFix(1)], length(avgNrmFine(1,:)));

        %-- Get radial axis 
        % Use xax to get L/D per 'pixel' value then rescale rvec by this value
        rax2a = (xax2a(2)-xax2a(1))*rvec2a;

        %-- Print Peak
        %Average for the donut ring with the peak coupling
        [~, distR] = max(radvg2a);
        
        for i = 1:length(BWs)
            %-- Get axes for image (relative to null point):
            [xax2a, yax2a]  = VFN_An_getAxes(nmFine, ind_SFix, an_params);
            %-- Calculate average radial profile
            [radvg2a,rvec2a] = VFN_An_radAverage(squeeze(nrmFine(i,:,:)), [ind_SFix(2),ind_SFix(1)], length(nrmFine(1,1,:)));

            %-- Get radial axis 
            % Use xax to get L/D per 'pixel' value then rescale rvec by this value
            rax2a = (xax2a(2)-xax2a(1))*rvec2a;

            %-- Print Peak
            %Average for the donut ring with the peak coupling
            eta_pFine(i) = radvg2a(distR);
            %figure; scatter(rax2a, radvg2a); ylim([0, 0.062]);
        end
    end
end
%% Section 6: Eta P linear plot
if Sec2Run(6)
	disp("______Section 6:______")
    
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,eta_pFine*100,50,'red', 'filled')
        xlabel('Bandwidth  [\Delta\lambda/\lambda_0]')
        ylabel('\eta_{p} [%]')
        if 1.2*max(eta_pFine)*100>8
            ylim([0 1.2*max(eta_pFine)*100])
        else
            ylim([0 8])
        end
        range = [xlim;ylim];
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
        set(gca,'fontsize',fontsize2)
        grid on;
        
        if isLinFit
            P = polyfit(xax,eta_pFine,1);
            yfit = P(1)*xax+P(2);
            hold on;
            plot(xax,yfit,'r-.');
            text(xax(1),eta_pFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        end
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
    
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(xax,log10(eta_pFine),50,'red', 'filled')
        xlabel('Bandwidth  [\Delta\lambda/\lambda_0]')
        ylabel('\eta_{p} [log]')
        set(gca,'fontsize',fontsize2)
        range = [xlim;ylim];
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
        if isLinFit
            P = polyfit(xax,log10(eta_pFine),1);
            yfit = P(1)*xax+P(2);
            hold on;
            plot(xax,yfit,'r-.');
            text(xax(1),eta_pFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        end
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_Log_etaP_CoupMap.png'],'-r300', '-painters');
        end
    end
    
       
    fprintf("______End Section 7\n\n")
end 
%% Section 8: Eta P xylog plot
%  NOTE: extracts specific frames from the scanned cube
if Sec2Run(8)
	disp("______Section 8:______")
    
    for j = 1:length(Wvls)
        %- Show Null frame
        figh = figure('color','white','units', 'inches', 'Position', [0 0 7 6]);
        xax = BWs ./ Wvls(1);
        scatter(log10(xax),log10(eta_pFine),50,'red', 'filled')
        xlabel(['Bandwidth about ' num2str(Wvls(j)) 'nm [log(\Delta\lambda/\lambda_0)]'])
        ylabel('\eta_{p} [log]')
        set(gca,'fontsize',fontsize2)
        range = [xlim;ylim];
        title(['Bandwidth Scan about ' num2str(Wvls(1)) 'nm'], 'VerticalAlignment', 'bottom','HorizontalAlignment', 'left', 'Position', [range(1,1) range(2,2)+.05*(range(2,2)-range(2,1)) 0], 'fontsize', fontsize2)
        if isLinFit
            P = polyfit(log10(xax),log10(eta_pFine),1);
            yfit = P(1)*log10(xax)+P(2);
            hold on;
            plot(log10(xax),yfit,'r-.');
            text(log10(xax(1)),eta_pFine(end),['y = ' num2str(P(1)) 'x + ' num2str(P(2))])
        end
        if isSaveFig
            [~, tmpnm, ~] = fileparts(nmFine);
            export_fig([savefld filesep tmpnm '_An_BWGraph_xyLog_etaP_CoupMap.png'],'-r300', '-painters');
        end
    end
    
       
    fprintf("______End Section 8\n\n")
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