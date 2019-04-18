% New, updated plotting script for the VFN MonoWrapUp data. It uses the new libraries 
% and has been cleaned of many of the unecessary/extraneous sections present in the 
% SAMPLE script. 

%% Start Fresh
close all;
clear;

%% Define Parameters and Load Data
%-- Path to VFN-Simulation Library (Polar Transform function)
addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib')


%-- Path to VFN-Lab Library (Analysis functions)
addpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\AnalysisLib')


%-- Data Folders
% Lab
fldnm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\PupilVFN\CouplingCurves';
% Simulations
flSimm = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\NewAnalysis\';
%-- Ouput Folder (where data is saved)
savefld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\NewAnalysis\';


%-- Define analysis parameters
% Scaling factor between Femto gain settings
an_params.gnFactor = 9.97;

% Final lens throughput (for red pm normalization)
an_params.lensTHPT = 0.9973;

% Conversion: Red thorlabs PM [uW] to fmto @10^6 [V]
an_params.uw2fmto  = 0.6318;        %P_femto[V]/P_red[uW]
an_params.rdOfmto  = 1/an_params.uw2fmto;    %ratio as defined by Nem (red[uW]/fmto[V])

% Conversion: microns to L/D in focal plane [L/D / micron]
an_params.um2LD = 1/(0.635*10.96/2.1);  %(1/lambda*F#)
                            % f=10.96 taken from Thorlabs A397 datasheet

% Conversion: piezo [V] to microns
an_params.V2um = 1/8.72;    %<-- Value from PZT characterization (spec = 1/10)

% Fresnel Reflections at fiber surfaces
an_params.FR = 0.9654;  % (1-loss)
% Propagation losses through fiber (fraction per meter)
an_params.PL = 0.9966;  % (1-loss)



%-- Read data
% Read all .fits files in the data directories
nmWpath = { [fldnm '\021619_FNM2\*.fits'];
            [flSimm '\*.fits'];
            [fldnm '\*_FTO3\*.fits']; 
            [fldnm '\*_FTO4\*.fits']};
an_params.STRNMS = VFN_An_getSTRNMS(nmWpath);

%-- Sections to Run 
Sec2Run = [
    true;       %Section 1 - FNM2 [23] Peak Coupling on PSF
    false;        %Section 2 - FNM2 [21] Null Coupling on Donut
    false;       %Section 3 - MNO1 [15] Peak Coupling on OLD FNUM PSF
    false;        %Section 4 - FTO3 [31] Null Coupling on OLD FNUM Donut
    ];

%-- Flags for save/print options
isPrintLoc  = false;    % Print null/planet locations
isSave      = false;    % Save (.mat) data
isSaveFig   = false;    % Save figures


% Default fonts for saved figures
fontsize1 = 14;
fontsize2 = 14;

%% Section 1: PSF scan for peak coupling - FNM2 [23]
seci    = 1;
if Sec2Run(seci)
    fprintf('\n---- Section %i ----\n',seci)
    
    
    %################################################################
    %################## Basic analysis (max in frame) 
    %################################################################
    fprintf('\n-- Basic analysis\n')
    
    %-- Get normalized data
    nm1a    = ['021619_FNM2' filesep '23PSF_PeakNewFNM2'];
    nrm1a   = VFN_An_fitsNormed(nm1a, [1,1,1], an_params);
    
    %-- Calculate peak coupling
    [eta_p1a, indP1a] = VFN_An_getEta_p(nrm1a);
    fprintf('\nPeak coupling on PSF [FNM2 #23]:    %f\n',eta_p1a);
    if isPrintLoc
        fprintf([' Peak Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 indP1a(2), indP1a(1), indP1a(4), indP1a(5), indP1a(6));
    end
    
    %-- Get axes for image (relative to peak point):
    [xax1a, yax1a]  = VFN_An_getAxes(nm1a, indP1a, an_params);
    
    %-- Show Frame
    figh1a = VFN_An_fitsDisp(nrm1a, xax1a, yax1a);
    title('PSF Coupling - [FNM2 #23]')
    
    
    %################################################################
    %################## Radial Average  
    %################################################################
    fprintf('\n-- Radial Average\n')
    %-- Calculate average radial profile
    [radvg1a,rvec1a] = VFN_An_radAverage(nrm1a, [indP1a(2), indP1a(1)]);
    
    %-- Get radial axis 
    % Use xax to get L/D per 'pixel' value then rescale rvec by this value
    rax1a = (xax1a(2)-xax1a(1))*rvec1a;
    
    %-- Display radial profile
    figh1a_rad = figure('color', 'white');
    z = plot(rax1a,radvg1a*100, 'LineWidth', 4);
    set(gca,'fontsize',fontsize2)
    title('Average radial coupling - [FNM2 #23]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    
    
	%################################################################
    %################## Compare to ideal simulation with WFE
    %################################################################
    fprintf('\n-- Analyze ideal simulation\n')
    %-- Get simulated data 
    % Get simulation filename
    inds  = contains(an_params.STRNMS, 'Sim_Gauss_PSF_DF1.4_FN5.22_WFEOn');
    nmSim1a = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim1a, {'ignore'});   
    nmSim1a = nmSim1a(inds);
        
    % NOTE: THE VORTEX POSITION FOR THESE SIMULATIONS WAS ACTUALLY INCORRECT
            % THE VORTEX WAS AT ~1.6147*apRad instead of ~5.3apRad.
            % However, this appears not to have a big effect on the PSF coupling
            % Results with the actual vortex position have 'VortPosFix' prefix.
    
    %-- Load simulation data
    sim1a = fitsread(nmSim1a);    
    
    %-- Find peak
    [sim1a_etap, sim1a_ind] = VFN_An_getEta_p(sim1a);
    fprintf('\nPeak coupling on PSF [ideal w/ WFE]:  %f\n',sim1a_etap);
       
    %-- Get axes in L/D    
    [sim1a_xax, sim1a_yax, sim1a_lambdaOverD] = ...
        VFN_An_simGetAxes(nmSim1a, sim1a_ind, an_params);
    
    %-- Display 2D map
    figh1a = VFN_An_fitsDisp(sim1a, sim1a_xax, sim1a_yax);
    axis([-1 1 -1 1]);
    title('PSF Coupling - [Sim ideal F#, WFE On]')    
    
    %-- Calculate average radial profile
    [sim1a_radvg,sim1a_rvec] = ...
        VFN_An_radAverage(sim1a, [sim1a_ind(2), sim1a_ind(1)]);
    
    %-- Get radial axis by rescaling rvec
    sim1a_rax = sim1a_rvec/sim1a_lambdaOverD;
    
    %-- Display radial profile
    % Simulation data
    fighsim1a_rad = figure('color', 'white');
    z = plot(sim1a_rax,sim1a_radvg*100, 'LineWidth', 4);
    xlim([0 1])
    set(gca,'fontsize',fontsize2)
    title('Average rad coup - [Sim ideal F#, WFE On]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    % Superimpose lab data
    hold on
    plot(rax1a, radvg1a*100, 'LineWidth', 4);
    legend({'Sim','Lab'}, 'FontSize', fontsize1);
    
    %################################################################
    %################## Compare to varying Df with WFE
    %################################################################
    fprintf('\n-- Analyze varying Df simulation\n')
    
    %-- Get simulated data 
    % Get simulation filenames
    inds  = contains(an_params.STRNMS, 'Sim_Gauss_PSF');
    nmSimDf = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSimDf, {'DF1.8','DF2.1','VortPosFix'});   
    nmSimDf = nmSimDf(inds);
    
    % NOTE: THE VORTEX POSITION FOR THESE SIMULATIONS WAS ACTUALLY INCORRECT
            % THE VORTEX WAS AT ~1.6147*apRad instead of ~5.3apRad.
            % However, this appears not to have a big effect on the PSF coupling
            % Results with the actual vortex position have 'VortPosFix' prefix.
    
    %-- Iterate through reading and processing simulated data
    simrpts = 4000;     % 4000 since ~2^13/2 radial samples
    % Preallocate arrays for results
    sim1b_radvg = nan(length(nmSimDf), simrpts); 
    sim1b_rax = nan(length(nmSimDf),simrpts);
    % Create figure for plotting
    fighsim1b = figure('color', 'white');
    hold on;
    for i = 1:length(nmSimDf)
        
        %-- Load simulation data
        sim1b = fitsread(nmSimDf(i));
        
        %-- Find peak
        [sim1b_etap, sim1b_ind] = VFN_An_getEta_p(sim1b);
        
        %-- Get just name part of filename
        [~, tmpnm, ~] = fileparts(nmSimDf(i));
        tmpnm = erase(tmpnm,[{'Sim_'},{'_FN5.22'},{'Gauss_PSF_'},{'WFEOn'}]);
        
        %-- Print peak coupling value
        fprintf('\n-%s Values\n', tmpnm);
        fprintf(' Peak coupling on PSF [Sim w/ WFE]:  %f\n',sim1b_etap);
        
        %-- Get lambda/D value
        [~, ~, lOverD] = ...
            VFN_An_simGetAxes(nmSimDf(i), sim1b_ind, an_params);
        
        %-- Calculate average radial profile
        [sim1b_radvg(i,:),tmpRvec] = ...
            VFN_An_radAverage(sim1b, [sim1b_ind(2), sim1b_ind(1)], simrpts);
        %sim1b_radvg(i,:) = tmpRadvg;
        
        %-- Get radial axis by rescaling rvec
        sim1b_rax(i,:) = tmpRvec/lOverD;
               
        %-- Plot results
        if contains(tmpnm, {'DF1.1', 'DF1.4', 'DF1.6'})
            linstyle = ':';
        else
            linstyle = '--';
        end
        plot(sim1b_rax(i,:),sim1b_radvg(i,:)*100, linstyle, 'LineWidth', 1.5);
        
        %-- Add legend label
        legLabels{i} = tmpnm;
    end
    
    %-- Finish formatting figure
    % Superimpose lab data
    plot(rax1a, radvg1a*100, 'LineWidth', 3);
    hold off
    set(gca,'fontsize',fontsize2)
    lablabel = sprintf('Lab: PZ = %0.2f', 1/an_params.V2um);
    legLabels = [legLabels {lablabel}];
    % Add legend
    legend(legLabels, 'Interpreter', 'none', 'FontSize', fontsize1);
    % Add title and axes
    title(sprintf('Peak Coupling - Simulated vs Lab PZGain = %0.2f', 1/an_params.V2um));
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1);
    xlim([0 1])
    
    if false
        svnm = sprintf('PSFCoup_SimVsLab_PZ%0.2f_w6.5FibMod.png', 1/an_params.V2um);
        export_fig([savefld filesep svnm],'-r300', '-painters');
    end
    
    fprintf('\n---- End Section %i \n\n',seci)
end

%% Section 2: Donut scan for Null - FNM2 [21]
seci    = 2;
if Sec2Run(seci)
    fprintf('\n---- Section %i ----\n',seci)
    
    
    %################################################################
    %################## Basic analysis (max in frame) 
    %################################################################
    fprintf('\n-- Basic analysis\n')
    
    %-- Get normalized data
    nm2a    = ['021619_FNM2' filesep '21Don_FineFull2'];
    nrm2a   = VFN_An_fitsNormed(nm2a, [1,1,1], an_params);
    
    %-- Calculate null coupling
    [eta_2a, ind2a] = VFN_An_getEta_s(nrm2a);
    fprintf('\nNull coupling on Donut [FNM2 #21]:    %f\n',eta_2a);
    if isPrintLoc
        fprintf([' Null Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 ind2a(2), ind2a(1), ind2a(4), ind2a(5), ind2a(6));
    end
    
    %-- Get axes for image (relative to null point):
    [xax2a, yax2a]  = VFN_An_getAxes(nm2a, ind2a, an_params);
    
    %-- Show Frame
    figh2a = VFN_An_fitsDisp(nrm2a, xax2a, yax2a);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Donut Coupling - [FNM2 #21]')
    
    
    %################################################################
    %################## Radial Average  
    %################################################################
    fprintf('\n-- Radial Average\n')
    %-- Calculate average radial profile
    [radvg2a,rvec2a] = VFN_An_radAverage(nrm2a, [ind2a(2), ind2a(1)]);
    
    %-- Get radial axis 
    % Use xax to get L/D per 'pixel' value then rescale rvec by this value
    rax2a = (xax2a(2)-xax2a(1))*rvec2a;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    radvgPK2a = max(radvg2a);
    fprintf('Peak average [FNM2 #21]: = %f\n', radvgPK2a)
    fprintf(' Radius of peak average: = %f\n', rax2a(radvg2a==radvgPK2a))
    
    %-- Display radial profile
    figh2a_rad = figure('color', 'white');
    z = plot(rax2a,radvg2a*100, 'LineWidth', 4);
    ylim([0 20])
    set(gca,'fontsize',fontsize2)
    title('Average radial coupling - [FNM2 #21]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    
    
    %################################################################
    %################## Compare to ideal sim: Vort Cent and WFE On
    %################################################################
    fprintf('\n-- Analyze ideal simulation (Vort cent + WFE on)\n')
    %-- Get simulated data 
    % Get simulation filename
    inds  = contains(an_params.STRNMS, 'Sim_Gauss_DON_DF1.40_FN5.22_WFEOn_VortCen');
    nmSim2a = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim2a, {'ignore'});   
    nmSim2a = nmSim2a(inds);    
    
    %-- Load simulation data
    sim2a = fitsread(nmSim2a);    
    
    %-- Find Null
        %Note: cropVal=2 helps find null in sim data with large N
    [sim2a_eta, sim2a_ind] = VFN_An_getEta_s(sim2a,2);
    fprintf('\nNull coupling on Don [ideal F#, w/ WFE+VortCen]:  %f\n',sim2a_eta);
       
    %-- Get axes in L/D    
    [sim2a_xax, sim2a_yax, sim2a_lambdaOverD] = ...
        VFN_An_simGetAxes(nmSim2a, sim2a_ind, an_params);
    
    %-- Display 2D map
    figh2a = VFN_An_fitsDisp(sim2a, sim2a_xax, sim2a_yax);
    axis([-1.3 1.3 -1.3 1.3]);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Don Coupling - [Sim ideal F#, WFE On + Vort Cent]')    
    
    %-- Calculate average radial profile
    [sim2a_radvg,sim2a_rvec] = ...
        VFN_An_radAverage(sim2a, [sim2a_ind(2), sim2a_ind(1)],2^13);
    
    %-- Get radial axis by rescaling rvec
    sim2a_rax = sim2a_rvec/sim2a_lambdaOverD;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    sim2a_radvgPK = max(sim2a_radvg);
    fprintf('Peak average [ideal F#, WFE On+Vort Cent]: = %f\n', sim2a_radvgPK)
    fprintf(' Radius of peak average: = %f\n', sim2a_rax(sim2a_radvg==sim2a_radvgPK))
    
    %-- Display radial profile
    % Simulation data
    fighsim2a_rad = figure('color', 'white');
    z = plot(sim2a_rax,sim2a_radvg*100, 'LineWidth', 4);
    xlim([0 1.2])
    set(gca,'fontsize',fontsize2)
    title('Average rad coup - [Sim ideal F#, WFE On + Vort Cent]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    % Superimpose lab data
    hold on
    plot(rax2a, radvg2a*100, 'LineWidth', 4);
    legend({'Sim','Lab'}, 'FontSize', fontsize1);
    
    
    %################################################################
    %################## Compare to ideal sim: Vort shift and WFE On
    %################################################################
    fprintf('\n-- Analyze ideal simulation (Vort shift + WFE on)\n')
    %-- Get simulated data 
    % Get simulation filename
    inds  = contains(an_params.STRNMS, 'Sim_Gauss_DON_DF1.40_FN5.22_WFEOn_VortShift');
    nmSim2b = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim2b, {'ignore'});   
    nmSim2b = nmSim2b(inds);
    
    % NOTE: THE VORTEX POSITION FOR THESE SIMULATIONS WAS ACTUALLY INCORRECT
            % THE VORTEX WAS AT ~0.0952*apRad instead of ~0.1905*apRad in X
                %           and ~0.0524*apRad instead of ~0.1048*apRad in Y
                
            % However, the old values recreate the 2D coupling map better so
            % they are probably more accurate anyway.
            % Results with the actual vortex position have 'VortPosFix' prefix.
    
    %-- Load simulation data
    sim2b = fitsread(nmSim2b);
    
    %-- Find Null
        %Note: cropVal=2 helps find null in sim data with large N
    [sim2b_eta, sim2b_ind] = VFN_An_getEta_s(sim2b,2.01);
    fprintf('\nNull coupling on Don [ideal F#, w/ WFE+VortShift]:  %f\n',sim2b_eta);
       
    %-- Get axes in L/D    
    [sim2b_xax, sim2b_yax, sim2b_lambdaOverD] = ...
        VFN_An_simGetAxes(nmSim2b, sim2b_ind, an_params);
    
    %-- Display 2D map
    figh2b = VFN_An_fitsDisp(sim2b, sim2b_xax, sim2b_yax);
    axis([-1.3 1.3 -1.3 1.3]);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Don Coupling - [Sim ideal F#, WFE On + Vort Shift]')    
    
    %-- Calculate average radial profile
    % Increase radpts sampling for sim data since very large
    [sim2b_radvg,sim2b_rvec] = ...
        VFN_An_radAverage(sim2b, [sim2b_ind(2), sim2b_ind(1)],2^13);
    
    %-- Get radial axis by rescaling rvec
    sim2b_rax = sim2b_rvec/sim2b_lambdaOverD;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    sim2b_radvgPK = max(sim2b_radvg);
    fprintf('Peak average [ideal F#, WFE On+Vort Shift]: = %f\n', sim2b_radvgPK)
    fprintf(' Radius of peak average: = %f\n', sim2b_rax(sim2b_radvg==sim2b_radvgPK))
    
    %-- Display radial profile
    % Simulation data
    fighsim2b_rad = figure('color', 'white');
    z = plot(sim2b_rax,sim2b_radvg*100, 'LineWidth', 4);
    xlim([0 1.2])
    set(gca,'fontsize',fontsize2)
    title('Average rad coup - [Sim ideal F#, WFE On + Vort Shift]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    % Superimpose lab data
    hold on
    plot(rax2a, radvg2a*100, 'LineWidth', 4);
    legend({'Sim','Lab'}, 'FontSize', fontsize1);
    
    
    %################################################################
    %################## Compare to varying Df with WFE + Vort Shift
    %################################################################
    fprintf('\n-- Analyze varying Df simulation: Vort Shift + WFE On\n')
    
    %-- Get simulated data
    % Get simulation filenames
    inds  = contains(an_params.STRNMS, 'Sim_Gauss_DON');
    nmSim2c = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim2c, {'VortCen', 'VortPosFix'});   %'DF1.96'
    nmSim2c = nmSim2c(inds);
    
    % NOTE: THE VORTEX POSITION FOR THESE SIMULATIONS WAS ACTUALLY INCORRECT
            % THE VORTEX WAS AT ~0.0952*apRad instead of ~0.1905*apRad in X
                %           and ~0.0524*apRad instead of ~0.1048*apRad in Y
                
            % However, the old values recreate the 2D coupling map better so
            % they are probably more accurate anyway.
            % Results with the actual vortex position have 'VortPosFix' prefix.
    
    %-- Iterate through reading and processing simulated data
    simrpts = 4000;     % 4000 since ~2^13/2 radial samples
    % Preallocate arrays for results
    sim2c_radvg = nan(length(nmSim2c), simrpts); 
    sim2c_rax = nan(length(nmSim2c),simrpts);
    % Create figure for plotting
    fighsim2c = figure('color', 'white');
    hold on;
    % Clear legLabels to avoid problems with old legend stuff
    legLabels = {};
    for i = 1:length(nmSim2c)
        
        %-- Load simulation data
        sim2c = fitsread(nmSim2c(i));
        
        %-- Find null
        %Note: cropVal=2 helps find null in sim data with large N
        [sim2c_eta, sim2c_ind] = VFN_An_getEta_s(sim2c,2.01);
        
        %-- Get just name part of filename
        [~, tmpnm, ~] = fileparts(nmSim2c(i));
        tmpnm = erase(tmpnm,[{'Sim_'},{'_FN5.22'},{'Gauss_DON_'},{'WFEOn_VortShift'}]);
        
        %-- Print null coupling value
        fprintf('\n-%s Values\n', tmpnm);
        fprintf(' Null coupling on Donut [Sim w/ WFE+VortShift]:  %f\n',sim2c_eta);
        
        %-- Get lambda/D value
        [sim2c_xax, sim2c_yax, lOverD] = ...
            VFN_An_simGetAxes(nmSim2c(i), sim2c_ind, an_params);
        
        %-- Calculate average radial profile
        [sim2c_radvg(i,:),tmpRvec] = ...
            VFN_An_radAverage(sim2c, [sim2c_ind(2), sim2c_ind(1)], simrpts);
        %sim1b_radvg(i,:) = tmpRadvg;
        
        %-- Get radial axis by rescaling rvec
        sim2c_rax(i,:) = tmpRvec/lOverD;
               
        %-- Print Peak
        %Average for the donut ring with the peak coupling
        sim2c_radvgPK = max(sim2c_radvg(i,:));
        fprintf(' Peak average [Sim w/ WFE+VortShift]: = %f\n', sim2c_radvgPK)
        tmprax = sim2c_rax(i,:);
        fprintf('  Radius of peak average: = %f\n', tmprax(sim2c_radvg(i,:)==sim2c_radvgPK))
        
%         %-- Display 2D map for each DF
%         VFN_An_fitsDisp(sim2c, sim2c_xax, sim2c_yax);
%         axis([-1.3 1.3 -1.3 1.3]);
%         % Use parula for this data since want to see nuances of asymmetry
%         colormap parula
%         title(strcat('Don Coupling - ',tmpnm), 'Interpreter', 'none');
%         
%         %-- Reselect figure for line profile
%         figure(fighsim2c);
        
        %-- Plot results
        if contains(tmpnm, {'DF1.09', 'DF1.4', 'DF1.6'})
            linstyle = ':';
        else
            linstyle = '--';
        end
        plot(sim2c_rax(i,:),sim2c_radvg(i,:)*100, linstyle, 'LineWidth', 1.5);
        
        %-- Add legend label
        legLabels{i} = tmpnm;
    end
    
    %-- Finish formatting figure
    % Superimpose lab data
    plot(rax2a, radvg2a*100, 'LineWidth', 3);
    hold off
    set(gca,'fontsize',fontsize2)
    lablabel = sprintf('Lab: PZ = %0.2f', 1/an_params.V2um);
    legLabels = [legLabels {lablabel}];
    % Add legend
    legend(legLabels, 'Interpreter', 'none', 'FontSize', fontsize1, 'Location','southeast');
    % Add title and axes
    title(sprintf('Coupling - Simulated (Shift) vs Lab PZGain = %0.2f', 1/an_params.V2um));
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1);
    xlim([0 1.2])
    
    if false
        svnm = sprintf('DONCoup_SimVortShiftVsLab_PZ%0.2f_w6.5FibMod.png', 1/an_params.V2um);
        export_fig([savefld filesep svnm],'-r300', '-painters');
    end
    
    fprintf('\n---- End Section %i \n\n',seci)
end

%% Section 4: OLD FNUM Donut scan for Null - FTO3 [31]
seci    = 4;
if Sec2Run(seci)
    fprintf('\n---- Section %i ----\n',seci)
    
    % SET FNUM TO OLD VALUE SINCE ANALYZING DATA WITH DIFF PUP DIAMETER
    um2LD_OLD = an_params.um2LD;
    an_params.um2LD = 1/(0.635*10.96/3.59);  %(1/lambda*F#)
    
    %################################################################
    %################## Basic analysis (Min in frame) 
    %################################################################
    fprintf('\n-- Basic analysis\n')
    
    %-- Get normalized data
    nm4a    = ['101718_FTO3' filesep '31Don_NullConf'];
    nrm4a   = VFN_An_fitsNormed(nm4a, [1,1,1], an_params);
    
    %-- Calculate null coupling
    [eta_4a, ind4a] = VFN_An_getEta_s(nrm4a);
    fprintf('\nNull coupling on Donut [FTO3 #31]:    %f\n',eta_4a);
    if isPrintLoc
        fprintf([' Null Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 ind4a(2), ind4a(1), ind4a(4), ind4a(5), ind4a(6));
    end
    
    %-- Get axes for image (relative to null point):
    [xax4a, yax4a]  = VFN_An_getAxes(nm4a, ind4a, an_params);
    
    %-- Show Frame
    figh4a = VFN_An_fitsDisp(nrm4a, xax4a, yax4a);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Donut Coupling - [FTO3 #31]')
        
    %################################################################
    %################## Radial Average  
    %################################################################
    fprintf('\n-- Radial Average\n')
    %-- Calculate average radial profile
    [radvg4a,rvec4a] = VFN_An_radAverage(nrm4a, [ind4a(2), ind4a(1)]);
    
    %-- Get radial axis 
    % Use xax to get L/D per 'pixel' value then rescale rvec by this value
    rax4a = (xax4a(2)-xax4a(1))*rvec4a;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    radvgPK4a = max(radvg4a);
    fprintf('Peak average [FTO3 #31]: = %f\n', radvgPK4a)
    fprintf(' Radius of peak average: = %f\n', rax4a(radvg4a==radvgPK4a))
    
    %-- Display radial profile
    figh4a_rad = figure('color', 'white');
    z = plot(rax4a,radvg4a*100, 'LineWidth', 4);
    ylim([0 20])
    set(gca,'fontsize',fontsize2)
    title('Average radial coupling - [FTO3 #31]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    
    %################################################################
    %################## Compare to ideal sim: Vort Cent and WFE On
    %################################################################
    fprintf('\n-- Analyze ideal simulation (Vort cent + WFE on)\n')
    %-- Get simulated data 
    % Get simulation filename
    inds  = contains(an_params.STRNMS, 'SimOLD_Gauss_DON_DF1.40_FN3.04_WFEOn_VortShift');
    nmSim4a = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim4a, {'ignore'});   
    nmSim4a = nmSim4a(inds);    
    
    %-- Load simulation data
    sim4a = fitsread(nmSim4a);    
    
    %-- Find Null
        %Note: cropVal=2 helps find null in sim data with large N
    [sim4a_eta, sim4a_ind] = VFN_An_getEta_s(sim4a,2.01);
    fprintf('\nNull coupling on Don [ideal F#, w/ WFE+VortCen]:  %f\n',sim4a_eta);
       
    %-- Get axes in L/D    
    [sim4a_xax, sim4a_yax, sim4a_lambdaOverD] = ...
        VFN_An_simGetAxes(nmSim4a, sim4a_ind, an_params);
    
    %-- Display 2D map
    figh4a = VFN_An_fitsDisp(sim4a, sim4a_xax, sim4a_yax);
    axis([-1.3 1.3 -1.3 1.3]);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Don Coupling - [Sim ideal F#, WFE On + Vort Cent]')    
    
    %-- Calculate average radial profile
    [sim4a_radvg,sim4a_rvec] = ...
        VFN_An_radAverage(sim4a, [sim4a_ind(2), sim4a_ind(1)],2^13);
    
    %-- Get radial axis by rescaling rvec
    sim4a_rax = sim4a_rvec/sim4a_lambdaOverD;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    sim4a_radvgPK = max(sim4a_radvg);
    fprintf('Peak average [ideal F#, WFE On+Vort Cent]: = %f\n', sim4a_radvgPK)
    fprintf(' Radius of peak average: = %f\n', sim4a_rax(sim4a_radvg==sim4a_radvgPK))
    
    %-- Display radial profile
    % Simulation data
    fighsim4a_rad = figure('color', 'white');
    z = plot(sim4a_rax,sim4a_radvg*100, 'LineWidth', 4);
    xlim([0 1.2])
    set(gca,'fontsize',fontsize2)
    title('Average rad coup - [Sim ideal F#, WFE On + Vort Cent]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    % Superimpose lab data
    hold on
    plot(rax4a, radvg4a*100, 'LineWidth', 4);
    legend({'Sim','Lab'}, 'FontSize', fontsize1);
    
    %################################################################
    %################## Compare to ideal sim: Vort shift and WFE On
    %################################################################
    fprintf('\n-- Analyze ideal simulation (Vort shift + WFE on)\n')
    %-- Get simulated data 
    % Get simulation filename
    inds  = contains(an_params.STRNMS, 'SimOLD_Gauss_DON_DF1.40_FN3.04_WFEOn_VortShift');
    nmSim4b = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim4b, {'ignore'});   
    nmSim4b = nmSim4b(inds);
    
    % NOTE: THE VORTEX WAS MOVED TO 1/2 WHAT THE ZABER POSITIONS WOULD SAY IT
        % SHOULD BE AT. THIS IS TO MATCH THE 1/2 VALUE USED FOR SECTION 2 WHICH
        % PRODUCED SUCH GREAT RESULTS.
        % THIS 1/2 AGAIN PRODUCES AN ASTOUNDINGLY ACCURATE 2D MAP
    
    %-- Load simulation data
    sim4b = fitsread(nmSim4b);
    
    %-- Find Null
        %Note: cropVal=2 helps find null in sim data with large N
    [sim4b_eta, sim4b_ind] = VFN_An_getEta_s(sim4b,2.01);
    fprintf('\nNull coupling on Don [ideal F#, w/ WFE+VortShift]:  %f\n',sim4b_eta);
       
    %-- Get axes in L/D    
    [sim4b_xax, sim4b_yax, sim4b_lambdaOverD] = ...
        VFN_An_simGetAxes(nmSim4b, sim4b_ind, an_params);
    
    %-- Display 2D map
    figh4b = VFN_An_fitsDisp(sim4b, sim4b_xax, sim4b_yax);
    axis([-1.3 1.3 -1.3 1.3]);
    % Use parula for this data since want to see nuances of asymmetry
    colormap parula
    title('Don Coupling - [Sim ideal F#, WFE On + Vort Shift]')    
    
    %-- Calculate average radial profile
    % Increase radpts sampling for sim data since very large
    [sim4b_radvg,sim4b_rvec] = ...
        VFN_An_radAverage(sim4b, [sim4b_ind(2), sim4b_ind(1)],2^13);
    
    %-- Get radial axis by rescaling rvec
    sim4b_rax = sim4b_rvec/sim4b_lambdaOverD;
    
    %-- Print Peak
    %Average for the donut ring with the peak coupling
    sim4b_radvgPK = max(sim4b_radvg);
    fprintf('Peak average [ideal F#, WFE On+Vort Shift]: = %f\n', sim4b_radvgPK)
    fprintf(' Radius of peak average: = %f\n', sim4b_rax(sim4b_radvg==sim4b_radvgPK))
    
    %-- Display radial profile
    % Simulation data
    fighsim4b_rad = figure('color', 'white');
    z = plot(sim4b_rax,sim4b_radvg*100, 'LineWidth', 4);
    xlim([0 1.2])
    set(gca,'fontsize',fontsize2)
    title('Average rad coup - [Sim ideal F#, WFE On + Vort Shift]')
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Coupling [%]', 'FontSize', fontsize1);
    % Superimpose lab data
    hold on
    plot(rax4a, radvg4a*100, 'LineWidth', 4);
    legend({'Sim','Lab'}, 'FontSize', fontsize1);
        
    %################################################################
    %################## Compare to varying Df with WFE + Vort Shift
    %################################################################
    fprintf('\n-- Analyze varying Df simulation: Vort Shift + WFE On\n')
    
    %-- Get simulated data
    % Get simulation filenames
    inds  = contains(an_params.STRNMS, 'SimOLD_Gauss_DON'); %  'SimOLD_BigConv_Gauss_DON'
        %NOTE: The BigConv data is more reliable for big fiber modes since a
        %larger cropping region was used for the convolution in the simulation.
        %See simulation code for more details.
    nmSim4c = an_params.STRNMS(inds);
    % Remove specified simulations  (note: something should always be in 
        % the char array else all names get removed. Use 'ignore' if needed
    inds  = ~contains(nmSim4c, {'VortCen'});  
    nmSim4c = nmSim4c(inds);
    
    % NOTE: THE VORTEX WAS MOVED TO 1/2 WHAT THE ZABER POSITIONS WOULD SAY IT
        % SHOULD BE AT. THIS IS TO MATCH THE 1/2 VALUE USED FOR SECTION 2 WHICH
        % PRODUCED SUCH GREAT RESULTS.
        % THIS 1/2 AGAIN PRODUCES AN ASTOUNDINGLY ACCURATE 2D MAP    
    
    %-- Iterate through reading and processing simulated data
    simrpts = 4000;     % 4000 since ~2^13/2 radial samples
    % Preallocate arrays for results
    sim4c_radvg = nan(length(nmSim4c), simrpts); 
    sim4c_rax = nan(length(nmSim4c),simrpts);
    % Create figure for plotting
    fighsim4c = figure('color', 'white');
    hold on;
    % Clear legLabels to avoid problems with old legend stuff
    legLabels = {};
    for i = 1:length(nmSim4c)
        
        %-- Load simulation data
        sim4c = fitsread(nmSim4c(i));
        
        %-- Find null
        %Note: cropVal=2 helps find null in sim data with large N
        [sim4c_eta, sim4c_ind] = VFN_An_getEta_s(sim4c,2.03);
        
        %-- Get just name part of filename
        [~, tmpnm, ~] = fileparts(nmSim4c(i));
        tmpnm = erase(tmpnm,[{'SimOLD_'},{'_FN3.04'},{'Gauss_DON_'},{'WFEOn_VortShift'}, {'BigConv_'}]);
        
        %-- Print null coupling value
        fprintf('\n-%s Values\n', tmpnm);
        fprintf(' Null coupling on Donut [Sim w/ WFE+VortShift]:  %f\n',sim4c_eta);
        
        %-- Get lambda/D value
        [sim4c_xax, sim4c_yax, lOverD] = ...
            VFN_An_simGetAxes(nmSim4c(i), sim4c_ind, an_params);
        
        %-- Calculate average radial profile
        [sim4c_radvg(i,:),tmpRvec] = ...
            VFN_An_radAverage(sim4c, [sim4c_ind(2), sim4c_ind(1)], simrpts);
        %sim1b_radvg(i,:) = tmpRadvg;
        
        %-- Get radial axis by rescaling rvec
        sim4c_rax(i,:) = tmpRvec/lOverD;
               
        %-- Print Peak
        %Average for the donut ring with the peak coupling
        sim4c_radvgPK = max(sim4c_radvg(i,:));
        fprintf(' Peak average [Sim w/ WFE+VortShift]: = %f\n', sim4c_radvgPK)
        tmprax = sim4c_rax(i,:);
        fprintf('  Radius of peak average: = %f\n', tmprax(sim4c_radvg(i,:)==sim4c_radvgPK))
        
        %-- Display 2D map for each DF
        VFN_An_fitsDisp(sim4c, sim4c_xax, sim4c_yax);
        axis([-1.3 1.3 -1.3 1.3]);
        % Use parula for this data since want to see nuances of asymmetry
        colormap parula
        title(strcat('Don Coupling - ',tmpnm), 'Interpreter', 'none');
        
        %-- Reselect figure for line profile
        figure(fighsim4c);
        
        %-- Plot results
        if contains(tmpnm, {'DF1.87', 'DF1.4', 'DF2.75'})
            linstyle = ':';
        else
            linstyle = '--';
        end
        plot(sim4c_rax(i,:),sim4c_radvg(i,:)*100, linstyle, 'LineWidth', 1.5);
        
        %-- Add legend label
        legLabels{i} = tmpnm;
    end
    
    %-- Finish formatting figure
    % Superimpose lab data
    plot(rax4a, radvg4a*100, 'LineWidth', 3);
    hold off
    set(gca,'fontsize',fontsize2)
    lablabel = sprintf('Lab: PZ = %0.2f', 1/an_params.V2um);
    legLabels = [legLabels {lablabel}];
    % Add legend
    legend(legLabels, 'Interpreter', 'none', 'FontSize', fontsize1, 'Location','northwest');
    % Add title and axes
    title(sprintf('Coupling - Simulated (Shift) vs Lab PZGain = %0.2f', 1/an_params.V2um));
    xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
    ylabel('Planet Throughput [%]', 'FontSize', fontsize1);
    xlim([0 1.7])
    
    if false
        svnm = sprintf('OLDDONCoup_SimVortShiftVsLab_PZ%0.2f_w6.5FibMod.png', 1/an_params.V2um);
        export_fig([savefld filesep svnm],'-r300', '-painters');
    end
    
    
    % SET FNUM BACK TO VALUE FOR NEW SYSTEM 
     an_params.um2LD = um2LD_OLD;
    
    fprintf('\n---- End Section %i \n\n',seci)
end