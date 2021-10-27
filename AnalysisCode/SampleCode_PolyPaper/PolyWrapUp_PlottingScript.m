%{
Plotting Script for analyzing the data for the Polychromatic VFN paper

This script is based off the previous MonoWrapUp_PlotingScript_Updated.m
file. It does a basic analysis of all the the donut scans performed for the
Polychromatic VFN paper. As such, it loads the raw data and  processes it using
the user-provided parameters to generate the normalize 2D coupling maps
(Donuts). Those maps are then used to compute the null depth and average
planet coupling as a function of separation. The code generates several
figures displaying the results and also reports some of the numbers in the
command window.

These are the figures used for the VFN Paper.

The dataset that was analyzed for this was taken on 07/20/2021 and is saved
in the COV8 folder. See the related DataNotes.txt file for details. The
scans are identical to each other in parameters, the only thing that
changes between each scan is that the light source (superK varia) was set
to increasing bandwidths for each scan. 

As such, the analysis done for all the scans is identical and a single
analysis function is defined at the end of the script. 

NOTE: this script was originally run on the VFN Server (Ubuntu).

Author:         Daniel Echeverri
Origin Date:    07/22/2021
%}
close all; clear all;
%% Setup Code Paths
% add path to intelligent polar transform function 
addpath('../../../VFN-Simulations/VFNlib')

% add path to analysis library
addpath('../AnalysisLib')
%% User Inputs (Part 1) - Figure/Save properties
global svpms

%-- Save Information
% Folder into which data should be saved
savefld = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/PublishedInPaper';

% Fyletype for most saved figures
svpms.fltype  = '.png';  

% Default fonts for saved figures
% Fontsizes are in 'normalized' units. That means the value should be
  % fractional [0 to 1]. Ex: 0.1 sets text to 10% of the figure height.
svpms.fontsize1 = 0.07;     % Title fontsize
svpms.fontsize2 = 0.05;     % Other text fontsize

% Line properties
svpms.linewidth = 2.3;


%-- Gif Properties
gifnm     = 'VFN_650xBWScan.gif';   % Filename for the gif
gifDelay  = 0.7;        % [s] Delay between frames in gif

%-- Line profile (single scan) properties
svpms.gridmodeLP  = 'on';       % Show the grid? Should be 'on', 'off', 'minor', 'major', etc.
svpms.figdimLP = [7.6,6.5];  % Figure dimensions [width, height] in inches

%-- All line profiles overlayed figure properties
svpms.gridmodeAllLP  = 'on';       % Show the grid? Should be 'on', 'off', 'minor', 'major', etc.
svpms.figdimAllLP = [3 3 9.5 6.5];  % Fig dim and size [xpos, ypos, width, height] in inches

%-- Null and Planet vs. properties
svpms.gridmodeNPvBW  = 'minor';    % Show the grid? Should be 'on', 'off', 'minor', 'major', etc.
svpms.figdimNPvBW    = [3, 3, 7, 6.5]; % Fig dim and size [xpos, ypos, width, height] in inches
svpms.y_limsNPvBW    = [0,10];  % This is the y-limit for both left and right axes (ie. covers the null and planet coup)

%-- 2D Coupling Map (donut) properties
svpms.figdimDON  = [7.6,6.5];  % Figure dimensions [width, height] in inches


%-- Flags to choose what to save
isSaveFig = true;      % Save figures?
isSaveMat = true;      % Save .mat with key results?
isSaveGif = true;       % Save gif scannig through donuts?

%% User Inputs (Part 2) - Analysis Parameters
%-- Folder with data (Should be the folder containing all the VFN data folders)
fldnm = '/media/Data_Drive/VFN/TestbedData';

%-- Data to consider
% Cell array containing sub-directories to consider
datadirs = {'*COV8*'};
% Identifier numbers for the datasets to consider
  % First number in the filename (ex: "292" in 292Don_650PapPolyRep3)
runnumbers = 292:302;
% Bandwidths for the chosen datasets (ordered to match runnumbers)
BWs = [6,10:10:100];    % [nm] 

%-- Analysis properties
% Scaling factor for converting between femto gain settings
an_params.gnFactor = 10.01;     % Value measured from Femto+superK calibration in 07/2021

% Fiber Losses
an_params.FR = 0.9654;      % Fresnel reflection term in Nem's equation
an_params.PL = 0.9966;      % Propagation loss term in nem's equation.

% System settings
an_params.lambda = 650e-9;  % [m] Central wavelength
an_params.foc    = 11e-3;   % [m] EFL @ design wavelength (670nm) from Thorlabs spec sheet for A397 lens
an_params.D      = 2.45e-3; % [m] VFN pupil diameter from measurements near [COV7 #18]

% Power meter scaling factor (redPM to Femto ratio)
  % Note: this value is slightly BW-dependent (shifts to 0.64 at BW>70nm)
  % but the difference is so small that for simplicity I'll use 0.63 for
  % all scans
an_params.rdOfmto = 1/0.63;   % [V/uW] Value measured from Femto+superK calibration in 07/2021

% Conversion factor from microns to lambda/D in focal plane:
  % Note: 1lam/D in anglular scale = 1lam*F# in spatial scale
an_params.Fnum   = an_params.foc/an_params.D;    % []  System F/#
an_params.LD2um  = (an_params.lambda*an_params.Fnum)*1e6;  % [um / lambda/D] microns PER lambda/D in image plane
an_params.um2LD  = 1/an_params.LD2um;           % [lambda/D / um] lambda/D per micron in image plane

%-- Display the System properties
disp(an_params)

%-- Find all .fits files in the desired directories
nmWpath = fullfile(fldnm,datadirs,'*.fits');
an_params.STRNMS = VFN_An_getSTRNMS(nmWpath);

%% Iterate through and analyze all bandwidths
ind = 1;
for frmnum = runnumbers
    fprintf("\n______%dnm BW [COV8 #%d]:______\n",BWs(ind),frmnum)

    %-- Identifiers
    runstr = sprintf('%dnm BW Scan',BWs(ind)); % string identifier for print statements
    runid  = sprintf('[COV8 #%d]',frmnum);  % Numeric identifier for figure title
    svtag  = sprintf('COV8_%d_%dnmBW',frmnum,BWs(ind));  % Tag to use for filename in saved figure

    %-- File to find
    frmnm   = ['COV8' filesep sprintf('%dDon',frmnum)];

    %-- Perform analysis
    results(ind) = basicAnalyzer(frmnm, an_params, runstr, runid, true);
    
    %-- Save figures if requested
    if isSaveFig
        export_fig(results(ind).figh_donut,fullfile(savefld,[svtag '_CoupMap' svpms.fltype]),'-r300', '-painters');
        export_fig(results(ind).figh_radprof,fullfile(savefld,[svtag '_CoupProf' svpms.fltype]),'-r300', '-painters');
    end
    
    %-- Increment index 
    ind = ind + 1;
end

%% Get ideal coupling curve

% Read Charge 1 Profile file
%CH1 = readmatrix('/media/Data_Drive/VFN/PSISIM_Work/ETC/VFN/Charge1_Ideal.txt');
% Read Charge 2 Profile file
CH2 = readmatrix('/media/Data_Drive/VFN/PSISIM_Work/ETC/VFN/Charge2_Ideal.txt');

%-- Get coupling curve from model w/ fibermismatch
%model = readmatrix('/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/VFN_PolyPaper_ModeledCoupProf_4.49Fnum_ModeAutoGauss4.3umCoreRad.csv');

%-- Get coupling curve for models with varying MFD
MFDs = 3.6:0.2:5.3;
ind = 1;
for mfd = MFDs
    % Folder with model files
    modnm = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/PresentedAtSPIE/';
    % Format filename (changing just the MFD for this dataset)
    modnm = [modnm 'VFN_PolyPaper_ModeledCoupProf_4.49Fnum_MFDManual_gaussian_' num2str(mfd,'%0.2f') 'umMFD.csv'];
    % Read the file (as a struct so we can deal with it easier)
    model_tmp.LP = readmatrix(modnm);
    % Add to struct array
    models_MFDscan(ind) = model_tmp;
    % Incremenet counter
    ind = ind +1;
end
%% Make plot overlaying coupling curves with ideal

% Preallocate legend labels (+1 so we can add ideal label afterwards)
labels = cell(length(results)+1,1);
figCoupCurvs = figure('color','white','units', 'inches', 'Position', svpms.figdimAllLP);
% Iterate through results, adding the radial profile to the plot
for ind = 1:length(results)
    % Plot the coupling curve
    plot(results(ind).rax,results(ind).radavg*100,'LineWidth',svpms.linewidth)
    % Make formatted label for legend
    %labels{ind} = sprintf(' %0.1f',(BWs(ind)*1e-9)/an_params.lambda*100);
    labels{ind} = sprintf(' %d', BWs(ind));
    if ind == 1
        % set figure to overlay only after first plot
        hold on;
    end
end

% Add the ideal coupling curve
plot(CH2(:,1), CH2(:,2)*100, '--', 'LineWidth', svpms.linewidth);%, 'Color', get(ztb, 'Color'));
labels{ind+1} = ' (Ideal)';
% Add the model
%plot(model(:,1),model(:,2)*100, '--', 'LineWidth', svpms.linewidth);
%labels{ind+2} = ' Model F/4.5';
% Format the figure
grid(gca, svpms.gridmodeAllLP);
xlim([0 2.5])
ylim([0 11])
xlabel('Angular Separation [\lambda/D]')
ylabel({'Coupling [%]'})
set(gca, 'FontUnits', 'normalized', 'fontsize', svpms.fontsize2)
leg = legend(labels, 'Location', 'eastoutside');
%title(leg, {'Bandwidth', '\Delta\lambda/\lambda_0 [%]'})
title(leg, {'Bandwidth [nm]'})
tit = {'VFN Charge 2 Coupling Curves', ['(\lambda_0 = ', num2str(an_params.lambda*1e9,'%d'), 'nm)']};
title(tit, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1)


if isSaveFig
    export_fig(figCoupCurvs,fullfile(savefld,['VFN_CoupCurve_vs_Bandwidth', svpms.fltype]),'-r300', '-painters');
end

%% Same overlayed coupling curve plot but in loglog scale

% Preallocate legend labels (+1 so we can add ideal label afterwards)
labels = cell(length(results)+1,1);
figCoupCurvesLog = figure('color','white','units', 'inches', 'Position', svpms.figdimAllLP);
% Iterate through results, adding the radial profile to the plot
for ind = 1:length(results)
    % Plot the coupling curve
    loglog(results(ind).rax,results(ind).radavg*100,'LineWidth',svpms.linewidth)
    % Make formatted label for legend
    %labels{ind} = sprintf(' %0.1f',(BWs(ind)*1e-9)/an_params.lambda*100);
    labels{ind} = sprintf(' %d', BWs(ind));
    if ind == 1
        % set figure to overlay only after first plot
        hold on;
    end
end

% Add the ideal coupling curve
loglog(CH2(:,1), CH2(:,2)*100, '--', 'LineWidth', svpms.linewidth);%, 'Color', get(ztb, 'Color'));
labels{ind+1} = ' (Ideal)';
% Add the model
%plot(model(:,1),model(:,2)*100, '--', 'LineWidth', svpms.linewidth);
%labels{ind+2} = ' Model F/4.5';
% Format the figure
grid(gca, svpms.gridmodeAllLP);
xlim([0 2.5])
ylim([0 11])
xlabel('Angular Separation [\lambda/D]')
ylabel({'Coupling [%]'})
set(gca, 'FontUnits', 'normalized', 'fontsize', svpms.fontsize2)
leg = legend(labels, 'Location', 'eastoutside');
%title(leg, {'Bandwidth', '\Delta\lambda/\lambda_0 [%]'})
title(leg, {'Bandwidth [nm]'})
tit = {'VFN Charge 2 Coupling Curves', ['(\lambda_0 = ', num2str(an_params.lambda*1e9,'%d'), 'nm)']};
title(tit, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1)

%% Make coupline profile plot with the 6nm scan and multiple MFD models

% Preallocate legend labels (+1 so we can add ideal label afterwards)
labels = cell(length(MFDs)+1,1);
figCoupCurvesModels = figure('color','white','units', 'inches', 'Position', svpms.figdimAllLP);
% Iterate through results, adding the radial profile to the plot
for ind = 1:length(MFDs)
    % Plot the coupling curve
    plot(models_MFDscan(ind).LP(:,1),models_MFDscan(ind).LP(:,2)*100,'LineWidth',svpms.linewidth)
    % Make formatted label for legend
    %labels{ind} = sprintf(' %0.1f',(BWs(ind)*1e-9)/an_params.lambda*100);
    labels{ind} = sprintf(' %0.2f', MFDs(ind));
    if ind == 1
        % set figure to overlay only after first plot
        hold on;
    end
end

% Add the measured experimental 6nm curve
plot(results(1).rax,results(1).radavg*100,'--','LineWidth',svpms.linewidth)
labels{ind+1} = ' (6nm BW Actual)';
% Format the figure
grid(gca, svpms.gridmodeAllLP);
xlim([0 2.5])
ylim([0 11])
xlabel('Angular Separation [\lambda/D]')
ylabel({'Coupling [%]'})
set(gca, 'FontUnits', 'normalized', 'fontsize', svpms.fontsize2)
leg = legend(labels, 'Location', 'eastoutside');
%title(leg, {'Bandwidth', '\Delta\lambda/\lambda_0 [%]'})
title(leg, {'MFD [um]'})
tit = {'VFN Charge 2 Model Matching', ['(\lambda_0 = ', num2str(an_params.lambda*1e9,'%d'), 'nm)']};
title(tit, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1)

%% Make plot showing null & planet coup vs Bandwidth

% Get all the nulls
nulls = [results.eta_s];
% Get all the (average) peak planet couplings
peaks = [results.radmax];

%-- Make the figure
figNullvBW = figure('color','white','units', 'inches', 'Position', svpms.figdimNPvBW);
% Plot the nulls
yyaxis left
plot((BWs*1e-9)/an_params.lambda, nulls*1e5, '-d', 'MarkerSize', 10, 'LineWidth', svpms.linewidth);
ylabel({'Null Depth [x10^{-5}]'})
ylim(svpms.y_limsNPvBW)
% Plot the peaks
yyaxis right
plot((BWs*1e-9)/an_params.lambda, peaks*100,'-s', 'MarkerSize', 10, 'LineWidth', svpms.linewidth);
ylim(svpms.y_limsNPvBW)
ylabel({'Planet Coupling [%]'})
xlabel('Bandwidth [\Delta\lambda/\lambda_0]')
grid(gca, svpms.gridmodeNPvBW)
set(gca, 'FontUnits', 'normalized', 'fontsize', svpms.fontsize2)
tit = {'VFN Broadband Performance', ['(\lambda_0 = ', num2str(an_params.lambda*1e9,'%d'), 'nm)']};
title(tit, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1)


if isSaveFig
    export_fig(figNullvBW,fullfile(savefld,['VFN_NullAndPlanetCoup_vs_Bandwidths',svpms.fltype]),'-r300', '-painters');
end

%% Save .mat file
% Do this before making the gif since that part of the code modifies the X
% and Y axes of the figures slightly.
if isSaveMat
    savenm = fullfile(savefld, ['AnalysisData_' datestr(now,'mmddyy_HH_MM') '.mat']);
    save(savenm, 'results', 'an_params', 'datadirs', 'runnumbers', 'BWs', 'figNullvBW', 'figCoupCurvs', 'savefld', 'svpms');
end

%% Make a gif showing the donut at each bandwidth
% Note: we already made the figures like we want them so this section just
% loads them one at a time so they can be saved into a gif

if isSaveGif
    for ind = 1:length(results)
        % Make sure that all tick markers are the same
        if ind == 1
            % For the first frame, set the tick markers to be integers and
            % then save the tick markers to an array we can access in
            % remaining frames.
            ax = results(ind).figh_donut.Children;
            
            ax.XTick = unique( round(ax.XTick) );
            ax.YTick = unique( round(ax.YTick) );
            
            gifXtick = ax.XTick; gifYtick = ax.YTick;
        else
            % For remaining frames, set the tick markers to match the first
            % frame.
            ax = results(ind).figh_donut.Children;
            ax.XTick = gifXtick; ax.YTick = gifYtick;
        end
        
        frame = getframe(results(ind).figh_donut);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if ind == 1
            imwrite(imind,cm,fullfile(savefld,gifnm),'gif', 'Loopcount',inf, 'DelayTime', gifDelay);
        else
            imwrite(imind,cm,fullfile(savefld,gifnm),'gif','WriteMode','append','DelayTime', gifDelay);
        end
    end
end
%%
% for ind = 1:length(results)
%     if ind == 1
%         % First iteration so let's load the full first figure
%         ax = imagesc(surfCube(:,:,1,aind));
%         axis image; colorbar
%         caxis(clim)
%         set(gca,'xtick',[], 'ytick', [])
%         titH = title({tit, sprintf(subtit,amplitudes(aind))});
%     else
%         set(ax, 'CData', surfCube(:,:,1,aind))
%         titH.String = {tit, sprintf(subtit,amplitudes(aind))};
%     end
%     %ax.Parent.CLim = [min(surfCube(:,:,ind),[],'all'), max(surfCube(:,:,ind),[],'all')];
%     drawnow;
% 
%     if isSaveGif
%         % Add frame
%         frame = getframe(figgif);
%         im = frame2im(frame);
%         [imind, cm] = rgb2ind(im, 256);
%         if frstflg
%             imwrite(imind,cm,fullfile(svfld,gifnm),'gif', 'Loopcount',inf, 'DelayTime', gifDelay);
%             frstflg = false;
%         else
%             imwrite(imind,cm,fullfile(svfld,gifnm),'gif','WriteMode','append','DelayTime', gifDelay);
%         end
%     end
% end

%% Function to perform a basic analysis 
function res = basicAnalyzer(frmnm, an_params, runstr, runid, verbose)
    global svpms
    %-- Get the normalized data (full eta map)
    res.frmNrmd = VFN_An_fitsNormed(frmnm, [1,1,1], an_params);

    %-- Extract eta_s and location
    [res.eta_s, res.eta_s_ind] = VFN_An_getEta_s(res.frmNrmd);

    %-- Extract max coupling and location
    % NOTE: this is just the global max of the frame. This isn't a reliable metric 
    % for the planet coupling since only a single point hit this value hence the planet 
    % would need to land at this particular point.
    [res.eta_max, res.eta_max_ind] = VFN_An_getEta_p(res.frmNrmd);

    %-- Calculate Average Raidal Coupling
    % This is a much better metric since it provides the average coupling as a function 
    % of radial distance from the star. Hence, given a known planet separation, this 
    % provides the average coupling we can expect for that planet without knowing 
    % it's azimuthal location.
    [res.radavg, res.rvec] = VFN_An_radAverage(res.frmNrmd, [res.eta_s_ind(2), res.eta_s_ind(1)]);
    [res.radmax, res.rad_max_ind] = max(res.radavg);

    %-- Get axes for for plotting
    % Note: This sets the origin of the axes at the null point.
    % First get the x,y axes for donut images
    [res.xax, res.yax] = VFN_An_getAxes(frmnm, res.eta_s_ind, an_params); 
    % Now get the x-axis for radial profile (coupling curve). "rvec", the result 
    % from VFN_An_radAverage, is in units of "pixels (ie. number fiber steps from 
    % the null point). Thus, we can get the pixel size, in units of L/D, from the 
    % properly-scaled xax returned by VFN_An_getAxes
      % (use xax to get L/D per 'pixel' then use it to rescale rvec)
    res.rax = (res.xax(2)-res.xax(1))*res.rvec;

    %-- Print results
    if verbose
        % Null and location
        fprintf('\nNull in %s:    %e\n',runstr,res.eta_s)
        fprintf([' eta_s Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 res.eta_s_ind(2), res.eta_s_ind(1), res.eta_s_ind(4), ...
                 res.eta_s_ind(5), res.eta_s_ind(6));

        % Peak and location
        fprintf('\nPeak in %s:    %f\n',runstr,res.eta_max);
        fprintf([' eta_p Location: ', ...
                 '\n  X index = %3i', ...
                 '\n  Y index = %3i', ...
                 '\n  Foc ind = %3i', ...
                 '\n  VX ind  = %3i', ...
                 '\n  VY ind  = %3i\n'], ...
                 res.eta_max_ind(2), res.eta_max_ind(1), res.eta_max_ind(4), ...
                 res.eta_max_ind(5), res.eta_max_ind(6));
        fprintf('\nPeak RadAvg in %s:    %f\n', runstr, res.radmax);
        fprintf(' Peak RadAvg at %f lam/D\n', res.rax(res.rad_max_ind));
        %-- Print relative integration time
        fprintf('\nRel. Tint in %s:    %f\n', runstr, res.eta_s/(res.radmax^2));
    end
    
    %-- Make figures
    %- Show Donut
    res.figh_donut = figure('color','white','units', 'inches', 'Position', [0 0 svpms.figdimDON]);
    res.figh_donut = VFN_An_fitsDisp(res.frmNrmd, res.xax, res.yax, res.figh_donut);
    %tit.Position(2) = tit.Position(2)*1.02;  % Lift up the title to make cropping easier
    xlabel('Spatial [\lambda/D]')
    ylabel('Spatial [\lambda/D]')
    set(gca,'fontunits', 'normalized', 'fontsize',svpms.fontsize2)
    % Define title as cell array w/ empty second element to lift title off plot slightly
    title({sprintf('%s: %s',runstr,runid), ' '}, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1); 

    %- Show Line Profile
    res.figh_radprof = figure('color','white','units', 'inches', 'Position', [0 0 svpms.figdimLP]);
    plot(res.rax, res.radavg*100, 'LineWidth', svpms.linewidth);
    grid(gca, svpms.gridmodeLP)
    xlabel('Angular Separation [\lambda/D]');
    ylabel('Coupling [%]');
    set(gca, 'fontunits', 'normalized', 'fontsize', svpms.fontsize2)
    % Define title as cell array w/ empty second element to lift title off plot slightly
    title({sprintf('%s: %s',runstr,runid), ' '}, 'FontUnits', 'normalized', 'FontSize', svpms.fontsize1); 
    %tit.Position(2) = tit.Position(2)*1.01;  % Lift up the title to make cropping easier
end