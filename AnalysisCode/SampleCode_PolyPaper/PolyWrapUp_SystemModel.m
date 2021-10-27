%{
Script for modeling the VFN Polychromatic system performance

This is based off the main_PupilBasedVFN_Polychromatic script. 
It uses images of the actual Pupil and HASO WFE measurements to input the
true system conditions. Images of the vortex centering in the pupil are
also used to determine the vortex offsets.
%}

clear; close all; 
addpath('../AnalysisLib')
addpath('../../../VFN-Simulations/VFNlib/');
addpath('../../../falco-matlab/lib/utils')

%% Input parameters 

% Define sampling info
N = 2^12; % Size of computational grid (NxN samples) 
apRad = 256/2; % Aperture radius in samples 

% Define wavelength info
lambda0 = 650e-9; %central wavelength
fracBW = 1e-9; %100/650;% \Delta\lambda/\lambda
numWavelengths = 1;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = lambda0./lambdas; % simple scalar model

% ############# Pupil Settings
% Two options:

% ------ Option 1: Use auto-generated circular pupil
% Nothing is needed for this option

% ------ Option 2: Use measured pupil from image
% Folder containing the pupil images
fldnm = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/20210726_PupilImages_PolyPaper';
% Image file to use (asterisk wildcard allowed
  % NOTE: must be a 1x1 cell array 
nmWpath = {fullfile(fldnm,'4_*.fits')};
% Threshold to apply to image for identifying the pupil
thres = 10;
% Kernel size for gaussian blurring applied to the edge of the pupil
gaussBlurSig = 0.3;

%-- Choose an option to use
pupmodel = 2; % scalar 1 or 2 matching option # above
% ############# Pupil Settings

% ############# WFE Settings
% Two options for defining the WFE:

% ------ Option 1: Create zernikes based on provided coeffs
%--Provide Zernike coefficients for WF reconstruction
  % Define wavefront error at the central wavelength
nolls = [6,5, 13,12, 24,23, 39,38];
%coeffs = [0.0029,0.0070, -0.0046,0.0004, 0.0015,0.0004, 0*-0.0001,0*-0.0003]*(635*1e-9)/lambda0;
coeffs = [5.9e-4,0.0014, -7.2732e-04,6.3246e-05, 2.0045e-04,5.3452e-05, 0*-0.0001,0*-0.0003]*(635*1e-9)/lambda0;

% ------ Option 2: Load a WF map from HASO
%-- Filename for WFE map from HASO
%flnm_wfe = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/HASO_Measurements/4_IrisFullyOpen_FinalAlignment_200mmlens_WF_TTFocSpheSubt.txt';
flnm_wfe = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/20210726_HasoData_PolyPaper/1_After200mmLens_VortexOffCenter_WF_TTFocSpherSubt.txt';
%-- Add some defocus on top of the HASO measurements
  % Set to 0 if you don't want any additional defocus
defoc = 0; 

%-- Choose an option to use
wfemodel = 2;  % scalar 1 or 2 matching option # above
% ############# WFE Settings

% Give offsets for the vortex mask
offsetX = 0.0118*apRad*2; %PSF = 2.54pupdiam --> 214*0.0118
offsetY = 0.0256*apRad*2;

% ############# Fiber and F/# settings
% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
    % Core diam from: https://www.fiberoptics4sale.com/products/sm600
fiber_props.core_rad = 3.65e-6/2;% Core radius [um] Measured on Microscope --> 3.48e-6/2
fiber_props.n_core = 1.460095;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.45667;% cladding index (interpolated from linear fit to 3 points)
fiber_props.type = 'bessel';

%-- Choose which which fiber MFD to use for computing the fibermode:
  % 'MFDComputed' will use the the fiber_props to generate the fibermode
  % 'MFDManual' will use the user-provided 'MFD' to create a gaussian fibermode
modemodel = 'MFDComputed';

%-- MFD for SM600 (from Thorlabs website)
MFD = 1.12*4.45e-6;%mean([3.6,5.3])*1e-6;      % [m] 3.6-5.3

%-- Provide actual system Fnum (this is what will be used in the sim)
Fnum = 11.0e-3/2.45e-3;

%-- Provide fiber tip x and y tilt
  % I empirically determined that based on the equation I am using to apply
  % the fiber-tip tilt, a value of 8 below corresponds to a pupil
  % translation of 1 pupil diameter for the fibermode.
tltX = 1/4.8825;    
tltY = 1/10.6;

% ############# Fiber and F/# settings

%-- Saving info
isSave = false; % flag choosing whether to save a .csv with the central-band coupling profile
% folder into which to save the .csv
svfld = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/PublishedInPaper';

% FILENAME to use for the saved .csv file
  % set blank ( [] ) to auto-generate the name
%svnm   = 'VFN_PolyPaper_ModeledCoupProf_4.49Fnum_ModeManuGauss3.6umMFD.csv'; 
svnm = [];

%-- Default fonts for saved figures
% Fontsizes are in 'normalized' units. That means the value should be
  % fractional [0 to 1]. Ex: 0.1 sets text to 10% of the figure height.
fontsize1 = 0.07;     % Title fontsize
fontsize2 = 0.05;     % Other text fontsize

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

if pupmodel == 1
    PUPIL = makeCircularPupil( apRad, N );
elseif pupmodel == 2
    %-- Load the data
    % Find all samples for the given set
    STRNMS = VFN_An_getSTRNMS(nmWpath);
    % Preallocate matrix 
    PUPIL = nan([fitsinfo(STRNMS(1)).PrimaryData.Size, length(STRNMS)]);
    % Iterate through samples loading each one
    for ind = 1:length(STRNMS)
        PUPIL(:,:,ind) = fitsread(STRNMS(1));
    end

    %-- Average samples
    PUPIL = mean(PUPIL,3);
    %figure(); imagesc(pupim); axis image

    %-- Set areas outside pupil to 0 and inside pupil to 1
    PUPIL = PUPIL >= thres;
    %figure(); imagesc(pupim); axis image

    %-- Sharpen Edges
    windowSize = 20;
    fprintf('\n-- NOTE: using windowSize = %0.1f for edge sharpening on PUPIL\n',windowSize)
    kernel = ones(windowSize) / windowSize ^ 2;
    PUPIL = conv2(single(PUPIL), kernel, 'same');
    PUPIL = PUPIL > 0.5; % Rethreshold
    
    %-- Crop image so it just barely fits the pupil
    % Auto-detect the pupil boundaries
    [~, row0] = find(PUPIL',1);
    [~, col0] = find(PUPIL,1);
    [~, row1] = find(PUPIL',1,'last');
    [~, col1] = find(PUPIL,1,'last');
    % Crop 
    PUPIL = PUPIL(row0:row1,col0:col1);
    % Make array square
    PUPIL = pad_crop(PUPIL,max(size(PUPIL)));
    
    %-- Resample to match desired apRad
    % First rotate the pupil so that one pair of flat edges aligns with the
    % matrix dimensions
    rotang = 22.5;
    fprintf('\n-- NOTE: rotating PUPIL by %0.2fdeg to get edges aligned with matrix edges\n',rotang)
    PUPIL = imrotate(PUPIL, rotang);
    % Crop to just fit the pupil again
      % NOTE: the previous crop was needed for centering the pupil before
      % rotation, this crop is for getting the sampling/size right
    [~, row0] = find(PUPIL',1);
    [~, col0] = find(PUPIL,1);
    [~, row1] = find(PUPIL',1,'last');
    [~, col1] = find(PUPIL,1,'last');
    % Crop 
    PUPIL = PUPIL(row0:row1,col0:col1);
    % Make array square
    PUPIL = pad_crop(PUPIL,max(size(PUPIL)));
    % Display the pre-resampled and pre-gaussfiltered pupil for reference
    figure(100); imagesc(PUPIL); axis image; title('Cleaned PUPIL before Resampling and Edge Smoothing');
    % Resample to desired size
    [X,Y] = meshgrid(linspace(1,size(PUPIL,1),apRad*2),linspace(1,size(PUPIL,2),apRad*2));
    PUPIL = interp2(PUPIL,X,Y);
    % Smoothen the edges slighty to remove sharp diffraction effects
      % Use simple gaussian blurring 
    fprintf('\n-- NOTE: using gaussian blur kernel of sigma=%0.4f\n',gaussBlurSig)
    PUPIL = imgaussfilt(PUPIL,gaussBlurSig);
    
    %-- Rotate back to original orientation
    PUPIL = imrotate(PUPIL, -rotang);

    %-- Resize to match simulation matrices
    PUPIL = pad_crop(PUPIL, N);
    
    %-- Transpose to match system orientation
      % The camera was clocked 90deg w.r.t the HASO and fiber axes
    PUPIL = PUPIL';
else
    error('pupmodel was not recognized')
end

[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

figure(1)
imagesc(xvals/apRad,yvals/apRad,PUPIL);
axis image; 
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Define pupil field

if wfemodel == 1
    phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coords); 
elseif wfemodel == 2
    %-- Load measured WF from HASO
    % Read the file
    HASO_WFE = readmatrix(flnm_wfe);
    % Replace nans with 0
    HASO_WFE(isnan(HASO_WFE)) = 0;
    
    %-- Crop image around the WF data
    % Auto-detect the boundaries
    [~, row0] = find(HASO_WFE',1);
    [~, col0] = find(HASO_WFE,1);
    [~, row1] = find(HASO_WFE',1,'last');
    [~, col1] = find(HASO_WFE,1,'last');
    % Crop 
    HASO_WFE = HASO_WFE(row0:row1,col0:col1);
    % Make array square
    HASO_WFE = pad_crop(HASO_WFE,max(size(HASO_WFE)));
    % Crop to pertinent region
    %HASO_WFE = HASO_WFE(1:19,24:41);
    
    %-- Interpolate to the pupil size
    wfpad = 0;     % Add a little bit of padding around the edges 
    [X,Y] = meshgrid(linspace(1,size(HASO_WFE,1),apRad*2+wfpad),linspace(1,size(HASO_WFE,2),apRad*2+wfpad));
    HASO_WFE = interp2(HASO_WFE,X,Y);
    % Pad to size of pupil image in simulation
    HASO_WFE = pad_crop(HASO_WFE,N);
    
    %-- Rescale from waves @635 to 650nm
    HASO_WFE = HASO_WFE*(635*1e-9)/lambda0;
     
    %-- Add user-specified defocus
    phz = HASO_WFE + (1/(2*pi))*generateZernike_fromList(4, defoc, PUPIL, apRad, coords);
    
     %-- Convert from waves to phase units
    phz = 2*pi*phz;
else
    error('wfemodel was not recognized')
end

figure(2);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask

% Get broadband PSF
iPSF_BB = mean(getPSF(Epup,lambda0,lambdas,normI,coords),3);

figure(3)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSF_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(hsv(256));
colorbar; 
drawnow;

%% Get PSF with vortex mask

iPSFv_BB = mean(getPSF(Epup.*EPM,lambda0,lambdas,normI,coords),3);

figure(5)
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/ vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate coupling maps for each wavelength
%-- Compute MFD from fiber_props above
MFD_FibProp = getMFD(fiber_props,lambda0);
fprintf('MFD from fiber_props: %f [um]\n',MFD_FibProp*1e6);

%-- Compute ideal F/# assuming the MFD from fiber_props above
Fnum_ideal_FipProp = MFD_FibProp/(lambda0*1.4); % focal ratio of the beam at the fiber
fprintf('Ideal F/# from fiber_props: %f\n',Fnum_ideal_FipProp);

%-- Compute ideal F/# using the input MFD provided above
Fnum_ideal_MFD = MFD/(lambda0*1.4); 
fprintf('Ideal F/# from input MFD: %f\n',Fnum_ideal_MFD);

% Preallocate fibmode matrix
fibmode = nan(N,N,numWavelengths);
for ch = 1:numWavelengths
    %----- ############ TWO OPTIONS FOR GETTING ETA_MAPS: 
    if strcmpi(modemodel, 'MFDManual')
        %----- Option 1: Compute fibermode using the input MFD 
            % (assuming constant MFD at all wavelengths)
        disp('Using fibermode based on input MFD')
        fibmode(:,:,ch) = generateSMFmode_gaussian((MFD/(lambdas(ch)*Fnum))*lambdaOverD,coords);
    elseif strcmpi(modemodel, 'MFDComputed')
        %----- Option 2: Compute fibermode using the MFD from the fiber_props
            % (Same method that generateCouplingMap_polychromatic would normally use)
        disp('Using fibermode based on fiber_props')
        fibmode(:,:,ch) = generateSMFmode( fiber_props, lambdas(ch), Fnum, lambdaOverD, coords);
    else 
        error('modemodel was not valid')
    end
    
    %-- Apply fiber tip/tilt 
    fibmode(:,:,ch) = fibmode(:,:,ch).*exp(1i*2*pi*coords.X*tltX*(lambda0/lambdas(ch))*lambdaOverD/N);
    fibmode(:,:,ch) = fibmode(:,:,ch).*exp(1i*2*pi*coords.Y*tltY*(lambda0/lambdas(ch))*lambdaOverD/N);
end

eta_maps = generateCouplingMap_polychromatic( Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 5*lambdaOverD, coords, fibmode);

figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end

%% Do a "mock" analysis using method from testbed data

% Crop the coupling map to not show massive FOV unnecessarily
crpW = 3;  %[lambda/D] FOV (radius wr.t center of image)
xcrop = xvals(abs(xvals/lambdaOverD) < crpW)+N/2+1;
ycrop = yvals(abs(yvals/lambdaOverD) < crpW)+N/2+1;
xvals_cropped = xvals(xcrop)/lambdaOverD; yvals_cropped = yvals(ycrop)/lambdaOverD;

for ch = 1:numWavelengths
    results(ch).donmap = eta_maps(xcrop,ycrop,ch);
    %figure(); imagesc(results(ch).donmap); axis image; 

    % Find/Measure null
    [results(ch).eta_s, results(ch).eta_s_ind] = VFN_An_getEta_s(results(ch).donmap,7);

    % Find/Measure peak (overall max)
    [results(ch).eta_max, results(ch).eta_max_ind] = VFN_An_getEta_p(results(ch).donmap);

    % Compute radial average (centered on null position)
    [results(ch).radavg, results(ch).rvec] = VFN_An_radAverage(results(ch).donmap, [results(ch).eta_s_ind(2), results(ch).eta_s_ind(1)]);
    % Get planet coupling as peak in radial average
    [results(ch).radmax, results(ch).rad_max_ind] = max(results(ch).radavg);

    % Get radial average axes
    results(ch).rax = results(ch).rvec/lambdaOverD;

    %-- Print results
    % Null and location
    fprintf('Null: %e\n',results(ch).eta_s)
    fprintf([' eta_s Location: ', ...
             '\n  X index = %3i', ...
             '\n  Y index = %3i', ...
             '\n  Foc ind = %3i', ...
             '\n  VX ind  = %3i', ...
             '\n  VY ind  = %3i\n'], ...
             results(ch).eta_s_ind(2), results(ch).eta_s_ind(1), results(ch).eta_s_ind(4), ...
             results(ch).eta_s_ind(5), results(ch).eta_s_ind(6));

    % Peak and location
    fprintf('Max:  %f\n',results(ch).eta_max)
    fprintf([' eta_p Location: ', ...
             '\n  X index = %3i', ...
             '\n  Y index = %3i', ...
             '\n  Foc ind = %3i', ...
             '\n  VX ind  = %3i', ...
             '\n  VY ind  = %3i\n'], ...
             results(ch).eta_max_ind(2), results(ch).eta_max_ind(1), results(ch).eta_max_ind(4), ...
             results(ch).eta_max_ind(5), results(ch).eta_max_ind(6));
    fprintf('\nPeak RadAvg:    %f\n', results(ch).radmax);
    fprintf(' Peak RadAvg at %f lam/D\n', results(ch).rax(results(ch).rad_max_ind));
    %-- Print relative integration time
    fprintf('\nRel. Tint:    %f\n', results(ch).eta_s/(results(ch).radmax^2));

    %-- Make figures
    %- Show Donut
    results(ch).figh_donut = figure(20+ch);
    set(results(ch).figh_donut, 'color','white','units', 'inches', 'Position', [0 0 7, 6]);
    results(ch).figh_donut = VFN_An_fitsDisp(results(ch).donmap, xvals_cropped, yvals_cropped, results(ch).figh_donut);
    axis([-2.5, 2.5 -2.5 2.5])
    %tit.Position(2) = tit.Position(2)*1.02;  % Lift up the title to make cropping easier
    xlabel('Spatial [\lambda/D]')
    ylabel('Spatial [\lambda/D]')
    set(gca,'fontunits', 'normalized', 'fontsize',fontsize2)
    % Define title as cell array w/ empty second element to lift title off plot slightly
    title({[num2str(lambdas(ch)*1e9,'%0.0f') 'nm Coupling Map'], [num2str(fracBW,'%0.2f') '\Delta\lambda/\lambda_0']}, 'FontUnits', 'normalized', 'FontSize', fontsize1); 

    %- Show Line Profile
    results(ch).figh_radprof = figure(40+ch);
    set(results(ch).figh_radprof, 'color','white','units', 'inches', 'Position', [0 0 7 6]);
    plot(results(ch).rax, results(ch).radavg*100, 'LineWidth', 2.3);
    xlim([0, 2.5])
    grid on;
    xlabel('Angular Separation [\lambda/D]');
    ylabel('Coupling [%]');
    set(gca, 'fontunits', 'normalized', 'fontsize', fontsize2)
    % Define title as cell array w/ empty second element to lift title off plot slightly
    title({[num2str(lambdas(ch)*1e9,'%0.0f') 'nm Line Profile'], [num2str(fracBW,'%0.2f') '\Delta\lambda/\lambda_0']}, 'FontUnits', 'normalized', 'FontSize', fontsize1); 
    
end
    
%% (OPTIONAL) Save the Line Profile (radavg) for the central wavelengths

if isSave
    %-- Make a 2-column matrix containing the position in col1 and coup in col2
    coupmat = [results(central_band_index).rax',results(central_band_index).radavg];
    
    %-- Format the filename
    if isempty(svnm)
        % No filename was provided so generate one automatically
        
        prefix = sprintf('VFN_PolyPaper_ModeledCoupProf_%0.2fFnum_', Fnum);
        
        if strcmpi(modemodel, 'MFDManual')
            tag = [modemodel '_gaussian_' num2str(MFD*1e6,'%0.2fumMFD')];
        else
            tag = [modemodel '_' fiber_props.type '_' num2str(fiber_props.core_rad*1e6,'%0.2fumCoreRad')];
        end
        svnmFull = fullfile(svfld,[prefix tag '.csv']);
    else
        svnmFull = fullfile(svfld,svnm);
    end

    %-- Save the matrix as an easy-to-load csv
    writematrix(coupmat,svnmFull);
    
    fprintf('\nSaved line profile to:\n  %s\n',svnmFull);
end

%% Analyze the full polychromatic result
%-- First, combine the coupling map at all the wavelenths
  % This can be done as an average since all maps are normalized by the
  % same total power. This average is the same as if we'd summed all the
  % powers through the fiber and divided by all the incident power.
donmapBB = nan([size(results(1).donmap),numWavelengths]);
for ch = 1:numWavelengths
    donmapBB(:,:,ch) = results(ch).donmap;
end
donmapBB = mean(donmapBB,3);

% Find/Measure null
[eta_s, eta_s_ind] = VFN_An_getEta_s(donmapBB,7);

% Find/Measure peak (overall max)
[eta_max, eta_max_ind] = VFN_An_getEta_p(donmapBB);

% Compute radial average (centered on null position)
[radavg, rvec] = VFN_An_radAverage(donmapBB, [eta_s_ind(2), eta_s_ind(1)]);
% Get planet coupling as peak in radial average
[radmax, rad_max_ind] = max(radavg);

% Get radial average axes
rax = rvec/lambdaOverD;

%-- Print results
% Null
fprintf('Null: %e\n', eta_s)
% Peak
fprintf('Max:  %f\n',eta_max)
fprintf('\nPeak RadAvg:    %f\n', radmax);
fprintf(' Peak RadAvg at %f lam/D\n', rax(rad_max_ind));


%-- Make figures
%- Show Donut
figh_donut = figure(50);
set(figh_donut, 'color','white','units', 'inches', 'Position', [0 0 7, 6]);
figh_donut = VFN_An_fitsDisp(donmapBB, xvals_cropped, yvals_cropped, figh_donut);
axis([-2.5, 2.5 -2.5 2.5])
%tit.Position(2) = tit.Position(2)*1.02;  % Lift up the title to make cropping easier
xlabel('Spatial [\lambda/D]')
ylabel('Spatial [\lambda/D]')
set(gca,'fontunits', 'normalized', 'fontsize',fontsize2)
% Define title as cell array w/ empty second element to lift title off plot slightly
title({'Broadband Coupling Map', [num2str(fracBW,'%0.2f') '\Delta\lambda/\lambda_0']}, 'FontUnits', 'normalized', 'FontSize', fontsize1); 

%- Show Line Profile
figh_radprof = figure(51);
set(figh_radprof, 'color','white','units', 'inches', 'Position', [0 0 7 6]);
plot(rax, radavg*100, 'LineWidth', 2.3);
xlim([0, 2.5])
grid on;
xlabel('Angular Separation [\lambda/D]');
ylabel('Coupling [%]');
set(gca, 'fontunits', 'normalized', 'fontsize', fontsize2)
% Define title as cell array w/ empty second element to lift title off plot slightly
title({'Broadband Line Profile', [num2str(fracBW,'%0.2f') '\Delta\lambda/\lambda_0']}, 'FontUnits', 'normalized', 'FontSize', fontsize1); 
    