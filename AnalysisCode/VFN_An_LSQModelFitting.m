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

%-- Ouput Folder (where data is saved)
savefld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\MonoWrapUp\NewAnalysis\';

%-- Define analysis parameters (for processing lab data)
% Scaling factor between Femto gain settings
an_params.gnFactor = 9.97;
% Final lens throughput (for red pm normalization)
an_params.lensTHPT = 0.9973;
% Conversion: Red thorlabs PM [uW] to fmto @10^6 [V]
an_params.uw2fmto  = 0.6318;        %P_femto[V]/P_red[uW]
an_params.rdOfmto  = 1/an_params.uw2fmto;    %ratio as defined by Nem (red[uW]/fmto[V])
% System F#
an_params.Fnum  = 10.96/2.1;    % f=10.96 from Thorlabs A397 datasheet
% Wavelength
an_params.lambda = 0.635;   % [microns]
% Conversion: microns to L/D in focal plane [L/D / micron]
an_params.um2LD = 1/(an_params.lambda*an_params.Fnum);                              
% Conversion: piezo [V] to microns
an_params.V2um = 1/10;%8.72;    %<-- Value from PZT characterization (spec = 1/10)
% Fresnel Reflections at fiber surfaces
an_params.FR = 0.9654;  % (1-loss)
% Propagation losses through fiber (fraction per meter)
an_params.PL = 0.9966;  % (1-loss)

%-- Read data filenames
% Read all .fits files in the data directories
nmWpath = { [fldnm '\021619_FNM2\*.fits'];
            [fldnm '\*_FTO3\*.fits']; 
            [fldnm '\*_FTO4\*.fits']};
an_params.STRNMS = VFN_An_getSTRNMS(nmWpath);

%-- Interpolation points for Data vectors
    % Use 50 since PSF has rmax=19 and DON has rmax = 110
datveclength = 50;

%-- Flags for save/print options
isPrintLoc  = false;    % Print null/planet locations
isSave      = false;    % Save (.mat) data
isSaveFig   = false;    % Save figures
isDispItr   = true;     % Print fmincon iteration statuses

%-- Option to use parallel computing for fmincon calculation
isparcomp   = true;

% Default fonts for saved figures
fontsize1 = 14;
fontsize2 = 14;

%% Prep PSF data    (mostly coppied from PlottingScript_Updated)
disp('--- Prepping PSF data')

%-- Get normalized data
nmPSF    = ['021619_FNM2' filesep '23PSF_PeakNewFNM2'];
nrmPSF   = VFN_An_fitsNormed(nmPSF, [1,1,1], an_params);

%-- Calculate peak coupling
[eta_PSF, indPSF] = VFN_An_getEta_p(nrmPSF);
fprintf('\nPeak coupling on PSF [FNM2 #23]:    %f\n',eta_PSF);
if isPrintLoc
    fprintf([' Peak Location: ', ...
             '\n  X index = %3i', ...
             '\n  Y index = %3i', ...
             '\n  Foc ind = %3i', ...
             '\n  VX ind  = %3i', ...
             '\n  VY ind  = %3i\n'], ...
             indPSF(2), indPSF(1), indPSF(4), indPSF(5), indPSF(6));
end

%-- Get axes for image (relative to peak point):
[xaxPSF, yaxPSF]  = VFN_An_getAxes(nmPSF, indPSF, an_params);

%-- Show Frame
figPSF = VFN_An_fitsDisp(nrmPSF, xaxPSF, yaxPSF);
title('PSF Coupling - [FNM2 #23]')

disp('  Getting radial average')

%-- Calculate average radial profile
[radvgPSF,rvecPSF] = VFN_An_radAverage(nrmPSF, [indPSF(2), indPSF(1)], datveclength);

%-- Get radial axis 
% Use xax to get L/D per 'pixel' value then rescale rvec by this value
raxPSF = (xaxPSF(2)-xaxPSF(1))*rvecPSF;

%-- Display radial profile
figPSF_rad = figure('color', 'white');
z = plot(raxPSF,radvgPSF*100, 'LineWidth', 4);
set(gca,'fontsize',fontsize2)
title('Average radial coupling - [FNM2 #23]')
xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
ylabel('Coupling [%]', 'FontSize', fontsize1);


%---- Place variables in requisite structures
datPSF.radvg = radvgPSF;     % Average radial profile of PSF data - vector
datPSF.rax = raxPSF;         % X-axis (radial) for PSF data - vector


%% Prep Donut data      (mostly coppied from PlottingScript_Updated)
disp('--- Prepping Donut data')

%-- Get normalized data
nmDON    = ['021619_FNM2' filesep '21Don_FineFull2'];
nrmDON   = VFN_An_fitsNormed(nmDON, [1,1,1], an_params);

%-- Calculate null coupling
[eta_DON, indDON] = VFN_An_getEta_s(nrmDON);
fprintf('\nNull coupling on Donut [FNM2 #21]:    %f\n',eta_DON);
if isPrintLoc
    fprintf([' Null Location: ', ...
             '\n  X index = %3i', ...
             '\n  Y index = %3i', ...
             '\n  Foc ind = %3i', ...
             '\n  VX ind  = %3i', ...
             '\n  VY ind  = %3i\n'], ...
             indDON(2), indDON(1), indDON(4), indDON(5), indDON(6));
end

%-- Get axes for image (relative to null point):
[xaxDON, yaxDON]  = VFN_An_getAxes(nmDON, indDON, an_params);

%-- Show Frame
fighDON = VFN_An_fitsDisp(nrmDON, xaxDON, yaxDON);
% Use parula for this data since want to see nuances of asymmetry
colormap parula
title('Donut Coupling - [FNM2 #21]')

disp('  Getting radial average')

%-- Calculate average radial profile
[radvgDON,rvecDON] = VFN_An_radAverage(nrmDON, [indDON(2), indDON(1)], datveclength);

%-- Get radial axis 
% Use xax to get L/D per 'pixel' value then rescale rvec by this value
raxDON = (xaxDON(2)-xaxDON(1))*rvecDON;

%-- Print Peak
%Average for the donut ring with the peak coupling
radvgPKDON = max(radvgDON);
fprintf('Peak average [FNM2 #21]: = %f\n', radvgPKDON)
fprintf(' Radius of peak average: = %f\n', raxDON(radvgDON==radvgPKDON))

%-- Display radial profile
fighDON_rad = figure('color', 'white');
z = plot(raxDON,radvgDON*100, 'LineWidth', 4);
ylim([0 20])
set(gca,'fontsize',fontsize2)
title('Average radial coupling - [FNM2 #21]')
xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
ylabel('Coupling [%]', 'FontSize', fontsize1);


%---- Place variables in requisite structures
datDON.radvg = radvgDON;     % Average radial profile of Donut data - vector
datDON.rax = raxDON;         % X-axis (radial) for Donut data - vector

%% Start PSF Model      (mostly coppied from main_pupilBasedVFN)
disp('--- Prepping PSF Model')
N = 2^12; % Size of computational grid (NxN samples) 
apRad = 150; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%-- Create array with pupil function
EP = makeCircularPupil( apRad, N );

%- Plot
% figure(1)
% imagesc(xvals/apRad,yvals/apRad,EP);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));

%-- Get PSF without vortex mask and normalization factors 
% Keeps the origin at the center of the array (keep N even)
myfft2 = @(x) fftshift(fft2(fftshift(x)));
myifft2 = @(x) fftshift(ifft2(fftshift(x)));

PSF = myfft2(EP); % PSF with a flat wavefront 

normI = max(max(abs(PSF).^2));% Normalization for PSF plots
totalPower0 = sum(sum(abs(PSF).^2));% Normalization for coupling fractions
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

%- Plot
% figure(2)
% imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSF).^2/normI);
% axis image; 
% axis([-3 3 -3 3]);
% title('PSF w/o vortex');
% colorbar; 
% colormap(parula(256));

%-- Make vortex mask 
offsetX = 1.6147*apRad; % <-- PSF % 0.0952*apRad; % <-- best null position % 
offsetY = 0.0011*apRad; % <-- PSF % 0.0524*apRad; % <-- best null position % 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

% figure(3)
% imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase');
% colormap(hsv(256));
% colorbar; 

%-- Add Zernike aberration 
disp('  Applying Zernikes to Pupil')

noll_index = 7; % Noll index (Coma - vertical 90*)
coeff7 = -0.0006; %<--VortOut % -0.0008; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = exp(1i*2*pi*coeff7*Z);

noll_index = 8; % Noll index (Coma - horizontal 0*)
coeff8 = 0.0030; %<--VortOut % 0.0033; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff8*Z);

noll_index = 5; % Noll index (Astig - Oblique 45*) 
coeff5 =  -0.0104; %<--VortOut % -0.0099; %<--VortIn % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff5*Z);

noll_index = 6; % Noll index (Astig - Vertical 0*) 
coeff6 =  -0.0095; %<--VortOut % -0.0095; %<--VortIn % %waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff6*Z);

noll_index = 11; % Noll index (Primary Spherical) 
coeff11 = 0;% -0.0194; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff11*Z);

noll_index = 4; % Noll index (Focus) 
coeff4 = 0;%0.01; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff4*Z);

% figure(3)
% imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP.*ABER));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase');
% colormap(hsv(256));
% colorbar; 

%-- Get PSF with vortex mask
PSFv = myfft2(EP.*EPM.*ABER); % PSF with vortex mask and aberration 

% figure(4)
% imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSFv).^2/normI);
% axis image; 
% axis([-3 3 -3 3]);
% title('log10(PSF) w/ vortex');
% colorbar; 
% colormap(parula(256));


%---- Place variables in requisite structures
modPSF.PSFv = PSFv;            % E-field at fiber tip
modPSF.totalPower0 = totalPower0;   % Total power on fibertip
modPSF.lambdaOverD = lambdaOverD;   
modPSF.coords = coords;
modPSF.lambda = an_params.lambda;   % Wavelength


%% Start Donut Model      (mostly coppied from main_pupilBasedVFN)
disp('--- Prepping Donut Model')

N = 2^12; % Size of computational grid (NxN samples) 
apRad = 150; % Aperture radius in samples 
charge = 1; % Charge of the vortex mask 

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%-- Create array with pupil function
EP = makeCircularPupil( apRad, N );

%- Plot
% figure(1)
% imagesc(xvals/apRad,yvals/apRad,EP);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));

%-- Get PSF without vortex mask and normalization factors 
% Keeps the origin at the center of the array (keep N even)
myfft2 = @(x) fftshift(fft2(fftshift(x)));
myifft2 = @(x) fftshift(ifft2(fftshift(x)));

PSF = myfft2(EP); % PSF with a flat wavefront 

normI = max(max(abs(PSF).^2));% Normalization for PSF plots
totalPower0 = sum(sum(abs(PSF).^2));% Normalization for coupling fractions
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

%- Plot
% figure(2)
% imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSF).^2/normI);
% axis image; 
% axis([-3 3 -3 3]);
% title('PSF w/o vortex');
% colorbar; 
% colormap(parula(256));

%-- Make vortex mask 
offsetX = 0.0952*apRad; % <-- best null position % 1.6147*apRad; % <-- PSF % 
offsetY = 0.0524*apRad; % <-- best null position % 0.0011*apRad; % <-- PSF % 

EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

% figure(3)
% imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase');
% colormap(hsv(256));
% colorbar; 

%-- Add Zernike aberration 
disp('  Applying Zernikes to Pupil')

noll_index = 7; % Noll index (Coma - vertical 90*)
coeff7 = -0.0008; %<--VortIn % -0.0006; %<--VortOut % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = exp(1i*2*pi*coeff7*Z);

noll_index = 8; % Noll index (Coma - horizontal 0*)
coeff8 = 0.0033; %<--VortIn % 0.0030; %<--VortOut % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff8*Z);

noll_index = 5; % Noll index (Astig - Oblique 45*) 
coeff5 =  -0.0099; %<--VortIn % -0.0104; %<--VortOut % % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff5*Z);

noll_index = 6; % Noll index (Astig - Vertical 0*) 
coeff6 =  -0.0095; %<--VortIn % -0.0095; %<--VortOut % %waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff6*Z);

noll_index = 11; % Noll index (Primary Spherical) 
coeff11 = 0;% -0.0194; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff11*Z);

noll_index = 4; % Noll index (Focus) 
coeff4 = 0;%0.01; % waves rms
[Z,n,m] = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
Z = Z/sqrt(mean(Z(logical(EP)).^2)); % Re-normalize (useful when pupil is not a circle)
ABER = ABER.*exp(1i*2*pi*coeff4*Z);

% figure(3)
% imagesc(xvals/apRad,yvals/apRad,angle(EPM.*EP.*ABER));
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil phase');
% colormap(hsv(256));
% colorbar; 

%-- Get PSF with vortex mask
PSFv = myfft2(EP.*EPM.*ABER); % PSF with vortex mask and aberration 

% figure(4)
% imagesc(xvals/lambdaOverD,yvals/lambdaOverD,abs(PSFv).^2/normI);
% axis image; 
% axis([-3 3 -3 3]);
% title('log10(PSF) w/ vortex');
% colorbar; 
% colormap(parula(256));


%---- Place variables in requisite structures
modDON.PSFv = PSFv;            % E-field at fiber tip
modDON.totalPower0 = totalPower0;   % Total power on fibertip
modDON.lambdaOverD = lambdaOverD;   
modDON.coords = coords;
modDON.lambda = an_params.lambda;   % Wavelength

%% Use fmincon to do LSQ fit
disp('--- Doing LSQ minimization')

%-- Define function LSQ calculating function
func = @(X) VFN_An_LSQfunc(X, datPSF, datDON, modPSF, modDON);

%-- Define starting parameters
iMFD    = 1.4*an_params.lambda*an_params.Fnum;  % ideal MFD for F# and lambda
X0  = [1;                       % thpt: Start without scaling throughput
       1;                       % pztg: Start without scaling radial axis
       iMFD/an_params.Fnum];    % mfdf: Start with ideal MFD and F# for system

%-- Provide no constraints
A   = [];
b   = [];
Aeq = [];
beq = [];

%-- Provide lower and upper bounds
% No negative scaling values or MFD/F#;
lb = [0.5;                      % min throughput scaling of 1/2           
      1/an_params.V2um/13.5;    % max piezocal of 13.5 V/micron
      1.5/an_params.Fnum];      % min MFD of 1.5 microns
% Provide no upper limits
ub = [2;                    % max throughput scaling of 2
      1/an_params.V2um/7;   % min piezocal of 7 V/micron
      8/an_params.Fnum];    % max MFD of 8 microns

%-- Set special minimization parameters
% Set nonlcon blank so that we can provide 'options'
nonlcon = [];
% Set options:
    % Display, iters = display status at each iterations
    % StepTolerance = minimum change in X before exit
    % MaxIterations = maximum number of iterations before exit
    % UseParallel = Numerically calculate the gradient using parallel computing
if isDispItr
    disptype = 'iter';  % Status at each iteration
else
    disptype = 'final'; % Status only at end
end
options = optimoptions('fmincon', ...
                       'Display', disptype, ...
                       'StepTolerance', 1e-8, ... %2e-3, ...
                       'MaxIterations', 750, ... % 50, ...
                       'UseParallel', isparcomp);

%-- Minimize LSQ
disp('  Performing minimization')
[X, fval, exitflag, output] = fmincon(func, X0, A, b, Aeq, beq, lb, ub, nonlcon, options);
disp('  Minimization complete')

%% Process results
disp('--- Processing results')

%-- Extract results from X vector
    % NOTE: X-vector the produces results within spec but not best LSQval:
    %   X = [0.8780, 1.0731, 0.9146]; %<-- Original solution w/ close-to-spec MFD
    %   X = [0.920 , 1.085 , 1.0669]; %<-- Solution with lowest LSQ value
thpt = X(1);        % Throughput scaling factor
pztg = X(2);        % piezo-gain (radial) scaling factor
mfdf = X(3);        % MFD/F# ratio

%-- Apply new piezo gain
an_params.V2um = an_params.V2um*pztg;

%-- Calculate new MFD
MFD = an_params.Fnum*mfdf;

%-- Print values in meaningful terms
fprintf('\n  Throughput correction:  %0.3f\n', thpt);
fprintf(  '  Piezo Calibration:      %0.3f [um/V]\n',1/an_params.V2um);
fprintf(  '  Fiber MFD:              %0.3f [um]\n', MFD);

%% Re-create models with new MFD value
disp('--- Generating models with new values')

%-- Generate fibermode using correct MFD
% Fiberdiam (MFD) [in lambda/D]
fibD_psf = MFD*an_params.um2LD;
fibD_don = MFD*an_params.um2LD;
% Fiberdiam (MFD) [in samples of model]
fibD_psf = fibD_psf*modPSF.lambdaOverD;
fibD_don = fibD_don*modDON.lambdaOverD;
% Fiber mode
fibmod_psf = generateSMFmode_gaussian(fibD_psf, modPSF.coords);
fibmod_don = generateSMFmode_gaussian(fibD_don, modDON.coords);

disp('  Calculating coupling maps')

%-- Generate coupling maps
coup_psf = generateCouplingMap(fibmod_psf, modPSF.PSFv, modPSF.totalPower0, 8*modPSF.lambdaOverD, modPSF.coords);
coup_don = generateCouplingMap(fibmod_don, modDON.PSFv, modDON.totalPower0, 8*modDON.lambdaOverD, modDON.coords);

%-- Get peak (PSF) and null (Donut)
[etaP_psf, ind_psf] = VFN_An_getEta_p(coup_psf);
[etaS_psf, ind_don] = VFN_An_getEta_s(coup_don,2.03);

%-- Get average radial profiles
cent = [ind_psf(2), ind_psf(1)];
% Interpolate at 2x model sampling
[radvg_psf,rvec_psf] = VFN_An_radAverage(coup_psf, cent, modPSF.coords.N*2);
cent = [ind_don(2), ind_don(1)];
% Interpolate at 2x model sampling
[radvg_don,rvec_don] = VFN_An_radAverage(coup_don, cent, modDON.coords.N*2);

%-- Get radial axes
rax_psf = rvec_psf/modPSF.lambdaOverD;
rax_don = rvec_don/modDON.lambdaOverD;

%% Create ideal models (F# matched to MFD from fit)
disp('--- Generating ideal F# models')

%-- Generate fibermode using correct MFD
% Ideal fiberdiam (MFD) [in lambda/D]
fibD_psfi = 1.4;
fibD_doni = 1.4;
% Fiberdiam (MFD) [in samples of model]
fibD_psfi = fibD_psfi*modPSF.lambdaOverD;
fibD_doni = fibD_doni*modDON.lambdaOverD;
% Fiber mode
fibmod_psfi = generateSMFmode_gaussian(fibD_psfi, modPSF.coords);
fibmod_doni = generateSMFmode_gaussian(fibD_doni, modDON.coords);

disp('  Calculating coupling maps')

%-- Generate coupling maps
coup_psfi = generateCouplingMap(fibmod_psfi, modPSF.PSFv, modPSF.totalPower0, 8*modPSF.lambdaOverD, modPSF.coords);
coup_doni = generateCouplingMap(fibmod_doni, modDON.PSFv, modDON.totalPower0, 8*modDON.lambdaOverD, modDON.coords);

%-- Get peak (PSF) and null (Donut)
[etaP_psfi, ind_psfi] = VFN_An_getEta_p(coup_psfi);
[etaS_psfi, ind_doni] = VFN_An_getEta_s(coup_doni,2.03);

%-- Get average radial profiles
cent = [ind_psfi(2), ind_psfi(1)];
% Interpolate at 2x model sampling
[radvg_psfi,rvec_psfi] = VFN_An_radAverage(coup_psfi, cent, modPSF.coords.N*2);
cent = [ind_doni(2), ind_doni(1)];
% Interpolate at 2x model sampling
[radvg_doni,rvec_doni] = VFN_An_radAverage(coup_doni, cent, modDON.coords.N*2);

%-- Get radial axes
rax_psfi = rvec_psfi/modPSF.lambdaOverD;
rax_doni = rvec_doni/modDON.lambdaOverD;

%% Re-Analyze Lab Data with new Parameters
disp('--- Re-interpolating Lab Data')

%------------------------ PSF -------------------------------

%-- Get axes for image (relative to peak point) [this includes new V2um]:
[xaxPSF, yaxPSF]  = VFN_An_getAxes(nmPSF, indPSF, an_params);

%-- Re-scale throughput by new value
nrmPSFSc = nrmPSF/thpt;

disp('  Getting radial average')

%-- Calculate average radial profile
    % no radpts since want optimal sampling
[radvgPSF,rvecPSF] = VFN_An_radAverage(nrmPSFSc, [indPSF(2), indPSF(1)]);

%-- Get radial axis 
% Use xax to get L/D per 'pixel' value then rescale rvec by this value
raxPSF = (xaxPSF(2)-xaxPSF(1))*rvecPSF;

%------------------------ DONUT -------------------------------

%-- Get axes for image (relative to null point) [this includes new V2um]:
[xaxDON, yaxDON]  = VFN_An_getAxes(nmDON, indDON, an_params);

%-- Re-scale throughput by new value
nrmDONSc = nrmDON/thpt;

disp('  Getting radial average')

%-- Calculate average radial profile
    % no radpts since want optimal sampling
[radvgDON,rvecDON] = VFN_An_radAverage(nrmDONSc, [indDON(2), indDON(1)]);

%-- Get radial axis 
% Use xax to get L/D per 'pixel' value then rescale rvec by this value
raxDON = (xaxDON(2)-xaxDON(1))*rvecDON;

%% Plot Results
nmtag = sprintf('MFD%05.2f_PZTG%05.2f_THPT%05.2f',MFD, 1/an_params.V2um, thpt);

%-- Display 2D Map of Lab Donut with corrections
fighDON = VFN_An_fitsDisp(nrmDONSc, xaxDON, yaxDON);
% Use parula for this data since want to see nuances of asymmetry
colormap parula
title('Donut Coupling (Lab - Rescaled)')
if isSaveFig
    export_fig([savefld filesep 'ModelFit_' nmtag '_LabDonutCoupMap.png'],'-r300', '-painters')
end

%-- Display 2D Map of Model Donut that matches Lab
figure('color', 'white');
xax_sim = ((1:size(coup_don,2))-ind_don(2))/modDON.lambdaOverD;
yax_sim = ((1:size(coup_don,1))-ind_don(1))/modDON.lambdaOverD;
imagesc(xax_sim, yax_sim, coup_don); 
axis([-1.3 1.3 -1.3 1.3]); 
% Use parula for this data since want to see nuances of asymmetry 
colormap parula 
title('Donut Coupling (Sim - Rescaled)')  
if isSaveFig
    export_fig([savefld filesep 'ModelFit_' nmtag '_ModelDonutCoupMap.png'],'-r300', '-painters')
end

%-- Display PSF profile
figPSFTot_rad = figure('color', 'white');
% Plot Corrected lab data
z = plot(raxPSF,radvgPSF*100, 'LineWidth', 4);
hold on
% Plot matching model
plot(rax_psf, radvg_psf*100, 'LineWidth', 4);
% Plot ideal model
plot(rax_psfi, radvg_psfi*100, 'LineWidth', 4);
set(gca,'fontsize',fontsize2)
xlim([0 1.5])
title('Average radial PSF coupling')
xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
ylabel('Coupling [%]', 'FontSize', fontsize1);
legend('Lab', 'Model', 'Ideal');
if isSaveFig
    export_fig([savefld filesep 'ModelFit_' nmtag '_PSFRadProfile.png'],'-r300', '-painters')
end

%-- Display Donut profile
figDONTot_rad = figure('color', 'white');
% Plot Corrected lab data
z = plot(raxDON,radvgDON*100, 'LineWidth', 4);
hold on
% Plot matching model
plot(rax_don, radvg_don*100, 'LineWidth', 4);
% Plot ideal model
plot(rax_doni, radvg_doni*100, 'LineWidth', 4);
set(gca,'fontsize',fontsize2)
xlim([0 1.5])
title('Average radial Donut coupling')
xlabel('Angular separation [\lambda/D]', 'FontSize', fontsize1);
ylabel('Coupling [%]', 'FontSize', fontsize1);
legend('Lab', 'Model', 'Ideal');
if isSaveFig
    export_fig([savefld filesep 'ModelFit_' nmtag '_DonutRadProfile.png'],'-r300', '-painters')
end