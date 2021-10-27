%{
%}

clear; close all; 
addpath(genpath(fullfile('..','..','..','VFN-Simulations/VFNlib/')));
addpath(genpath(fullfile('..','..','..','falco-matlab')))

%% Input parameters 

%-- Provide regular parameters
% Define smapling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

%-- Define wavelength info
lambda0 = 650e-9; %central wavelength
fracBW = 0.2;%100/650;% \Delta\lambda/\lambda
numWavelengths = 3;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Parameters for Thorlabs SM600
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 3.65e-6/2;% Core radius [um]
fiber_props.n_core = 1.460095;% core index (interp @650nm)
fiber_props.n_clad = 1.45667;% cladding index (interp @650nm)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.45); %ideal for bess: 1.238 % focal ratio of the beam at the fiber

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 2.45e-3;    %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
foc = Fnum*DPup;
fprintf('Focus in use: %f [mm]\n',foc*1e3)

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 20; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 15;    %[lam0/D] field of view in final image plane

%-- Define charge of the vortex mask at each wavelength
charge = 0*ones(1,numWavelengths); % achromatic

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvalsPP = coordsPP.xvals;% Helpful for plotting
yvalsPP = coordsPP.yvals;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
coordsFP = generateCoordinates( Nxi );
xvalsFP = coordsFP.xvals;% Helpful for plotting
yvalsFP = coordsFP.yvals;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in radians at final focal plane 
lambda0Fnum_meters = lambda0*foc/DPup;     
% "pixel" size in pupil and image planes
dx  = DPup/(2*apRad);   % meters/sample
dxi = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample


%% Create array with pupil function

PUPIL = makeCircularPupil( apRad, N );

%-- Get normalization factors
% Norm for coupling fractions (simple sum since using MFT propagator)
totalPower0 = sum(abs(PUPIL(:)));
% Norm for PSF plots (peak of ideal PSF)
PSF = propcustom_mft_PtoF(PUPIL, foc, lambda0, dx, dxi, Nxi, dxi, Nxi);
normI = max(abs(PSF(:)).^2);

figure(1); 
imagesc(xvalsPP/apRad,yvalsPP/apRad,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Define pupil field

nolls = 2;
coeffs = 0;
phz = generateZernike_fromList( nolls, coeffs, PUPIL, apRad, coordsPP); 

figure(2);
Epup = nan(N,N,numWavelengths);
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    subplot(1,numWavelengths,ch);
    imagesc(xvalsPP/apRad,yvalsPP/apRad,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask
PSF = nan(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths    
    % get PSF 
    PSF(:,:,ch) = propcustom_mft_PtoF(Epup(:,:,ch),foc, lambdas(ch), dx, dxi, Nxi, dxi, Nxi);
end

figure;
for ch = 1:numWavelengths    
    
    PSF_lam = PSF(:,:,ch);
    radPro = PSF_lam(Nxi/2+1,Nxi/2+1:end);
    semilogy(abs(radPro).^2);
    hold on;
    
end


% Calculate intensity (peak normalized)
%iPSF_tmp = abs(PSF).^2/normI;
% Take mean and store in image matrix

iPSF_BB = mean(abs(PSF).^2,3);
iPSF_BB = iPSF_BB/max(iPSF_BB(:));

figure(3)
imagesc(xvalsFP/lambda0Fnum_samp,yvalsFP/lambda0Fnum_samp,iPSF_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask
offsetX = 0;
offsetY = 0;

EPM = generateVortexMask( charge, coordsPP, [offsetX offsetY] );

central_band_index = ceil(numWavelengths/2);

figure(4)
imagesc(xvalsPP/apRad,yvalsPP/apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
axis image; 
axis([-1 1 -1 1]);
title('Pupil phase at \lambda_0');
colormap(hsv(256));
colorbar; 
drawnow;

%% Get PSF with vortex mask
PSFv = nan(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths    
    % get PSF 
    PSFv(:,:,ch) = propcustom_mft_PtoF(Epup(:,:,ch).*EPM(:,:,ch),foc, lambdas(ch), dx, dxi, Nxi, dxi, Nxi);
end

% Calculate intensity (peak normalized)
iPSFv_tmp = abs(PSFv).^2/normI;
% Take mean and store in image matrix
iPSFv_BB = mean(iPSFv_tmp,3);


figure(5)
imagesc(xvalsFP/lambda0Fnum_samp,yvalsFP/lambda0Fnum_samp,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Generate fibermode at each lambda     
%-- Iterate through wavelengths generating modes
fibmode = nan(Nxi, Nxi, numWavelengths);
for ch = 1:numWavelengths
    %fibmode(:,:,ch) = generateSMFmode_gaussian(1.4*(lambdas(ch)/lambda0)*lambda0Fnum_samp,coordsFP);
    %fibmode(:,:,ch) = generateSMFmode_gaussian(1.4*lambda0Fnum_samp,coordsFP);
    
    %-- Fiber modes not provided; generate them
    % Generate the fiber mode for this wavelength with proper scaling
%     fibmode(:,:,ch) = generateSMFmode( fiber_props, lambdas(ch), Fnum, lambdas(ch)/lambda0*lambda0Fnum_samp, coordsFP);

    %fiber_props.NA = sqrt(fiber_props.n_core^2-fiber_props.n_clad^2); 
    %fibmode(:,:,ch) = generateSMFmode_bessel( fiber_props.NA, fiber_props.core_rad, lambdas(ch), dxi, coordsFP );
    
    % Recale to sampling for image as set by falco
    %lam_frac = lambdas(ch)/lambda0;
    %fibmode(:,:,ch) = interp2(coordsfp.X,coordsfp.Y,fibmode(:,:,ch),coordsfp.X/lam_frac,coordsfp.Y/lam_frac,'linear',0);
    
    
    % Generate the fiber mode for this wavelength with proper scaling
	fibmode(:,:,ch) = generateSMFmode_mft( fiber_props, lambdas(ch), dxi, coordsFP);
    
end

%% Calculate Throughputs
%-- Get null depth
for ch = 1:numWavelengths
    eta_ss(ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end
eta_ss
%%
%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFv(:,:,ch), totalPower0, 7*lambda0Fnum_samp, coordsFP);

    % @GARY DO I NEED TO RESCALE SOMETHING LIKE THE NORMAL COUPLING MAP FUNC DOES?

end

%eta_maps = generateCouplingMap_polychromatic(PUPIL.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 5*lambdaOverD, coords, fibmode);

%% Display coupling maps
%-- Linear scale
figure(6);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP/lambda0Fnum_samp,yvalsFP/lambda0Fnum_samp,eta_maps(:,:,ch));
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale
figure(7);
for ch = 1:numWavelengths
    subplot(1,numWavelengths,ch);
    imagesc(xvalsFP/lambda0Fnum_samp,yvalsFP/lambda0Fnum_samp,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

disp('Key Coupling Points:')
for ch = 1:numWavelengths
    fprintf('lambda = %f nm,    on-axis null = %e,    max = %f %%\n',lambdas(ch)*1e9, eta_ss(ch), max(eta_maps(:,:,ch),[],'all')*100);
end