%{
%}

clear; close all; 
addpath(genpath(fullfile('..','..','..','VFN-Simulations/VFNlib/')));
addpath(genpath(fullfile('..','..','..','falco-matlab')))

tic
%% Input parameters 

%-- Provide regular parameters
% Define smapling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % Aperture radius in samples 

%-- Define wavelength info
lambda0 = 650e-9; %central wavelength
fracBW = 1e-15;%100/650;% \Delta\lambda/\lambda
numWavelengths = 1;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%-- Parameters for Thorlabs SM600
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 4.3e-6/2;% Core radius [um]
fiber_props.n_core = 1.4601;% core index (interp @650nm)
fiber_props.n_clad = 1.4567;% cladding index (interp @650nm)
fiber_props.type = 'bessel';
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.42); %ideal for bess: 1.238 % focal ratio of the beam at the fiber

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
lambdaOverD_samp = 26; %[samples/ (lam/D)] in units of samples in the image plane
im_size = 15;    %[lam0/D] field of view in final image plane

%-- Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

%-- Define wavefront error at the central wavelength
  % RMS WFE per noll index
nolls = 3;
coeffs = 0.0;

%-- Give offsets for the vortex mask
offsetX = 0*0.0118*apRad*2;  
offsetY = 0*0.0256*apRad*2;

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvalsPP = coordsPP.xvals;% Helpful for plotting
yvalsPP = coordsPP.yvals;

%-- Coordinates in the pupil plane
Nxi = im_size*lambdaOverD_samp;
coordsFP = generateCoordinates( Nxi );
xvalsFP = coordsFP.xvals;% Helpful for plotting
yvalsFP = coordsFP.yvals;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in radians at final focal plane 
lambdaOverD_rad = lambda0*foc/DPup;     
% "pixel" size in pupil and image planes
dx  = DPup/(2*apRad);
dxi = lambdaOverD_rad/lambdaOverD_samp;


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

% Calculate intensity (peak normalized)
iPSF_tmp = abs(PSF).^2/normI;
% Take mean and store in image matrix
iPSF_BB = mean(iPSF_tmp,3);


figure(3)
imagesc(xvalsFP/lambdaOverD_samp,yvalsFP/lambdaOverD_samp,iPSF_BB);
axis image; 
axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Make vortex mask
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
imagesc(xvalsFP/lambdaOverD_samp,yvalsFP/lambdaOverD_samp,iPSFv_BB);
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
    %fibmode(:,:,ch) = generateSMFmode_gaussian(1.4*(lambdas(ch)/lambda0)*lambdaOverD_samp,coordsFP);
    
    %-- Fiber modes not provided; generate them
    % Generate the fiber mode for this wavelength with proper scaling
    fibmode(:,:,ch) = generateSMFmode( fiber_props, lambdas(ch), Fnum, lambdaOverD_samp, coordsFP);
    
    % Recale to sampling for image as set by falco
    %lam_frac = lambdas(ch)/lambda0;
    %fibmode(:,:,ch) = interp2(coordsfp.X,coordsfp.Y,fibmode(:,:,ch),coordsfp.X/lam_frac,coordsfp.Y/lam_frac,'linear',0);
end

%% Display Overlap Integral Visualizer
%-- Get null depth
for ch = 1:numWavelengths
    eta_ss(ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,numWavelengths);
for ch = 1:numWavelengths
    
%     crpW = 5;  %[lambda/D] FOV (radius wr.t center of image)
%     xcrop = xvalsFP(abs(xvalsFP/lambdaOverD_samp) < crpW)+coordsFP.N/2+1;
%     ycrop = yvalsFP(abs(yvalsFP/lambdaOverD_samp) < crpW)+coordsFP.N/2+1;
%     xax = xvalsFP(xcrop)/lambdaOverD_samp; yax = yvalsFP(ycrop)/lambdaOverD_samp;
%     PSFvNRM = PSFv/sqrt(totalPower0);
%     PSFvNRMCRP = PSFvNRM(xcrop,ycrop);
%     fibmodeCRP = fibmode(xcrop,ycrop);
%     figure(); surf(xax,yax,abs(fibmodeCRP),angle(fibmodeCRP),'FaceAlpha', 0.9,'EdgeColor','black','EdgeAlpha', 0.5); hold on;
%     surf(xax,yax,abs(PSFvNRMCRP),angle(PSFvNRMCRP),'FaceAlpha', 0.5,'EdgeColor','none');
%     colormap('hsv')
%     OLap = fibmodeCRP.*PSFvNRMCRP;
%     OLap(abs(OLap) < 1e-8) = 0;
%     figure(); surf(xax,yax,abs(OLap), angle(OLap),'EdgeColor','none'); colormap('hsv');
%     figure(); imagesc(xax,yax,angle(OLap)); colormap('hsv'); axis square; axis xy;
    

    %% Attempt 1
    % Define crop window for most displays
    crpW = 5;  %[lambda/D] FOV (radius w.r.t center of image)
    xax = xvalsFP/lambdaOverD_samp; yax = yvalsFP/lambdaOverD_samp;
    
    % Define crop window for overlayed display
    crpWO = 5;  %[lambda/D] 
    
    % Define crop window (circular) for functions themselves
    crpWF = 3*lambdaOverD_samp;  %[lambda/D]
    
    % Normalize the PSF (norm by total power)
    PSFvNRM = PSFv/sqrt(totalPower0);
    
    % Crop the functions to the inner 
    PSFvNRMCRP = PSFvNRM;
    %PSFvNRMCRP(coordsFP.RHO > crpWF) = nan;
    fibmodeCRP = fibmode;
    %fibmodeCRP(coordsFP.RHO > crpWF) = nan;
        
    % Display the PSF and Fibermode overlayed
    figure('Color', 'white', 'units', 'normalized', 'Position', [0 0 1, 1]); 
    surf(xax,yax,abs(fibmodeCRP),angle(fibmodeCRP),'FaceAlpha', 0.9,'EdgeColor','black','EdgeAlpha', 0.2); hold on;
    surf(xax,yax,abs(PSFvNRMCRP),angle(PSFvNRMCRP),'FaceAlpha', 0.7,'EdgeColor','none');
    colormap('hsv')
    axis([-crpWO, crpWO,-crpWO, crpWO]);
    grid minor;
    set(gca,'FontSize', 24)
    % Get zlimits so we can re-use them in the next plots
    zax = zlim;
    drawnow
    
    % Display just the fibmode
    figure('Color', 'white', 'units', 'normalized', 'Position', [0 0 1, 1]); 
    surf(xax,yax,abs(fibmodeCRP),angle(fibmodeCRP),'FaceAlpha', 0.9,'EdgeColor','black','EdgeAlpha', 0.2);
    colormap('hsv')
    axis([-crpW, crpW,-crpW, crpW, zax]);
    grid minor;
    set(gca,'FontSize', 24)
    drawnow
    
    % Display just the PSF
    figure('Color', 'white', 'units', 'normalized', 'Position', [0 0 1, 1]); 
    surf(xax,yax,abs(PSFvNRMCRP),angle(PSFvNRMCRP),'FaceAlpha', 0.5,'EdgeColor','none');
    colormap('hsv')
    axis([-crpW, crpW,-crpW, crpW, zax]);
    grid minor;
    set(gca,'FontSize', 24)
    drawnow
    
    % Compute the argument for the overlap integral
    OLap = fibmodeCRP.*PSFvNRMCRP;
    % Set argument to 0 wherever it is negligably small
    OLap(abs(OLap) < 1e-7) = 0;
    
    % Display the overlap integral in 3D
    figure('Color', 'white', 'units', 'normalized', 'Position', [0 0 1, 1]); 
    surf(xax,yax,abs(OLap), angle(OLap),'EdgeColor','none'); 
    colormap('hsv');
    axis([-crpW, crpW,-crpW, crpW]);
    grid minor;
    set(gca,'FontSize', 24)
    drawnow
    
    %  Display the overlap integral in 2D
    figure('Color', 'white', 'units', 'normalized', 'Position', [0 0 1, 1]); 
    imagesc(xax,yax,angle(OLap),'AlphaData',~isnan(OLap)); 
    colormap('hsv'); axis square; axis xy;
    axis([-crpW, crpW,-crpW, crpW]);
    grid on;
    set(gca,'FontSize', 24)
    drawnow

    %% Attempt 2
%     % Define crop window for most displays
%     crpW = 5;  %[lambda/D] FOV (radius w.r.t center of image)
%     xax = xvalsFP/lambdaOverD_samp; yax = yvalsFP/lambdaOverD_samp;
%     
%     % Define crop window for overlayed display
%     crpWO = 5;  %[lambda/D] 
%     
%     % Normalize the PSF (norm by total power)
%     PSFvNRM = PSFv/sqrt(totalPower0);
%         
%     % Display the PSF and Fibermode overlayed
%     figure(); 
%     surf(xax,yax,-abs(fibmode),angle(fibmode),'FaceAlpha', 0.9,'EdgeColor','black','EdgeAlpha', 0.5); hold on;
%     surf(xax,yax,abs(PSFvNRM),angle(PSFvNRM),'FaceAlpha', 0.5,'EdgeColor','none');
%     colormap('hsv')
%     axis([-crpWO, crpWO,-crpWO, crpWO]);
%     grid minor;
%     % Get zlimits so we can re-use them in the next plots
%     zax = zlim;
%     drawnow
%     
%     % Display just the fibmode
%     figure(); 
%     surf(xax,yax,-abs(fibmode),angle(fibmode),'FaceAlpha', 0.9,'EdgeColor','black','EdgeAlpha', 0.5);
%     colormap('hsv')
%     axis([-crpW, crpW,-crpW, crpW, zax]);
%     grid minor;
%     drawnow
%     
%     % Display just the PSF
%     figure(); 
%     surf(xax,yax,abs(PSFvNRM),angle(PSFvNRM),'FaceAlpha', 0.5,'EdgeColor','none');
%     colormap('hsv')
%     axis([-crpW, crpW,-crpW, crpW, zax]);
%     grid minor;
%     drawnow
%     
%     % Compute the argument for the overlap integral
%     OLap = fibmode.*PSFvNRM;
%     % Set argument to 0 wherever it is negligably small
%     OLap(abs(OLap) < 1e-7) = 0;
%     
%     % Display the overlap integral in 3D
%     figure(); surf(xax,yax,abs(OLap), angle(OLap),'EdgeColor','none'); 
%     colormap('hsv');
%     axis([-crpW, crpW,-crpW, crpW]);
%     grid minor;
%     drawnow
%     
%     %  Display the overlap integral in 2D
%     figure(); imagesc(xax,yax,angle(OLap)); 
%     colormap('hsv'); axis square; axis xy;
%     axis([-crpW, crpW,-crpW, crpW]);
%     grid on;
%     drawnow
%     
    
    %% Misc.
    % (side note: display PSF similar to how Gary did)
      % Note, need to set values outside rho to 1e-30, not nan
    %figure(); imagesc(angle(PSFvNRMCRP),'AlphaData',abs(PSFvNRMCRP)/max(abs(PSFvNRMCRP(:)))); colormap(hsv(256));set(gca,'Color','k')
    
    %% Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFv(:,:,ch), totalPower0, 7*lambdaOverD_samp, coordsFP);
    
    figure('Color', 'White');
    imagesc(xvalsFP/lambdaOverD_samp, yvalsFP/lambdaOverD_samp, eta_maps(:,:,ch));
    colormap(gray); axis image;
    axis([-2.5, 2.5, -2.5, 2.5]);
    xlabel('Angular Separation [\lambda/D]')
    ylabel('Angular Separation [\lambda/D]')
    title('2D Coupling Map')
    set(gca,'TickDir','out')
    set(gca,'FontSize', 24)
end

%eta_maps = generateCouplingMap_polychromatic(PUPIL.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 5*lambdaOverD, coords, fibmode);
toc