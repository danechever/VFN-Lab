%{
Script to determine the vortex centering in the pupil using pupil images 
taken of the VFN system with the vortex at the same position used in the 
actual Poly Paper scans.

This script also shows some other pupil images in this dataset.
%}
close all; clear all;
% add path to analysis library
addpath('../AnalysisLib')

%% Folder and file locations
% Folder containing the pupil images
fldnm = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/20210726_PupilImages_PolyPaper';
% Image file to use (asterisk wildcard allowed
  % NOTE: must be a 1x1 cell array 
nmWpath = {fullfile(fldnm,'4_*.fits')};


%% Load the images and average into 1
% Find all samples for the given set
STRNMS = VFN_An_getSTRNMS(nmWpath);
% Preallocate matrix 
PUPIL = nan([fitsinfo(STRNMS(1)).PrimaryData.Size, length(STRNMS)]);
% Iterate through samples loading each one
for ind = 1:length(STRNMS)
    PUPIL(:,:,ind) = fitsread(STRNMS(1));
end

%-- Average samples
PUPIL_O = mean(PUPIL,3);


%-- Display the original image
figure(); imagesc(PUPIL_O); axis image
title({'VFN Pupil Image - Original', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;

%% Crop to just barely fit the pupil

%-- Set areas outside pupil to 0 
%PUPIL(PUPIL < 13) = 0;
PUPIL = PUPIL_O > 28;

%-- Sharpen Edges
%windowSize = 20;
%fprintf('\n-- NOTE: using windowSize = %0.1f for edge sharpening on PUPIL\n',windowSize)
%kernel = ones(windowSize) / windowSize ^ 2;
%PUPIL = conv2(single(PUPIL), kernel, 'same');
%PUPIL = PUPIL > 0.5; % Rethreshold

%-- Crop image so it just barely fits the pupil
% Auto-detect the pupil boundaries
[~, row0] = find(PUPIL',1);
[~, col0] = find(PUPIL,1);
[~, row1] = find(PUPIL',1,'last');
[~, col1] = find(PUPIL,1,'last');
% Crop 
PUPIL = PUPIL(row0:row1,col0:col1);
PUPIL_O = PUPIL_O(row0:row1,col0:col1);

%% Display the thresholded and cropped pupil image
figure(); imagesc(PUPIL); axis image
drawnow;

%% Crop the image to the central region
  % We'll search for circles in this cropped image since there's less noise
  % towards the saturated center 
cropVal = 6;
%-- Crop matrix to central region in case the full donut is shown
rowmin = max(floor(size(PUPIL, 1)/cropVal),1);
colmin = max(floor(size(PUPIL, 2)/cropVal),1);
rowmax = (floor((cropVal-1)*size(PUPIL,1)/cropVal+1));
colmax = (floor((cropVal-1)*size(PUPIL,2)/cropVal+1));
PUPILCR = PUPIL(rowmin:rowmax,colmin:colmax,:,:,:,:);

%-- Show the cropped image
figure(); imagesc(PUPILCR); axis image
drawnow;

%% Search for circles in the cropped image
  % We'll assume the vortex is towards the center and that the hole it
  % causes is roughly circular

[cents, radii] = imfindcircles(PUPILCR, [30, 100], 'ObjectPolarity','dark','Sensitivity',0.95);


%-- Translate submatrix indexes to full matrix ones
cents(2)   = cents(2)+rowmin-1;
cents(1)   = cents(1)+colmin-1;

%% Display the final image with overlays

% Compute the pupil radius in each dimension
xpuprad = (size(PUPIL_O,2)/2);
ypuprad = (size(PUPIL_O,1)/2);

% Show the cropped but not thresholded full-pupil image
figure(); imagesc(PUPIL_O); axis image; axis xy
hold on

% Add the recognized circle
viscircles(cents, radii);

% Add vertical line marking center of pupil
xline(xpuprad+1, 'r', 'LineWidth', 2);

% Add horizontal line marking center of pupil
yline(ypuprad+1, 'r', 'LineWidth', 2);

% Add X-marker to the center of the vortex location
scatter(cents(1),cents(2),200,'rX', 'LineWidth', 4)

title({'VFN Pupil Image - Centered', 'Vortex in circle'})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')

drawnow;
%% Print the vortex offsets

%-- Compute offset in pixels
offsetX = cents(1) - (xpuprad +1);
offsetY = cents(2) - (ypuprad +1);

fprintf('\nVortex X-offset [pixels]: %f',offsetX)
fprintf('\nVortex Y-offset [pixels]: %f',offsetY)

%-- Convert to fractions of the pupil
offsetX_frac = offsetX / (2*xpuprad);
offsetY_frac = offsetY / (2*ypuprad);

fprintf('\nVortex X-offset [pup. frac.]: %f',offsetX_frac)
fprintf('\nVortex Y-offset [pup. frac.]: %f\n',offsetY_frac)

%% Display the final image with overlays (in fractions of pupil units)
% Set axes to be centered on the image (and hence on the pupil since the
% pupil is centered on the image after cropping). Also rescale axes to be
% to be in fractions of pupil size.
xax = ((1:size(PUPIL_O,2)) - xpuprad)/(2*xpuprad);
yax = ((1:size(PUPIL_O,1)) - ypuprad)/(2*ypuprad);

% Show the cropped but not thresholded full-pupil image
figure(); imagesc(xax,yax,PUPIL_O); axis image; axis xy
hold on

% Add the recognized circle
viscircles([offsetX_frac,offsetY_frac], radii/(2*xpuprad));

% Add vertical line marking center of pupil
xline(0, 'r', 'LineWidth', 2);

% Add horizontal line marking center of pupil
yline(0, 'r', 'LineWidth', 2);

% Add X-marker to the center of the vortex location
scatter(offsetX_frac,offsetY_frac,200,'rX', 'LineWidth', 4)

title({'VFN Pupil Image - Centered', 'Vortex in circle'})
xlabel('Spatial [Fraction of Pupil]')
ylabel('Spatial [Fraction of Pupil]')

drawnow;

%% Display an image of the pupil with the vortex out of the beam (Saturated)
% Image file to use (asterisk wildcard allowed)
  % NOTE: must be a 1x1 cell array 
nmWpathNoVort = {fullfile(fldnm,'7_*.fits')};

%----------------- Load the images and average into 1
% Find all samples for the given set
STRNMSNoVort = VFN_An_getSTRNMS(nmWpathNoVort);
% Preallocate matrix 
PUPILNoVort = nan([fitsinfo(STRNMSNoVort(1)).PrimaryData.Size, length(STRNMSNoVort)]);
% Iterate through samples loading each one
for ind = 1:length(STRNMSNoVort)
    PUPILNoVort(:,:,ind) = fitsread(STRNMSNoVort(1));
end

%-- Average samples
PUPIL_ONoVort = mean(PUPILNoVort,3);


%-- Display the original image
figure(); imagesc(PUPIL_ONoVort); axis image
title({'VFN Pupil Image - Vort Out - Original', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;

%---------------- Crop to just fit the pupil

%-- Set areas outside pupil to 0 
%PUPIL(PUPIL < 13) = 0;
PUPILNoVort = PUPIL_ONoVort > 28;

%-- Crop image so it just barely fits the pupil
% Auto-detect the pupil boundaries
[~, row0] = find(PUPILNoVort',1);
[~, col0] = find(PUPILNoVort,1);
[~, row1] = find(PUPILNoVort',1,'last');
[~, col1] = find(PUPILNoVort,1,'last');
% Crop 
PUPIL_ONoVort = PUPIL_ONoVort(row0:row1,col0:col1);

%----------------- Display
%-- Display the original image
figure(); imagesc(PUPIL_ONoVort); axis image
title({'VFN Pupil Image - Vort Out - Cropped', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;

%% Display an image of the pupil with the vortex out of the beam (Not Sat)
% Image file to use (asterisk wildcard allowed)
  % NOTE: must be a 1x1 cell array 
nmWpathNoVortNoSat = {fullfile(fldnm,'9_*.fits')};

%----------------- Load the images and average into 1
% Find all samples for the given set
STRNMSNoVortNoSat = VFN_An_getSTRNMS(nmWpathNoVortNoSat);
% Preallocate matrix 
PUPILNoVortNoSat = nan([fitsinfo(STRNMSNoVortNoSat(1)).PrimaryData.Size, length(STRNMSNoVortNoSat)]);
% Iterate through samples loading each one
for ind = 1:length(STRNMSNoVortNoSat)
    PUPILNoVortNoSat(:,:,ind) = fitsread(STRNMSNoVortNoSat(1));
end

%-- Average samples
PUPIL_ONoVortNoSat = mean(PUPILNoVortNoSat,3);


%-- Display the original image
figure(); imagesc(PUPIL_ONoVortNoSat); axis image
title({'VFN Pupil Image - Vort Out - Not Saturated - Original', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;

%---------------- Crop to just fit the pupil

%-- Set areas outside pupil to 0 
%PUPIL(PUPIL < 13) = 0;
PUPILNoVortNoSat = PUPIL_ONoVortNoSat > 5;

%-- Crop image so it just barely fits the pupil
% Auto-detect the pupil boundaries
[~, row0] = find(PUPILNoVortNoSat',1);
[~, col0] = find(PUPILNoVortNoSat,1);
[~, row1] = find(PUPILNoVortNoSat',1,'last');
[~, col1] = find(PUPILNoVortNoSat,1,'last');
% Crop 
PUPIL_ONoVortNoSat = PUPIL_ONoVortNoSat(row0:row1,col0:col1);

%----------------- Display
%-- Display the original image
figure(); imagesc(PUPIL_ONoVortNoSat); axis image
title({'VFN Pupil Image - Vort Out - Not Saturated - Cropped', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;
