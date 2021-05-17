%{
Script used to do a Zernike analysis of DM Data. 
** This script acts as an example for doing Zernike Analyses of Zygo Data

This script is a simplified version of previous Zernike Analysis scripts I 
have written before. The goal is to make it more readable and accessible to
users by removing extraneous features.

* Based on GenZernAnalysis_IterableAndSelectable_KPIC35mmOAP.m

This example script loads a fits file, applies a circular pupil per user 
inputs, subtracts unwanted zernikes, and then displays the resulting surface
and zernike coefficients. Other values are also reported and the script 
provies the option to save the resulting figures as .png and values as an
excel (.xlsx) file. The analysis is done over a user-provided "main" 
region of interest (ROI) and alternate ROI - this allows you to look at 
the WFE in a subsection of the optic as well as the full clear aperture.


NOTE: Zygo frames have different P/T/T values between individual samples
due to the way the measurements are taken. As such, any averaging or
subtraction of individual frames needs to happen after P/T/T correction.
This code does that. 

******************************************************
- Dependencies:
    zernfun2_DE     Function to calculate Zernike modes
******************************************************

Initial Version:    Daniel Echeverri, 10/20/2017
Last Edit:          Daniel Echeverri
Last Modified:      05/16/2021
%}


%% Define scaling constant
close all; clear all
m2nm    = 1e9;      % meter to nanometer conversion (zygo gives results in m)

%% Parameters to set
%-- Naming Parameters
MainSize = 13.6; % [mm] Size of region defined by "radius" variable below
AltSize= 15;    % [mm] Alternate region to analyze. (same center just larger or smaller radius)
% Filenames for saved results
fignmMain   = 'Main_Pup';   % main surface plot
fignmAlt    = 'Alt_Pup';    % alt. surface plot
fignmZern   = 'ZernCoeff';  % Zernike decomp bar plot
xlnm        = 'ZernCoeff';  % Excel file of zernike decomp
isSave      = true;      % Flag to choose whether results should be saved

%-- File/Analysis Parameters
% File Properties
fldnm   = '/home/nfiudev/dev/dechever/DMCharac/Data/UniformVoltage';  %Folder containing data
flnm    = 'UniformTest.fits';       %Filename of fits file to read
    % NOTE: file dims should be (imgX, imgY, samples, tes condition) 
ind2an  = 40;  % Index (along 4th dim) of frames to analyse (test condition)
    % NOTE: all frames along 3rd dim will be averaged (ie. those should all
    % be samples of the same test condition)
% Pupil Properties
centCoord = [120, 132]; % Center of pupil 
radius    = 56;
% Zernike Analysis
ApproxOrd = 36;      %Number of zernike modes used in decomposition/approximation
modRemove = [1,2,3];  %Vector representing modes to remove. DOES NOT FOLLOW REGULAR ORDER
        %First few modes: 1 = piston, 2 = tip, 3 = tilt, 4 = Astig Y, 5 = Focus, 6 = Astig 45

%% Load the fits file
fprintf('Loading file...\n')
fits    = fitsread(fullfile(fldnm,flnm)); 
fprintf('-- extracting ind %d\n',ind2an)
% Extract requested frames 
fits    = fits(:,:,:,ind2an);  % test condition

%Display 1 sample from imported data
samp = fits(:,:,1);
figure(1); imagesc(samp); axis image; title('Sample Test Frame')

%% Calculate Pupil
fprintf('Getting Region of Interest...\n')
altFrac = AltSize/MainSize;       %fraction of main aperture radius for alt pupil

xc = centCoord(1); yc = centCoord(2);
fprintf('-- xC = %d, yC = %d, rad = %d\n',xc, yc, radius)
figure(1);
viscircles([xc, yc], radius);  % Main aperature
viscircles([xc, yc], radius*altFrac,'LineStyle','--'); % Alt aperature
drawnow

%Define crop region
%Cropped region will be row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
rcrop = yc-radius;        %Row to start crop
ccrop = xc-radius;        %Column to start crop
croplength = radius*2;   %Size of crop

%-- ALTERNATE METHOD (commented below):
  % (automatically finds a circular pupil based on user click)
%{
% Params to change
edge    = 8;        %amount of edge to remove. Crops the "full" image by shrinking the radius to remove some edge effects 
rSearch = 25;       %Starting radius leeway. Radius search will be calculated value +- rSearch
                        %^^ Very useful if code fails with "Index exceeds matrix dimensions." becuase of "rcrop = centers(2)-radi"
                        %25 is usually good but 50 can be better for larger pupils
%prompt user to click on center of mirror
fprintf('Please click near the center of the circle to analyze\n');
[x, y]  = ginput(1);
fprintf('Calculating Region...\n');
rows    = samp(:,round(x));         %extract elements of center column
row0    = find(~rows);              %find elements that are == zero
radi    = find(row0 <= y, 1, 'last');   %get bottom edge of mirror
radi    = round(y - radi);          %get aproximate radius

%find mirror center and radius
%resize subtraction term if needed to prevent negative radius searches
tmp = rSearch;
if radi <= rSearch
    tmp = radi-1;
end
[centers, radi] = imfindcircles(samp(:,:,1)~=0,[radi-tmp, radi+rSearch],'Sensitivity',0.91);

%If multiple circles were found, use circle with center closest to user click
if numel(radi)>1
    mindist = 1000;
    for n = 1:numel(radi)
        tmp = [centers(n,:); x,y];
        tmp = pdist(tmp);
        if tmp<mindist
            mindist = tmp;
            ind     = n;
        end
    end
    centers = centers(ind,:);
    radi    = radi(ind);
end

centers = round(centers); radi = round(radi-edge);   %round results for appropriate indexing
viscircles(centers, radi);          %plot calculated circle
viscircles(centers, radi*altFrac,'LineStyle','--'); % Alt aperature
drawnow

%Cropped region will be row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
%MANUALLY OVERWRITE CROP PARAMETERS HERE TO ENTER YOUR OWN VALUES
rcrop = centers(2)-radi;        %Row to start crop
ccrop = centers(1)-radi;        %Column to start crop
croplength = radi*2;   %Size of crop
%}

%% Crop and prepare data
fprintf('Cropping data...\n')
%Calculate alt-pupil crop parameters
AltCropL = floor(croplength*altFrac);               %width of alt aperture
Altrcrop = rcrop+(floor((croplength-AltCropL)/2));  %row crop start for alt ap
Altccrop = ccrop+(floor((croplength-AltCropL)/2));  %col crop start for alt ap

% Crop 
Altfits = fits(Altrcrop:Altrcrop+AltCropL, Altccrop:Altccrop+AltCropL,:);
fits = fits(rcrop:rcrop+croplength, ccrop:ccrop+croplength,:); 

% Zygo sets points it couldn't read to 0. Thus, assume all 0 values aren't 
% real data and set to NaN to remove from analysis.
    %odds of a real point being at exactly 0 are really slim...
fits(~fits) = NaN;
Altfits(~Altfits) = NaN;
                                
%Scale appropriately (convert all phase data to nm)
fits = fits*m2nm;
Altfits = Altfits*m2nm;

%% Perform Zernike Correction
fprintf('Creating Unit Zernikes...\n')
%Calculate the unit zernikes over both main and alt pupils
[zernM, inCircM] = unitZern(fits, ApproxOrd);      
[zernA, inCircA] = unitZern(Altfits, ApproxOrd);   

fprintf('Decomposing and Subtracting Zernikes...\n');
nSamp = size(fits,3);
%loop through Zer. Decomp. and mode subtraction               
for n = 1:nSamp
    fits(:,:,n) = zernSubt(fits(:,:,n), zernM, inCircM, modRemove);         %Corrected main image surface
    Altfits(:,:,n)  = zernSubt(Altfits(:,:,n), zernA, inCircA, modRemove);    %Corrected alt aperture surface
end

%% Final processing (averaging across samples)
fprintf('Analyzing the Results...\n')
%Calculate the mean of the images
  %"omitnan" is useful for averaging out holes within maps.
meanSurf = mean(fits,3,'omitnan');  
meanSurfA = mean(Altfits,3,'omitnan');

%% Display Results

%Calculate max and min for a uniform colorbar. 
  % One of the surfs  might saturate with this colorbar but this
  % lets us plot both on the same colorbar for comparison
cmin = min(min(meanSurf));
cmax = max(max(meanSurf));

% Plot the main Surface
figM = figure();
imagesc(meanSurf,[cmin cmax]); colorbar; axis image
title(['Main Pupil (' num2str(MainSize, '%4.2f') 'mm)'])

%Plot alt aperture
figA = figure();
imagesc(meanSurfA, [cmin cmax]); colorbar; axis image
title(['Alt Pupil (' num2str(AltSize, '%4.2f') 'mm)'])

%Calculate rms error over full and sub pupils
vals = ~isnan(meanSurf);    %remove nan values from calculation
altvals = ~isnan(meanSurfA);
rmsMain = sqrt(mean(meanSurf(vals).^2));    %full rms
rmsAltAp = sqrt(mean(meanSurfA(altvals).^2));    %subAperture rms
fprintf('\nResults:\nRMS Main Pupil: %6.4f',rmsMain); 
fprintf('\nRMS Alt Pupil : %6.4f',rmsAltAp); 

%Calculate remaining Zernike terms over the avg. and corrected main aperture
[~, cof] = zernSubt(meanSurf, zernM, inCircM, modRemove);
%Calculate remaining Zernike terms over the avg. and corrected alt aperture
[~, cofA] = zernSubt(meanSurfA, zernA, inCircA, modRemove);
%Use coefficients to calculate the RMS error as provided by zernike approximation
    %Sum in quadrature
rmsCof = sqrt(sum(cof.^2));
fprintf('\nRMS Reconstructed: %6.4f',rmsCof); 

%Calculate RMS coma totals
com1Cof = sqrt(sum(cof(8:9).^2));
com2Cof = sqrt(sum(cof(18:19).^2));
fprintf('\nRMS 1st Order Coma Total: %6.4f',com1Cof); 
fprintf('\nRMS 2nd Order Coma Total: %6.4f \n',com2Cof); 

%Uncomment section below to display remaining Zernike Coefficients
%Plot bar graph of first 21 Zernikes
figZ = figure();
bar(cof(1:21))
title('First 21 Zernike Coefficients Over the Main Pupil')
xlabel({'Coefficient','*does not follow regular order'})
ylabel('RMS Contribution (nm)')
xlim([0 22])
grid on

%% Create a table into which results can be saved

%Create table of first 15 Zernikes + 2nd order coma
cofThres = .1;
cofT    = [cof(1:15); cof(18:19)];
cofi    = abs(cofT)<cofThres;
cofT    = cellstr(num2str(cofT, '%+06.2f'));
cofT(cofi)    = {'------'};
%Provide formatting for table
cofT1   = {'N index' '0' '1' '1' '2' '2' '2' '3' '3' '3' '3' '4' '4' '4' '4' '4' '5' '5'};
cofT2   = {'M index' '0' '1' '-1' '2' '0' '-2' '3' '1' '-1' '-3' '4' '2' '0' '-2' '-4' '1' '-1'};
cofT3   = {'Name' 'Piston' 'Tip' 'Tilt' 'Astig - Y' 'Focus' 'Astig - 45' 'Trefoil - 45' ...
    'Coma - X' 'Coma - Y' 'Trefoil - Y' 'QdFoil - Y' '2nd Astig - Y' 'Spherical' ...
    '2nd Astig - 45' 'QdFoil - 45' '2nd Coma - X' '2nd Coma -Y'};
cofTab = [cofT1' cofT2' cofT3' ['RMS Contribution (nm)'; cofT]];

%Provide RMS Values
anaT    = cell(13,4);
anaT(1,1:2)   = ['RMS Flatness: Calculated Main (nm)' cellstr(num2str(rmsMain))];
anaT(2,1:2)   = ['RMS Flatness: Calculated Alt. (nm)' cellstr(num2str(rmsAltAp))];
anaT(3,1:2)   = ['RMS Flatness: Zern Approx. (nm)' cellstr(num2str(rmsCof))];
anaT(4,1:2)   = ['RMS Flatness: 1st Coma Total (nm)' cellstr(num2str(com1Cof))];
anaT(5,1:2)   = ['RMS Flatness: 2nd Coma Total (nm)' cellstr(num2str(com2Cof))];

%Provide Analysis Parameters
anaT(6,1:2)   = ['Sample File' cellstr(flnm)];
anaT(7,1:2)   = ['# Samples' cellstr(num2str(nSamp))];
anaT(8,1:2)   = ['Crop Regions:' {'_____'}];
anaT(9,1:2)   = ['   rcrop' cellstr(num2str(rcrop))];
anaT(10,1:2)   = ['   ccrop' cellstr(num2str(ccrop))];
anaT(11,1:2)   = ['   croplength' cellstr(num2str(croplength))];
anaT(12,1:2)  = ['Alt Aperture' cellstr(num2str(altFrac))];
anaT(13,1:4)    = [{'Zernike Coefficients'} {'_____'} {'_____'} {'_____'}];

%% Save results

if isSave
    saveas(figM, fullfile(fldnm, [fignmMain '.png']));
    saveas(figA, fullfile(fldnm, [fignmAlt '.png']));
    saveas(figZ, fullfile(fldnm, [fignmZern '.png']));
    writecell([anaT; cofTab], fullfile(fldnm, [xlnm '.xlsx']), 'UseExcel', false)
end

%% Helper Functions
%function to calculate the unit Zernikes over the given area (rawSurf)
function [zerns, inCirc] = unitZern(rawSurf, approxOrder)
    %Create grid coordinate matrices in polar of size of rawSurf
    L = size(rawSurf,1);            
    X = -1:2/(L-1):1;               %center origin of coord. system
    [x,y] = meshgrid(X);            %create cartesian grid
    x = x(:); y = y(:);             %variable for polar translation: cartesian coordinates in plane
    [theta, r] = cart2pol(x,y);     %convert to polar coordinates

    %Compute Zernike decomp of rawSurf on unit circle
    inCirc = (r<=1);    %boolean matrix to indicate unit circle area for analysis
    zerns  = zernfun2_DE(0:approxOrder, r(inCirc), theta(inCirc),'norm');  %Calculate Zernikes approximating to approxOrder (these are the standard, unscaled Zernike modes)
end

%function that decomposes and subtracts out modes identified by "subModes"
        %subModes should be a vector of integer corresponding to
            %coefficients IN ORDER PROVIDED BY ZERNFUN2_DE() FUNCTION
function [correctedSurf, coeffs] = zernSubt(rawSurf, zerns, inCirc, subModes)
    %Remove NaN values from the analysis region:
    circ    = rawSurf(inCirc);  %vector of numbers within the pupil
    circI   = ~isnan(circ);     %vector of indexes of non-nan (ie. number) values within the circular pupil
                                    %Allows extraction of number values from the region s.t. analysis ignores NaN values
    numSurf = circ(circI);      %Vector of number values within circular pupil (NaN values now excluded)
    numzern = zerns(circI,:);     %Vector of values on zernike map corresponding to non-NaN values of surface
    
    %Calculate coefficients for Zernikes by dividing rawSurf with the unit Zernikes provided
    coeffs = numzern\numSurf;   %Ignores NaN values in RawSurf
    
    %Sum the Zernike orders that will be subtracted from the original image
    subt = NaN(size(rawSurf));
    subt(inCirc) = zerns(:, subModes)*coeffs(subModes);     %Calculate by multiplying the coeffs with the unit Zernikes
    
    %Subtract the undesired modes from the original image
        %This subtraction implicitly removes values outside the analysis
        %circle since: number - NaN = NaN
    correctedSurf = rawSurf - subt;
end
