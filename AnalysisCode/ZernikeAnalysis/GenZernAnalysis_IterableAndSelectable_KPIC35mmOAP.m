%{
**THIS SECTION NOT UPDATED**

Script used to do a Zernike analysis of any data. Provides three pupil
modes: autocenter and autosize, autocenter and user-defined size,
user-defined center and user-defined size. Can remove any zernike modes
user defines. Performs analysis of full and user-defined sub-pupil.
Computes the RMS error over both pupils. Creates a plot of both the full 
and the sub-pupil. 

*** Based on OAP[...]_4 Code:
    *Prints RMS values in the comand window.
    - Controllable input parameters have all been brought to the top:___
        s.t. all user modifications between data sets occur in one seciton.
    - Code provides status and feedback in the command window in case of
        long runtimes to give user input on current process.

*** Change for From OAP2 to this version
    Prints a third RMS value, which is the RMS value as calculated by
    summing the zernike coefficients in quadrature. This allows comparison
    between this value and the rmsSubAp which is calculated by taking the
    rms directly from the points within the unit circle of the image.

*** Change from version _2_ to _2_3
    - Normalizes the Zernikes over their individual RMS values to return the RMS
    contribution of each.
    - Optimization of code by only calculating the Zernikes once.
        Previous code recalculated the unit Zernikes every time the zernike
        subtraction function was called. This was unnecessary since the
        same fundamental modes exist for each one, all that needs to be
        done is the division to determine the coefficients for each sample.
        As such, this code has a function to calculate the unit zernikes
        and then these are used repeatedly from then on.

*** MAJOR DIFFERENCE FROM PREVIOUS CODES (versions _2 and beyond):
    - This version does the sub pupil analysis completely separately
    from the full analysis. The old codes simply cropped the results of the
    full image down to the sub aperture. However, this means any local
    piston/T/T (primarily piston) will still be in the value. Thise will
    not actually be seen by the beam so the zernike decomposition should be
    performed locally on the sub aperture to remove P/T/T from this region
    alone. 

%Crops image to specific size and location

%Allows for circular mask to be applied for subaperture analysis. Circle
origin will be at center of full pupil.

%NEED TO CONFIRM: ***[CONFIRMED - v. _4]*** 
    the vector returned by zernfun2() is assumed to only include values 
    within the unit circle. That is, it is the values of the zernike modes 
    within r<=1, in some order. This was determined by looking at the code 
    and seeing that it only deals with inCirc whenever working with the 
    columns of the zernike matrix. 
******************************************************
- Arguments:
    fldnm           Folder containing the data 
    fnm             Name of files to analyze.
    isPupil         Defines type of pupil to use:
                        0 = Uses a Pupil that is centered and inscribed on the provided image
                        1 = Uses a Pupil that is centered on the image but has a radius defined by the user
                        2 = Uses a user-defined Pupil where user provides center and radius 
    subAper         Fraction of full aperture radius over which sub analysis 
                        is performed on. 
    ApproxOrd       Number of zernike modes used in decomposition. 
                        Actual decomposition will be ApproxOrd+1 modes.
                        36 recommended, higher orders increase runtime.
    modRemove       Vector representing modes to remove. DOES NOT FOLLOW REGULAR ORDER
                        First few modes: 1 = piston, 2 = tip, 3 = tilt, 4 = Astig Y, 5 = Focus, 6 = Astig 45

- Returns:
    NONE            Plots a figure showing a single full image and the
                        corresponding circle/crop region
                    Plots the resulting WF over the full and sub pupils

            Readily available variables:
                - RMS over full and Sub pupil, RMS over sub pupil calculated from the zernike approximation and coefficients (rms...)
                - Raw and subtracted surface maps (...Surf...)
                - Unit zernike map over calculated region (zern...)
                - Zernike Coefficients of sub pupil (cof). Easy enough to modify to obtain other coefficients.
                - Unit circle region indexes (inCirc...)
- Dependencies:
    zernfun2_DE     Function to calculate Zernike modes
                        Can be downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/7687-zernike-polynomials
******************************************************

Initial Version:    Daniel Echeverri, 10/20/2017
Last Edit:          Daniel Echeverri
Last Modified:      03/12/2018
%}

%% Parameters to set
%Naming Parameters
subSize = 6;    %Size of subaperture (in mm)
fullSize= 23;    %Size of full optic (in mm)
fignmFull   = 'Full_Pup';
fignmSub    = 'Sub_Pup';
fignmZern   = 'ZernCoeff';
xlnm        = 'ZernCoeff';

%File/Analysis Parameters
isItr   = 1;        %Flag to mark if code should iterate
                        %iterates over all fits files in given directory and
                        %averages them into 1 [ie. NOT FOR OAP ANALYSIS]
isPupil = 2;        %Defines type of pupil to use:
                        %0 = Uses a Pupil that is centered and inscribed on the provided image
                        %1 = Uses a Pupil that is centered on the image but has a radius defined by the user
                        %2 = Uses a user-defined Pupil where user provides center and radius 
                        %3 = Automatically find Pupil based on user clicking near center of circle
%IF USING isPupil = 3:
edge    = 8;        %amount of edge to remove. Crops the "full" image by shrinking the radius to remove some edge effects 
rSearch = 25;       %Starting radius leeway. Radius search will be calculated value +- rSearch
                        %^^ Very useful if code fails with "Index exceeds matrix dimensions." becuase of "rcrop = centers(2)-radi"
                        %25 is usually good but 50 can be better for larger pupils
                        
fldnm   = 'C:\Users\Daniel\Documents\Mawet201718\VortexFiberNuller_VFN\PoRT_Design\KPIC_35mmOAP_Analysis\Flatness_ball4';  %Folder containing data
fnm     = 'Flatness_Reflection_RForward_001.fits';   %Filename if not meant to iterate over all files
subAper = subSize/fullSize;       %fraction of full aperture radius over which to calculate rms
ApproxOrd = 36;      %Number of zernike modes used in decomposition/approximation
modRemove = [1,2,3];  %Vector representing modes to remove. DOES NOT FOLLOW REGULAR ORDER
        %First few modes: 1 = piston, 2 = tip, 3 = tilt, 4 = Astig Y, 5 = Focus, 6 = Astig 45
        
%% Begin Code_______________
addpath('C:\Users\Daniel\Documents\Mawet201718')
m2nm    = 1e9;      %micrometer to nanometer conversion 

%% Calculate Analysis Region
%Load all necessary filenames
if isItr
    fls = dir([fldnm '\*.fits']);
    fnms = strings(length(fls), 1);     %There's definitely a better way to deal with the strings vs. char here but oh well
    for i = 1:length(fls)
        fnms(i) = fls(i).name;
    end
    fnms = char(fnms);
else
    fnms = fnm;    %Filename for data when not iterating
end
nfiles = size(fnms,1);
%Import file. 'fits' name leftover from older versions, did not change to maintain readability
fits    = fitsread([fldnm filesep fnms(1,:)]);   %returns .textdata (header) and .data
%Display imported data
surf(fits)
view(0,90)
axis equal
shading interp

%Calculate Pupil
if isPupil == 0       %Use centered, inscribed pupil
    fprintf('Auto-centered Pupil Mode\n');
    [ldim, ind] = min(size(fits)); 
    if ind == 1
        rcrop = 1;
        ccrop = round((size(fits,2)-ldim)/2);
    else
        ccrop = 1;
        rcrop = round((size(fits,1)-ldim)/2);
    end
    croplength = ldim-1;
    viscircles([ccrop+(croplength-1)/2, rcrop+(croplength-1)/2], (croplength-1)/2, 'LineStyle','--');          %plot calculated circle
    drawnow
elseif isPupil == 1
    fprintf('Centered Pupil with user-defined Radius\n');
    xc = round(size(fits,2)/2); yc = round(size(fits,1)/2);
    fprintf('Please click at the edge of the analysis pupil\n')
    [xr, yr]  = ginput(1);
    xr = round(xr); yr = round(yr);
    %calculate radius as distance between both clicks
    radi    = round(pdist([xc,yc;xr,yr]));
    viscircles([xc, yc], radi);          %plot calculated circle
    drawnow

    %Check that user pupil does not clip edges and correct as needed
    if radi>=xc || radi>=yc || radi+yc>size(fits,1) || radi+xc>size(fits,2)
        radi = min([xc-1,yc-1, size(fits,1)-yc, size(fits,2)-xc]);
    end
    %Display final pupil in dotted red line
    viscircles([xc, yc], radi, 'LineStyle','--');          %plot calculated circle
    drawnow
    
    %Define resulting crop region
    %Cropped region will be row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
    rcrop = yc-radi;        %Row to start crop
    ccrop = xc-radi;        %Column to start crop
    croplength = radi*2;   %Size of crop
    
elseif isPupil ==2                %Use user-defined pupil
    %prompt user to define pupil center
    fprintf('User-defined Pupil Mode\n');
    fprintf('Please click on the center of the analysis pupil \n');
    [xc, yc]  = ginput(1);
    xc = round(xc); yc = round(yc);
    fprintf('Please click at the edge of the analysis pupil\n')
    [xr, yr]  = ginput(1);
    xr = round(xr); yr = round(yr);
    %calculate radius as distance between both clicks
    radi    = round(pdist([xc,yc;xr,yr]));
    viscircles([xc, yc], radi);          %plot calculated circle
    drawnow

    %Check that user pupil does not clip edges and correct as needed
    if radi>=xc || radi>=yc || radi+yc>size(fits,1) || radi+xc>size(fits,2)
        radi = min([xc-1,yc-1, size(fits,1)-yc, size(fits,2)-xc]);
    end
    %Display final pupil in dotted red line
    viscircles([xc, yc], radi, 'LineStyle','--');          %plot calculated circle
    drawnow
    
    %Define resulting crop region
    %Cropped region will be row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
    rcrop = yc-radi;        %Row to start crop
    ccrop = xc-radi;        %Column to start crop
    croplength = radi*2;   %Size of crop
elseif isPupil ==3
    %prompt user to click on center of mirror
    fprintf('Please click near the center of the circle to analyze\n');
    [x, y]  = ginput(1);
    fprintf('Calculating Region...\n');
    rows    = fits(:,round(x));         %extract elements of center column
    row0    = find(~rows);              %find elements that are == zero
    radi    = find(row0 <= y, 1, 'last');   %get bottom edge of mirror
    radi    = round(y - radi);          %get aproximate radius

    %find mirror center and radius
    %resize subtraction term if needed to prevent negative radius searches
    tmp = rSearch;
    if radi <= rSearch
        tmp = radi-1;
    end
    [centers, radi] = imfindcircles(fits~=0,[radi-tmp, radi+rSearch],'Sensitivity',0.91);

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
    drawnow

    %Cropped region will be row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
    %MANUALLY OVERWRITE CROP PARAMETERS HERE TO ENTER YOUR OWN VALUES
    rcrop = centers(2)-radi;        %Row to start crop
    ccrop = centers(1)-radi;        %Column to start crop
    croplength = radi*2;   %Size of crop
else
    error('Invalid Pupil Mode: please correct "isPupil" variable in "Parameters" section');
end
fprintf('Importing Files...\n')
%% Crop and prepare data
%preallocate space for imported files
files = zeros(croplength+1,croplength+1,nfiles);            %Full image files
%Calculate sub-pupil crop parameters
subCropL = floor(croplength*subAper);                              %width of sub aperture
subrcrop = rcrop+(floor((croplength-subCropL)/2));                 %row crop start for sub ap
subccrop = ccrop+(floor((croplength-subCropL)/2));                 %col crop start for sub ap
subApfiles = zeros(subCropL+1,subCropL+1,nfiles);       %Sub Aperture files
%loop to iterate through opening and cropping the files
for n = 1:nfiles
    fits = fitsread([fldnm filesep fnms(n,:)]);   %returns .textdata (header) and .data
    %Crop file and save it into the existing files stack
    files(:,:,n) = fits(rcrop:rcrop+croplength, ccrop:ccrop+croplength);    
    subApfiles(:,:,n) = fits(subrcrop:subrcrop+subCropL, subccrop:subccrop+subCropL);
end

files(~files) = NaN;        %Assume all 0 values are not real data and set to NaN to remove from analysis
subApfiles(~subApfiles) = NaN;
                                %odds of a point being at exactly 0 are really slim...
%Scale appropriately
files = files*m2nm ;        %convert all phase data to nm
subApfiles = subApfiles*m2nm;

%% Perform Zernike Analysis
fprintf('Creating Unit Zernikes...\n')
%Calculate the unit zernikes over both full and sub pupils
[zernF, inCircF] = unitZern(files, ApproxOrd);         %zernF and inCircF are the values for the full pupil
[zernS, inCircS] = unitZern(subApfiles, ApproxOrd);    %zernS and inCircS are the values for the sub pupil

fprintf('Decomposing and Subtracting Zernikes...\n');
%Perform subtraction and Averaging________________
newSurf = NaN(croplength+1, croplength+1,nfiles);      %preallocates space for P/T/T corrected images
newSurfS = NaN(subCropL+1, subCropL+1,nfiles);
%loop through Zer. Decomp. and subtraction               
for n = 1:nfiles
    newSurf(:,:,n) = zernSubt(files(:,:,n), zernF, inCircF, modRemove);            %Corrected full image surface
    newSurfS(:,:,n)  = zernSubt(subApfiles(:,:,n), zernS, inCircS, modRemove);    %Corrected sub aperture surface
end
 
%% Process, Display, and Save the Data
fprintf('Analyzing the Results...\n')
%Calculate the mean of the images
meanSurf = mean(newSurf,3,'omitnan');       %ignores nan values in averaging. Useful for averaging out holes within maps.
meanSurfS = mean(newSurfS,3,'omitnan');

%show the meanSurface
figure
%mesh(meanSurf)
surf(meanSurf)      %Use surf since mesh is too sparse for this set

%Calculate max and min of subApSurf for a uniform colorbar. 
        % MeanSurf will likely saturate with this colorbar but it allows visualization of small pupil
cmin = min(min(meanSurfS));
cmax = max(max(meanSurfS));
caxis manual
caxis([cmin cmax]);
colorbar
zl = zlim;
title(['Full Pupil (' num2str(fullSize, '%4.2f') 'mm)'])
view(0,90)      %reorient image to be face-on for saving
axis equal
%Modify margins so that saved image has minimal white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
shading interp          %Brightens surf by smoothing it out
saveas(gcf, [fldnm fignmFull '.png']) %Save figure if needed

%Plot masked aperture
figure
%mesh(meanSurfS)
surf(meanSurfS)
zlim(zl);
caxis manual
caxis([cmin cmax]);
colorbar
title(['Sub Pupil (' num2str(subSize, '%4.2f') 'mm)'])
view(0,90)      %reorient image to be face-on for saving
axis equal
%Modify margins so that saved image has minimal white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
shading interp          %Brightens surf by smoothing it out
saveas(gcf, [fldnm fignmSub '.png'])  %Save figure if needed

%Calculate rms error over full and sub pupils
vals = ~isnan(meanSurf);    %remove nan values from calculation
subvals = ~isnan(meanSurfS);
rmsFull = sqrt(mean(meanSurf(vals).^2));    %full rms
rmsSubAp = sqrt(mean(meanSurfS(subvals).^2));    %subAperture rms
fprintf('\nResults:\nRMS Full Pupil: %6.4f',rmsFull); 
fprintf('\nRMS Sub Pupil : %6.4f',rmsSubAp); 

%Calculate remaining Zernike terms over the averaged full aperture
[~, cofF] = zernSubt(meanSurf, zernF, inCircF, modRemove);
%Calculate remaining Zernike terms over the averaged sub aperture
[~, cof] = zernSubt(meanSurfS, zernS, inCircS, modRemove);
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
figure
bar(cof(1:21))
title('First 21 Zernike Coefficients Over the Sub Pupil')
xlabel({'Coefficient','*does not follow regular order'})
ylabel('RMS Contribution (nm)')
xlim([0 22])
grid on
saveas(gcf, [fldnm fignmZern '.png']);

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
anaT(1,1:2)   = ['RMS Flatness: Calculated Full (nm)' cellstr(num2str(rmsFull))];
anaT(2,1:2)   = ['RMS Flatness: Calculated Sub (nm)' cellstr(num2str(rmsSubAp))];
anaT(3,1:2)   = ['RMS Flatness: Zern Approx. (nm)' cellstr(num2str(rmsCof))];
anaT(4,1:2)   = ['RMS Flatness: 1st Coma Total (nm)' cellstr(num2str(com1Cof))];
anaT(5,1:2)   = ['RMS Flatness: 2nd Coma Total (nm)' cellstr(num2str(com2Cof))];

%Provide Analysis Parameters
anaT(6,1:2)   = ['Sample File' cellstr(fnms(1,1:end))];
anaT(7,1:2)   = ['# Samples' cellstr(num2str(nfiles))];
anaT(8,1:2)   = ['Crop Regions:' {'_____'}];
anaT(9,1:2)   = ['   rcrop' cellstr(num2str(rcrop))];
anaT(10,1:2)   = ['   ccrop' cellstr(num2str(ccrop))];
anaT(11,1:2)   = ['   croplength' cellstr(num2str(croplength))];
anaT(12,1:2)  = ['Sub Aperture' cellstr(num2str(subAper))];
anaT(13,1:4)    = [{'Zernike Coefficients'} {'_____'} {'_____'} {'_____'}];

xlswrite([fldnm xlnm '.xlsx'], [anaT; cofTab], 'A1:D31')

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
    circ    = rawSurf(inCirc);  %vector of numbers within the pupel
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