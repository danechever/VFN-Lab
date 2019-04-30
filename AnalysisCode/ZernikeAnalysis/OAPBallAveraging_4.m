%{
Script used to process the Zygo Data. Removes Piston, Tip, and Tilt from
the data. Also averages multiple samples and computes the RMS over the full
area and over a subaperture. Creates a plot of both the full and the
subaperture. Finally, returns the Zernike coefficients over the subaperture
in a bar graph as well as in an excel file with pertinent analysis 
information and results. The coefficients have been normalized such that 
they represent the RMS contribution of each mode.

*** MAJOR OVERHAUL IN VERSION _4:
    - Controllable input parameters have all been brought to the top:___
        s.t. all user modifications between data sets occur in one
        section. User now has easy control over the data names, cropping 
        approximations, aperture size, approximation order, and output file
        names. Code also shifts easily between ball averaging and flat
        averaging modes using isOAP parameter.
    * Code automatically finds the circle over which to analyze: ____
        and usesthis to define the crop region. User clicks near the center 
        of the optic in the image and code uses imfindcircl() to determine the
        center and radius. A figure displaying the resulting crop region is
        displayed. This is robust enough to find optics that are partially
        obscured or not completely round as well as images where multiple
        optics are present.
    * Non-reflective regions in the image are completely ignored: ____
        If the optic is partially obscured or the data has missing sections
        within the reflective region due to poor zygo results, these
        regions are completely ignored by the calculation. This was
        PURPOSELY INCLUDED TO IMPROVE THE QUALITY OF THE OF THE
        APPROXIMATION. It was found that these regions would lead to poor
        zernike decomposition before since they represented a 0 in the data
        and thus a huge local piston. The new code ignores this.
       -This also provides better maps when there are non-reflective
        regions in the optic surface as seen in the first OAP4 data.
       -The practical effect is that the user need not be as careful about
        cropping the region smaller than the surface to reduce the effect
        of extraneous points and outliers
    - An "edge" variable was added: ____
        which allows the user to crop down the analysis region slightly along 
        the whole circle. This is useful for data where the edges show sharp 
        discontinuity as is often seen in Zygo data.
    * All output data is now saved into an XLS file: ______
        Instead of printing results into the command window, they are
        stored in an XLS file automatically, for easy access later
       -The XLS file is now unique to each data set, instead of using the
        one for all analyses. It is stored in the same directory as the figures.
       -The XLS file also contains information about the analysis such as
        the crop region, size, filename, etc. to allow easy identification
        and replication of results.
    - Code provides status and feedback in the command window in case of
        long runtimes to give user input on current process.

*** Change for From OAP2 to this version
    ** [no longer prints this value but saves it in xls file - v. _4] **
    [Prints] a third value, 'ans', which is the RMS value as calculated by
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
    - This version does the sub aperture analysis completely separately
    from the full analysis. The old codes simply cropped the results of the
    full image down to the sub aperture. However, this means any local
    piston/T/T (primarily piston) will still be in the value. Thise will
    not actually be seen by the beam so the zernike decomposition should be
    performed locally on the sub aperture to remove P/T/T from this region
    alone. 

%ASSUMES FITS NAMES FOLLOW FORMAT (check in code for exact format if uncertain)
        Code was changing too much between versions to ensure that this
        section was up to date...

%Crops image to specific size and location

%Allows for circular mask to be applied for subaperture analysis. Circle
origin will be at center of analysis image.

%NEED TO CONFIRM: ***[CONFIRMED - v. _4]*** 
    the vector returned by zernfun2() is assumed to only include values 
    within the unit circle. That is, it is the values of the zernike modes 
    within r<=1, in some order. This was determined by looking at the code 
    and seeing that it only deals with inCirc whenever working with the 
    columns of the zernike matrix. 
******************************************************
- Arguments:
    fldnm           Folder containing the data 
                        DO NOT INCLUDE "SEMI_REDU\" ADDRESS IN THIS NAME.
    fnm             Name of files to analyze. Should look like: 
                        "Semi_redu\Final_Flatnes..."
    isOAP           Flag Identifying whether ball indexing should be performed
                        0 : No ball indexing (Flats being analyzed)
                                "_00#" will be appended to filename and # 
                                will increase directly between iterations. 
                                "_001, _002, _003,..."
                        1 : Ball Indexing (OAP's being analyzed)
                                "$_00#" will be appended to fnm. $ will
                                increase when # resets. # will alternate
                                between 1 and 2.
                                "1_001, 1_002, 2_001, 2_002, 3_001,..."
    nfiles          Total number of files to analyze
    subAper         Fraction of full aperture radius over which sub analysis 
                        is performed on. 
    edge            Fudge factor on circle radius. Allows user to shrink
                        radius by "edge" pixels and therefore the crop
                        region by 2*"edge" pixels. This is good for
                        removing the sharp peaks often present at the edge
                        of an optic.
    rSearch         Fudge factor on search region of circle radius.
                        EXTREMELY USEFUL IF CODE THROWS AN INDEXING ERROR
                        AT "rcrop = centers(2)-radi" while calculating the region...
                        - Increase the value and run code again.
    ApproxOrd       Number of zernike modes used in decomposition. 
                        Actual decomposition will be ApproxOrd+1 modes.
                        36 recommended, higher orders increase runtime.
    fullSize        Actual size of the full optic. This will be displayed
                        in the full map figure and saved file.
    subSize         Same as fullSize but for the subaperture figure
    fignmFull       Name of saved .png file containing full map
    fignmSub        Name of saved .png file containing sub aperture map
    fignmZern       Name of saved .png file containing the bar graph of
                        zernike coefficients
    xlnm            Name for the .xls file containing computational results
    
    *** CROP CAN BE MANUALLY ADJUSTED TO OVERRIDE AUTOMATIC CIRCLE CALCULATION
          *must be a square centered over the circular apetrure to analyze. 
                    Cropped region will be: 
                        row = [rcrop:rcrop+croplength], col = [ccrop:ccrop+croplength]
    rcrop           Row to start cropping on
    ccrop           column to start cropping on
    croplength      number of rows/columns to include

- Returns:
    NONE            Plots a figure showing a single full image and the
                        corresponding circle/crop region
                    Plots the resulting OAP flatness over the full and sub
                        apertures
                    Plots a bar graph containing the RMS contribution of
                        each zernike mode
                    Saves a .xls file containing RMS values and analysis
                        parameters. RMS values in nm.

            Readily available variables:
                - RMS over full and Sub ap., RMS over sub ap. calculated from the zernike approximation and coefficients (rms...)
                - Raw and subtracted surface maps, average maps (...Surf...)
                - Unit zernike map over calculated region (zern...)
                - Zernike Coefficients of sub ap (cof). Easy enough to modify to obtain other coefficients.
                - Unit circle region indexes (inCirc...)
- Dependencies:
    zernfun2_DE     Function to calculate Zernike modes
                        Can be downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/7687-zernike-polynomials
******************************************************

Initial Version:    Daniel Echeverri, 10/20/2017
Last Edit:          Daniel Echeverri
Last Modified:      11/17/2017
%}

%% Parameters to set
%File/Analysis Parameters
fldnm   = 'C:\Users\Daniel\Desktop\Mawet201718\HCST_optics\Flat_Mirror_1\';  %Folder containing data
fnm     = 'Semi_redu\Final_flatness';    %Filename for data
    %Assumes all filenames are index-1. If data is index-0, minor changes will need to be made in the code itself.
isOAP   = 0;        %Flag to mark that the fnm has "_ball" in it and needs to iterate through balls as well as samples
                        %0: adds only a 00# at the end; 1: adds "#_00#"
nfiles  = 10;        %Number of files to analyze
subAper = 36.33/50.00;       %fraction of full aperture radius over which to calculate rms
edge    = 5;        %amount of edge to remove. Crops the "full" image by shrinking the radius to remove some edge effects 
rSearch = 50;       %Starting radius leeway. Radius search will be calculated value +- rSearch
                        %^^ Very useful if code fails with "Index exceeds matrix dimensions." becuase of "rcrop = centers(2)-radi"
                        %25 is usually good but 50 can be better for larger pupils
ApproxOrd = 36;      %Number of zernike modes used in decomposition/approximation

%Naming Parameters
fullSize= 36.33;    %Size of full optic (in mm)
subSize = 50.00;    %Size of subaperture (in mm)
fignmFull   = 'Full_Aperture_Tst4';
fignmSub    = 'Sub_Aperture_Tst4';
fignmZern   = 'ZernCoeff_Tst4';
xlnm        = 'ZernCoeff_Tst4';
%IF NEEDED, YOU CAN MODIFY THE CROP PARAMETERS MANUALLY FURTHER DOWN

%% Begin Code_______________
m2nm    = 1e9;      %meter to nanometer conversion 

%% Calculate Analysis Region
%Import/display Data
if ~isOAP
    fnm1    = [fldnm fnm '_' num2str(1, '%03.0f') '.fits'];
else
    fnm1    = [fldnm fnm num2str(1) '_' num2str(1, '%03.0f') '.fits'];
end
fits    = fitsread(fnm1);
mesh(fits)
view(0,90)
axis equal

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

fprintf('Importing Files...\n')
%% Load Files
%preallocates space for imported files
files = zeros(croplength+1,croplength+1,nfiles);            %Full image files
subCropL = floor(croplength*subAper);                              %width of sub aperture
subrcrop = rcrop+(floor((croplength-subCropL)/2));                 %row crop start for sub ap
subccrop = ccrop+(floor((croplength-subCropL)/2));                 %col crop start for sub ap
subApfiles = zeros(subCropL+1,subCropL+1,nfiles);       %Sub Aperture files
%loop to iterate through opening and cropping the files
for n = 1:nfiles
    %Create string for file name; CHANGE TO MATCH LOCATION OF FILE
    if ~isOAP
        fnm1    = [fldnm fnm '_' num2str(n, '%03.0f') '.fits'];
    else
        i = floor((n+1)/2);     %ball index; corresponds to each ball rotation
        j = mod(n,2)+1;         %sample index; corresponds to each sample in a ball rotation
        fnm1    = [fldnm fnm num2str(i) '_' num2str(j, '%03.0f') '.fits'];
    end
    fits = fitsread(fnm1);
    %Crop file and save it into the existing files stack
    files(:,:,n) = fits(rcrop:rcrop+croplength, ccrop:ccrop+croplength);    
    subApfiles(:,:,n) = fits(subrcrop:subrcrop+subCropL, subccrop:subccrop+subCropL);
end

files(~files) = NaN;        %Assume all 0 values are not real data and set to NaN to remove from analysis
subApfiles(~subApfiles) = NaN;
                                %odds of a point being at exactly 0 are really slim...
files = files*m2nm ;        %convert all phase data to nm
subApfiles = subApfiles*m2nm;

%% Perform Zernike Analysis
fprintf('Creating Unit Zernikes...\n')
%Calculate the unit zernikes over both full and sub apertures
[zernF, inCircF] = unitZern(files, ApproxOrd);         %zernF and inCircF are the values for the full aperture
[zernS, inCircS] = unitZern(subApfiles, ApproxOrd);    %zernS and inCircS are the values for the sub aperture

fprintf('Decomposing and Subtracting Zernikes...\n');
%Perform subtraction and Averaging________________
newSurf = NaN(croplength+1, croplength+1,nfiles);      %preallocates space for P/T/T corrected images
newSurfS = NaN(subCropL+1, subCropL+1,nfiles);
%loop through Zer. Decomp. and subtraction of 1st 3 modes               
 for n = 1:nfiles
     newSurf(:,:,n) = zernSubt(files(:,:,n), zernF, inCircF, 3);            %Corrected full image surface
     newSurfS(:,:,n)  = zernSubt(subApfiles(:,:,n), zernS, inCircS, 3);    %Corrected sub aperture surface
 end

%% Process, Display, and Save the Data
fprintf('Analyzing the Results...\n')
%Calculate the mean of the images
meanSurf = mean(newSurf,3,'omitnan');       %ignores nan values in averaging. Useful for averaging out holes within maps.
meanSurfS = mean(newSurfS,3,'omitnan');

%show the meanSurface
figure
mesh(meanSurf)
%surf(meanSurf)      %Use surf since mesh is too sparse for this set

%Calculate max and min of subApSurf for a uniform colorbar. 
        % MeanSurf will likely saturate with this colorbar but it allows visualization of small pupil
cmin = min(min(meanSurfS));
cmax = max(max(meanSurfS));
caxis manual
caxis([cmin cmax]);
colorbar
zl = zlim;
title(['Full Aperture (' num2str(fullSize, '%4.2f') 'mm)'])
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
%shading interp          %Brightens surf by smoothing it out
saveas(gcf, [fldnm fignmFull '.png'])

%Plot masked aperture
figure
mesh(meanSurfS)
%surf(meanSurfS)
zlim(zl);
caxis manual
caxis([cmin cmax]);
colorbar
title(['Effective Aperture (' num2str(subSize, '%4.2f') 'mm)'])
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
%shading interp          %Brightens surf by smoothing it out
saveas(gcf, [fldnm fignmSub '.png'])

%Calculate rms error over full and masked aperture
vals = ~isnan(meanSurf);    %remove nan values from calculation
subvals = ~isnan(meanSurfS);
rmsFull = sqrt(mean(meanSurf(vals).^2));    %full rms
rmsSubAp = sqrt(mean(meanSurfS(subvals).^2));    %subAperture rms

%Calculate remaining Zernike terms over the averaged sub aperture
[~, cof] = zernSubt(meanSurfS, zernS, inCircS, 3);
%Use coefficients to calculate the RMS error as provided by zernike approximation
    %Sum in quadrature
rmsCof = sqrt(sum(cof.^2));


%Plot bar graph of first 21 Zernikes
figure
bar(cof(1:21))
title('First 21 Zernike Coefficients Over the Effective Aperture')
xlabel({'Coefficient','*does not follow regular order'})
ylabel('RMS Contribution (nm)')
xlim([0 22])
grid on
saveas(gcf, [fldnm fignmZern '.png']);

%Create table of first 15 Zernikes
cofThres = .5;
cofT    = cof(1:15);
cofi    = abs(cofT)<cofThres;
cofT    = cellstr(num2str(cofT, '%+06.2f'));
cofT(cofi)    = {'------'};
%Provide formatting for table
cofT1   = {'N index' '0' '1' '1' '2' '2' '2' '3' '3' '3' '3' '4' '4' '4' '4' '4'};
cofT2   = {'M index' '0' '1' '-1' '2' '0' '-2' '3' '1' '-1' '-3' '4' '2' '0' '-2' '-4'};
cofT3   = {'Name' 'Piston' 'Tip' 'Tilt' 'Astig - Y' 'Focus' 'Astig - 45' 'Trefoil - 45' ...
    'Coma - X' 'Coma - Y' 'Trefoil - Y' 'QdFoil - Y' '2nd Astig - Y' 'Spherical' ...
    '2nd Astig - 45' 'QdFoil - 45'};
cofTab = [cofT1' cofT2' cofT3' ['RMS Contribution (nm)'; cofT]];

%Provide RMS Values
anaT    = cell(11,4);
anaT(1,1:2)   = ['RMS Flatness: Calculated Full (nm)' cellstr(num2str(rmsFull))];
anaT(2,1:2)   = ['RMS Flatness: Calculated Sub (nm)' cellstr(num2str(rmsSubAp))];
anaT(3,1:2)   = ['RMS Flatness: Zern Approx. (nm)' cellstr(num2str(rmsCof))];

%Provide Analysis Parameters
anaT(4,1:2)   = ['Sample File' cellstr(fnm1)];
anaT(5,1:2)   = ['# Samples' cellstr(num2str(nfiles))];
anaT(6,1:2)   = ['Crop Regions:' {'_____'}];
anaT(7,1:2)   = ['   rcrop' cellstr(num2str(rcrop))];
anaT(8,1:2)   = ['   ccrop' cellstr(num2str(ccrop))];
anaT(9,1:2)   = ['   croplength' cellstr(num2str(croplength))];
anaT(10,1:2)  = ['Sub Aperture' cellstr(num2str(subAper))];
anaT(11,1:4)    = [{'Zernike Coefficients'} {'_____'} {'_____'} {'_____'}];

xlswrite([fldnm xlnm '.xlsx'], [anaT; cofTab], 'A1:D27')

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

%function that decomposes and subtracts out first a given number of modes
function [correctedSurf, coeffs] = zernSubt(rawSurf, zerns, inCirc, subOrder)
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
    subt(inCirc) = zerns(:, 1:subOrder)*coeffs(1:subOrder);     %Calculate by multiplying the coeffs with the unit Zernikes
    
    %Subtract the undesired modes from the original image
        %This subtraction implicitly removes values outside the analysis
        %circle since: number - NaN = NaN 
    correctedSurf = rawSurf - subt;
end