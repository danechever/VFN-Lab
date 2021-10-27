%{
Script to compute the core diameter of a fiber imaged on the microscope


This script automatically detects the imaged cladding and core as circles
in the picture. It assumes the vendor spec for the cladding diameter is
correct and uses this to define the pixelscale for the image. This
pixelscale is then used to compute the core diameter.

%}

close all; clear all;

%% User inputs

%-- File to consider
fldnm = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/FiberTipImaging/VFN_Dan_FiberTipImaging_SaveAs';
flnm = 'AutoWhiteBalance_DFDarkField_FiberCapped.bmp';

%-- Vendor-specific cladding diameter
diam_spec_clad = 125;   % [um] spec for SM600 for thorlabs website (+/- 1um according to website)

%-- Value to use for threshold when isolating cladding in the image
thres_clad = 170;

%-- (OPTIONAL) crop the top N rows to remove text in image
isTxtTopCrop = true;   % Flag to select if cropping should occur
txtTopCrop   = 60;     % Number of rows to crop from the top

%-- (OPTIONAL) crop the bottom N rows to remove text in image
isTxtBotCrop = true;   % Flag to select if cropping should occur
txtBotCrop   = 60;     % Number of rows to crop from the top

%-- Subwindow size to use when finding the core
corecrop_sz = 300;  % [pix] subwindow will be [corecrop_sz X corecrop_sz]

%-- Value to use for threshold when isolating core in the subwindowed image
thres_core = 130;

%% Load the file
fibtipim = imread(fullfile(fldnm,flnm));

%-- (ONLY IF RELEVANT) crop the top of the image which contains some text
if isTxtTopCrop
    fibtipim = fibtipim(txtTopCrop:end,:,:);
end

%-- (ONLY IF RELEVANT) crop the bottom of the image which contains some text
if isTxtBotCrop
    fibtipim = fibtipim(1:end-txtBotCrop,:,:);
end


%% Clean the image to make it easier to find the cladding
%-- Convert image to grayscale to simplify things
fibtipim_clean = rgb2gray(fibtipim);

%-- Threshold so that circle finder works better
fibtipim_binary = fibtipim_clean < thres_clad;

%-- Display the thresholded image
%figure(); imagesc(fibtipim_binary); axis image;

%% Find the fiber cladding in the image

%-- Method 1: imfindcircles
  % This method was very finicky and required lots of tuning to get it to
  % work correctly. The keys to success were: (1) set 'Sensitivity' very
  % high (>0.99), (2) Setting the right ObjectPolarity, and (3) making sure
  % to provide radii NOT diameter values for second parameter.
% [cents_clad, radii_clad] = imfindcircles(fibtipim_binary,[700,800],'Sensitivity',0.992,'ObjectPolarity','bright');

%-- Method 2: regionprops
  % This method was much more robust and work almost immediately. However,
  % it is also much more sensitive and will find many more circles than
  % desired. keys to getting a good result with this one were: (1)
  % thresholding the image so that the circle to find is 1 while outside
  % the circle is 0. 
% Find circles
stats_clad = regionprops('table', fibtipim_binary, 'Centroid', 'Eccentricity', 'EquivDiameter');
% Filter for found objects that are more circular using eccentricity
thres_ecc = .3;
stats_clad( stats_clad.Eccentricity > thres_ecc , : ) = [];
% Filter for found objects with large enough diameter to be our fiber
thres_dia = 1000;
stats_clad( stats_clad.EquivDiameter < thres_dia , : ) = [];
% Place the results in aptly-named variables
cents_clad = stats_clad.Centroid;
radii_clad = stats_clad.EquivDiameter/2;

%-- Method 3: having user draw a line marking the diameter of the cladding
  % This method is very straightforward but relies on the user succesfully
  % selecting the edges of the circle and properly bisecting the circle (so
  % that the line is the diameter and not just a chord).
% figure(); imagesc(fibtipim_clean); axis image; 
% title('Click&drag to draw line marking the diameter of the fiber cladding');
% d = drawline;
% disp('Click&drag to draw line marking the diameter of the fiber cladding')
% pos = d.Position;
% del_pos = diff(pos);
% cents_clad = [(pos(1,1) + pos(2,1))/2, (pos(1,2) + pos(2,2))/2];
% radii_clad = hypot(del_pos(1),del_pos(2))/2;


%-- Compute the pixel scale
pixscale = (radii_clad*2)/diam_spec_clad;  % [pix / um]

%% Display the cladding circle overlayed on the image

%-- Rescale the axes to um and center on the cladding
xax = ((1:size(fibtipim_clean,2)) - cents_clad(1)) / pixscale;
yax = ((1:size(fibtipim_clean,1)) - cents_clad(2)) / pixscale;

%-- Make figure
% Show the image
figure(); imagesc(xax,yax,fibtipim_clean); axis image;
% Add the recognized circle 
  % (use 0,0 since the axes are now centered on the cladding)
viscircles([0,0], radii_clad/pixscale);
% Format figure
xlabel('Spatial [um]')
ylabel('Spatial [um]')
tit_ln1 = sprintf('Cladding Diameter: %0.2f pix (%0.2f um spec)',radii_clad*2,diam_spec_clad);
tit_ln2 = sprintf('PixScale = %0.4f [pix/um]',pixscale);
title({tit_ln1,tit_ln2});

drawnow;
%% Clean the image to make it easier to find the core
%-- Subwindow around the core
  % I assume that the core is near the center of the cladding s.t. we
  % subwindow around that region
rowcrop = round([cents_clad(2)-corecrop_sz/2, cents_clad(2)+corecrop_sz/2]);
colcrop = round([cents_clad(1)-corecrop_sz/2, cents_clad(1)+corecrop_sz/2]);
fibtipim_core = fibtipim_clean(rowcrop(1):rowcrop(2) , colcrop(1):colcrop(2));

%-- Threshold so that circle finder works better
fibtipim_core_binary = fibtipim_core < thres_core ;

%-- Display the thresholded image
%figure(); imagesc(fibtipim_core_binary); axis image;

%% Find the fiber core in the subwindowed image
% Find circles
stats_core = regionprops('table', fibtipim_core_binary, 'Centroid', 'Eccentricity', 'EquivDiameter');
% Filter for found objects that are more circular using eccentricity
thres_ecc = .5;  % Use 0.5 here to be more liberal on the circularity requirement since the core is poorly resolved
stats_core( stats_core.Eccentricity > thres_ecc , : ) = [];
% Filter for found objects with large enough diameter to be our fiber
thres_dia = 5;
stats_core( stats_core.EquivDiameter < thres_dia , : ) = [];
% Place the results in aptly-named variables
cents_core = stats_core.Centroid;
radii_core = stats_core.EquivDiameter/2;

%% Compute the core diameter

%-- Core diameter
daim_core = radii_core*2 / pixscale;

%-- Compute center location in larger cladding image
cents_core_mainim = [cents_core(1) + colcrop(1) - cents_clad(1), cents_core(2) + rowcrop(1) - cents_clad(2)];


%% Display the core circle overlayed on the subwindowed image
figure(); imagesc(fibtipim_core); axis image;
viscircles(cents_core, radii_core);
title('Core in subwindowed image')

drawnow;
%% Display the core circle overlayed on the main image
% Show the image
figure(); imagesc(xax,yax,fibtipim_clean); axis image;
% Add the core circle 
viscircles(cents_core_mainim/pixscale, radii_core/pixscale,'LineWidth', 0.5);
% Add the cladding circle 
  % (use 0,0 since the axes are now centered on the cladding)
viscircles([0,0], radii_clad/pixscale,'LineWidth', 2);
% Format figure
xlabel('Spatial [um]')
ylabel('Spatial [um]')
tit_ln1 = sprintf('Cladding Diameter: %0.2f pix (%0.2f um spec)',radii_clad*2,diam_spec_clad);
tit_ln2 = sprintf('Core Diameter: %0.4f pix (%0.4f um)',radii_core*2, daim_core);
tit_ln3 = sprintf('PixScale = %0.4f [pix/um]',pixscale);
title({tit_ln1,tit_ln2,tit_ln3});


%tit_ln1 = sprintf('Cladding Diameter: %0.2f pix (%0.2f um spec)',radii_clad*2,diam_spec_clad);
%tit_ln2 = sprintf('PixScale = %0.4f [pix/um]',(radii_clad*2)/diam_spec_clad);
%title({tit_ln1,tit_ln2});