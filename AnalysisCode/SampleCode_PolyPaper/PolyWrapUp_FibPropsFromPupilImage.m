%{

Script for getting the fiber NA and modefield (in pupil-plane) using the
pupil images I took of the reverse-fed system.

%} 
close all; clear all;
% add path to analysis library
addpath('../AnalysisLib')

%% Folder and file locations
% Folder containing the pupil images
fldnm = '/media/Data_Drive/VFN/TestbedData_Processed/PolyPaper/20210726_PupilImages_PolyPaper';


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

%% Do gaussian fit to fiber mode
[params, residuals] = act2DGaussFit(PUPIL_ONoVortNoSat,true);

%% Analyze the results
sigma2D = mean([params(3), params(5)]); % 1STD size of fibermode
% Derivation for diameter from sigma:
  % Refs: https://en.wikipedia.org/wiki/Beam_diameter
  % Refs: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
  % FWHM = 2*sqrt(2*ln(2))*sigma
  % 2w = sqrt(2)*(2*sqrt(2*ln(2)))*sigma / sqrt(ln(2))
  % ---> 2w = 4sigma  where w is the radius of the 1/e2 point
Dfibe2  = 4*sigma2D;                % fibermode diameter to the 1/e2 point
cent0   = [params(2),params(4)];    % center of fibermode relative to origin at center of image
centim  = cent0 + (size(PUPIL_ONoVortNoSat)/2); % center of fibermode relative to corner of image (normal image coordinates)

fprintf('\nMode Diam to 1/e2 point: %f [pix]\n',Dfibe2);

% Get pupil diameter and use it to compute the pixelscale
  % Treat the average cropped image size as the pupil diameter
Dpuppix  = mean(size(PUPIL_ONoVortNoSat));
% Known pupil diameter from knife edge and other measurements
Dpupreal = 2.45;    % [mm]
pixscale = Dpupreal / Dpuppix;  % [mm / pix]

fprintf('Image Pixel scale: %f [mm / pix]\n',pixscale);

% Compute fiber NA. Follow Nem's definition from the OAP coup doc
  % NA = [ (1/e2 diameter) / focal_len ] / 2
foc_len = 11.00;                % [mm]
Dfibe2mm = Dfibe2 * pixscale;   % [mm] fiber diameter in mm
fibNA = (Dfibe2mm / foc_len) / 2;     % 

fprintf('Mode Diam to 1/e2 point: %f [mm]\n',Dfibe2mm)
fprintf('Fiber NA: %f\n', fibNA)

% Compute the ideal pupil diameter based on the fiber mode size
  % Nem's OAP coup doc suggests that the ideal for their fiber's 12.24mm
  % fiber diameter at the pupil is a pupil stop of 13.6mm. Let's assume
  % that ratio holds for us as well.
fib2pupRatio = 12.24/13.6;  % [fiber / pupil]
DpupIdeal = Dfibe2mm / fib2pupRatio;

fprintf('Ideal pupil diameter would be: %f [mm]\n',DpupIdeal);

% Compute the F/# mismatch in the VFN system
FnumIdeal = foc_len / DpupIdeal;
FnumActl  = foc_len / Dpupreal; 

fprintf('Real  Fnum is: %f\n',FnumActl);
fprintf('Ideal Fnum is: %f\n',FnumIdeal);
fprintf('==> The real F# is   %f   times the ideal size\n', DpupIdeal/Dpupreal)
fprintf('==> Our pupil is     %f   times larger than it should be\n', Dpupreal/DpupIdeal)

%-- Display fibermode diameter over image
figure(); imagesc(PUPIL_ONoVortNoSat); axis image; hold on;
viscircles(centim, Dfibe2/2);
scatter(centim(1),centim(2),300,'r+')
%title({'VFN Pupil Image - Vort Out - Not Saturated - Cropped', ' '})
xlabel('Spatial [pix]')
ylabel('Spatial [pix]')
drawnow;

%%
function [params, residuals] = act2DGaussFit(frame, isPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Fit a 2D Gaussian Function to Data
%
%  Adapted from: https://www.mathworks.com/matlabcentral/fileexchange/55033-fit-1d-and-2d-gaussian-to-noisy-data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:  
%   Fit a 2D gaussian centroid to a DM actuator poke
% 
% INPUT:
%   frame: two-dimensional array with poke actuator surface.
%
% OUTPUT: 
%   Gaussian function parameters.
%
% NOTE:
%   1.This routine uses Matlab's 'lsqcurvefit' function to fit.
%   2.The initial values in x0 must be close to x in order for the fit
%   to converge to the values of x (especially if noise is added).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---Fitting Functions---
%
% Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, offset, Angle(in Radians)]
%
% X-data convention:
%	X(:,:,1) : x-coordinates,  
%	X(:,:,2) : y-coordinates.
%
% In this numerical test we use two-dimensional fitting functions:
% 1. 2D Gaussian function ( A requires 5 coefs ).
g = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) ) + A(6);
% 2. 2D Rotated Gaussian function ( A requires 6 coefs ).
f = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(7))-X(:,:,2)*sin(A(7)) - A(2)*cos(A(7))+A(4)*sin(A(7)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(7))+X(:,:,2)*cos(A(7)) - A(2)*sin(A(7))-A(4)*cos(A(7)) ).^2/(2*A(5)^2) ) ) + A(6);
%% ---Parameters---
n = size(frame,2)-1; m = size(frame,1)-1;   % n x m pixels area/data matrix
% assume that max in frame is due to the poked actuator --> Amp = max
% assume "frame" is cropped with the actuator at the center --> x0,y0 = 0
% random guess at width --> guess 3 pixels
params0 = [max(frame(:)),0,1100,0,1100,0,0];   % Inital (guess) parameters
InterpMethod='linear'; % 'nearest','linear','spline','cubic'
FitOrientation='fit';	% 'fit': fit for orientation, 'dont' fit for orientation
%% ---Build numerical Grids---
% Numerical Grid
[x,y]=meshgrid(-n/2:n/2,-m/2:m/2); X=zeros(m+1,n+1,2); X(:,:,1)=x; X(:,:,2)=y;
% High Resolution Grid
h=3; [xh,yh]=meshgrid(-n/2:1/h:n/2,-m/2:1/h:m/2); Xh=zeros(h*m+1,h*n+1,2); Xh(:,:,1)=xh; Xh(:,:,2)=yh;
%% ---Fit---
% Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
lb = [0,-n/2,0,-n/2,0,-inf,0];
ub = [realmax('double'),n/2,(n/2)^2,n/2,(n/2)^2,inf,pi];
opts = optimset('Display','off');
% Fit sample data
switch FitOrientation
    case 'dont', [params,~,residuals,~,output] = lsqcurvefit(g,params0(1:6),X,frame,lb(1:6),ub(1:6),opts);
    case 'fit',  [params,~,residuals,~,output] = lsqcurvefit(f,params0,X,frame,lb,ub,opts);
    otherwise, error('invalid entry');
end
%disp(output); % display summary of LSQ algorithm
%% ---Plot Data---
if isPlot
%    % Plot 3D Data and Fitted curve
%    hf1=figure(); set(hf1,'Position',[1000 600 800 500]); 
%    switch FitOrientation
%        case 'dont', C=del2(g(params,Xh)); mesh(xh,yh,g(params,Xh),C); hold on
%        case 'fit',  C=del2(f(params,Xh)); mesh(xh,yh,f(params,Xh),C); hold on
%    end
%    surface(x,y,frame,'EdgeColor','none'); alpha(0.5); 
%    colormap('pink'); view(-60,20); grid on; hold off

    % Plot Sample Pixels data
    hf2=figure(); set(hf2,'Position',[20 20 800 800]); 
    subplot(4,4,[5,6,7,9,10,11,13,14,15]); imagesc(x(1,:),y(:,1),frame); 
    colormap('hot');
    % Output and compare data and fitted function coefs
    text(-n/2-n/10,m/2+m/10,sprintf('      Amplitude | X-Coord | X-Width | Y-Coord | Y-Width | Offset | Angle'),'Color','black');
    text(-n/2-n/10,m/2+m/7,sprintf('Fit | % 9.4f | % 7.4f | % 7.4f | % 7.4f | % 7.4f | % 7.4f | % 7.4f',params),'Color','red');
    % Plot vertical and horizontal axis
    vx_h=x(1,:); vy_v=y(:,1);
    switch FitOrientation
        case 'fit', M=-tan(params(7));
            % generate points along _horizontal & _vertical axis
            vy_h=M*(vx_h-params(2))+params(4); hPoints = interp2(x,y,frame,vx_h,vy_h,InterpMethod);
            vx_v=M*(params(4)-vy_v)+params(2); vPoints = interp2(x,y,frame,vx_v,vy_v,InterpMethod);
        case 'dont', params(7)=0; 
            % generate points along _horizontal & _vertical axis
            vy_h=params(4)*ones(size(vx_h)); hPoints = interp2(x,y,frame,vx_h,vy_h,InterpMethod);
            vx_v=params(2)*ones(size(vy_v)); vPoints = interp2(x,y,frame,vx_v,vy_v,InterpMethod);
    end
    % plot lines 
    hold on; plot(params(2),params(4),'+b',vx_h,vy_h,'.r',vx_v,vy_v,'.g'); hold off;
    % Plot cross sections 
    dmin=1.1*min(frame(:)); xfit=xh(1,:); hfit=params(1)*exp(-(xfit-params(2)).^2/(2*params(3)^2))+params(6);
    dmax=1.1*max(frame(:)); yfit=yh(:,1); vfit=params(1)*exp(-(yfit-params(4)).^2/(2*params(5)^2))+params(6);
    subplot(4,4,[1,2,3]); xposh = (vx_h-params(2))/cos(params(7))+params(2); 
    plot(xposh,hPoints,'r.',xfit,hfit,'black'); grid on; axis([-n/2,n/2,dmin,dmax]);
    subplot(4,4,[8,12,16]); xposv = (vy_v-params(4))/cos(params(7))+params(4); 
    plot(vPoints,xposv,'g.',vfit,yfit,'black'); grid on; axis([dmin,dmax,-m/2,m/2]); 
    set(gca,'YDir','reverse');
end
end