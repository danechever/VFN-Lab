%{
This version of the code is used for doing the full, multi-wavelength scan.
  It builds off the tests done with the _FiberScan.m code.

This version uses a slightly different method than the last:
  - The focus is determined by scanning for the peak coupling
  - Then only the fiber X,Y are scanned to find the optimal null

At this point, it does not implement the varia stuff but it will soon.
%}

%% add the necessary functions
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

close all
clear all

global s min_hist FMTO_scale

%% Zaber Setup
fprintf('\n-- Initializing Zabers')
VFN_setUpZabers; % Instantiate the zabers

% Relabel the variables for clarity
if exist('ax14', 'var')
    vortX = ax14;
    clear ax14
end
if exist('ax63', 'var')
    vortY = ax63;
    clear ax63
end
if exist('ax93', 'var')
    fibZ = ax93;
    clear ax93
end    

%% MDT (piezo) Setup
fprintf('\n-- Initializing Piezos')
%Create serial port and define communication properties
    %COM# changes between computers; must check in Device Manager
MDT     = serial('COM4', 'BaudRate', 115200, 'InputBufferSize', 1500, ...
    'Terminator', {'CR' 'LF/CR'});
fopen(MDT);           %Connect port to device

%% Femto setup
fprintf('\n-- Initializing Femto')
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
%%% NOTE:::: The value below is hardcoded into MinFunc_FiberScan instead of here
gnFact  = 9.97;

% Starting gain
FMTO_scale = 6;     % The n in: 10^n for the gain setting
% Other settings
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]

% Setup Femto 
VFN_setUpFMTO;

s = VFN_FMTO_setGain(s, FMTO_scale);

fprintf('\nCurrent Gain setting: %i\n', FMTO_scale)

%% Define General (Constant) Search Parameters
%-- Normalization factor for zaber rescaling
nrmFac = 100;

%-- Constraints
A   = [];
b   = [];
Aeq = [];
beq = [];
% Set nonlcon blank so that we can provide 'options'
nonlcon = [];

%-- Set tolerance and step size values
sTol = .01;                     %Smallest step size

%-- Set Options
opts = optimoptions('patternsearch', 'StepTolerance', sTol, 'Display', 'iter');
%opts = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter', 'StepTolerance', sTol);

% -------- Below is from VFN_An_LSQModelFitting ---------- %
% %-- Set special minimization parameters
% % Set nonlcon blank so that we can provide 'options'
% nonlcon = [];
% % Set options:
%     % Display, iters = display status at each iterations
%     % StepTolerance = minimum change in X before exit
%     % MaxIterations = maximum number of iterations before exit
%     % UseParallel = Numerically calculate the gradient using parallel computing
% if isDispItr
%     disptype = 'iter';  % Status at each iteration
% else
%     disptype = 'final'; % Status only at end
% end
% options = optimoptions('fmincon', ...
%                        'Display', disptype, ...
%                        'StepTolerance', 1e-8, ... %2e-3, ...
%                        'MaxIterations', 750, ... % 50, ...
%                        'UseParallel', isparcomp);

% ----- Below is based on https://www.mathworks.com/help/gads/pattern-search-options.html
% -- Plot best function value vs. iteration.
% options = optimoptions('patternsearch','PlotFcn','psplotbestf');

%--- Define raster scan properties
% Raster step size
rstStp = 10;
% Scan close-to full field
distX=5:rstStp:145; distY = distX;
dist = [distX' distY'];

%% Course raster to see donut in general
measScl = CourseRaster_FullScan(MDT, dist);

%% Minimize to Find Peak 
%%--- Define Search Parameters for finding peak
%-- Use peak in course scan as starting point
[pkVal, pkInd] = max(measScl(:));
[pkX, pkY] = ind2sub(size(meas1), pkInd);

%-- Set initial params: 
fibX0 = distX(pkX);             % Fiber X position [in V]
fibY0 = distY(pkY);             % Fiber Y position [in V]
    %NOTE!!: The zaber functions use mm but I rescale here so that the a scalar 
    % step tolerance can be used. (We want ~0.3V step tol for the piezos and ~
    % 3 microns for zabers. By scaling x100, 3 microns = 0.3 and thus we can
    % sharing same step tolerance for both actuator types.
focZ0 = 3.64450*nrmFac;      % Fiber Z position [in mm converted 10s of microns]

%-- Vectorize initial params
X0 = [fibX0, fibY0, focZ0];

%-- Set Bounds and Constraints
% Bounds
fibBnd = 15;                                %Max distance from fibX0/Y0 [in V]
focBnd = 10;                                %Max distance from focZ0 [in 10s of microns]
LB = [fibX0-fibBnd;                          %Lower bounds
      fibY0-fibBnd; 
      focZ0-focBnd];
UB = [fibX0+fibBnd;                          %Upper bounds
      fibY0+fibBnd;
      focZ0+focBnd];
% Make sure bounds don't go beyond range
LB = max(LB, [0;0;0]);
UB = min(UB, [150,150,7.8]);

%%--- Perform Peak Search
fprintf('\n-- Searching for Peak\n')
%-- Define history vector to contain iteration information
min_hist.X      = [];
min_hist.PWR    = [];
min_hist.scales = [];

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFuncPeak_FullScan(X, fibZ, MDT, nrmFac);

[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);

%%--- Print results
fprintf('\n      - Peak found =  %f', fval);
fprintf('\n      - Ideal pos  =  (%5.2f V, %5.2f, %f mm)\n', X(1), X(2), X(3)/nrmFac);

figure; plot(min_hist.X(:,1)); title('X(1)');
figure; plot(min_hist.X(:,2)); title('X(2)');
figure; plot(min_hist.X(:,3)); title('X(3)');
figure; plot(min_hist.PWR); title('PWR');

%% Course raster at ideal focus to get fresh donut
%-- Set focZ to optimal position
VFN_Zab_move(fibZ, X(3)/nrmFac);
%-- Do raster
measScl = CourseRaster_FullScan(MDT, dist);

%% Find donut in image
%fl = fitsread('5Don_780NullNew1_20x20_SEMIREDU.fits');
%[pkVal, pkInd] = max(fl(:));
%[pkX, pkY] = ind2sub(size(fl), pkInd);
%flThs = zeros(size(fl));
%flThs(fl>pkVal/3) = fl(fl>pkVal/3);
%[cent, radi] = imfindcircles(flThs,[4,15], 'Sensitivity',0.95)
%-- Search for circles in the scan
[cent, radi] = imfindcircles(measScl,[4,13], 'Sensitivity',0.97);
%-- Plot results
figure; imagesc(measScl);
viscircles(cent,radi);
%-- Translate center coords from pixels to voltages
cent = [distX(round(cent(1))) distY(round(cent(2)))];
%-- Convert to polar for circular bounds
[pX, pY] = meshgrid(distX-cent(1), distY-cent(2));
[TH, R]  = cart2pol(pX, pY);
%-- Find null within circular bounds to use as starting point for scan
measIN = measScl;
measIN(R>radi*rstStp) = nan;
[nlVal, nlInd] = min(measIN(:));
[nlX, nlY] = ind2sub(size(measIN),nlInd);
nlX = distX(nlX); nlY = distY(nlY);     % Translate to voltages

%% Minimize to Find null 
%%--- Define Search Parameters for finding null
%-- Set initial params (based on polar coords): 
fibR0 = sqrt(nlX^2 + nlY^2);             % Fiber R position [in V]
fibT0 = acos(nlX/fibR0);             % Fiber Theta position [in rad]

%-- Vectorize initial params
X0 = [fibR0, fibT0];

%-- Set Bounds and Constraints
% Bounds
RBnd = rstStp*radi;                    %Max distance from fibR0 [in V]
LB = [0;                          %Lower bounds
      0];
UB = [fibR0+RBnd;                 %Upper bounds
      2*pi];

%%--- Perform Null Search
fprintf('\n-- Searching for Null\n')
%-- Define history vector to contain iteration information
min_hist.X      = [];
min_hist.PWR    = [];
min_hist.scales = [];

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFuncNullPolar_FullScan(X, MDT, cent);

[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);

%%--- Print results
fprintf('\n      - Peak found =  %f', fval);
fprintf('\n      - Ideal pos  =  (R=%5.2f V, Th=%5.2f rad)\n', X(1), X(2));
fprintf('\n      - Ideal pos  =  (X=%5.2f V, Y=%5.2f V)\n', X(1)*cos(X(2)), X(1)*sin(X(2)));

figure; plot(min_hist.X(:,1)); title('fibX');
figure; plot(min_hist.X(:,2)); title('fibY');
figure; plot(min_hist.X(:,3)); title('fibR');
figure; plot(min_hist.X(:,4)); title('fibT');
figure; plot(min_hist.PWR); title('PWR');

%% Close the connections to all devices
%-- Set Femto back to 10^6 
FMTO_scale = 6;
s = VFN_FMTO_setGain(s, FMTO_scale);
fprintf('\n Femto Gain set to %i\n',FMTO_scale);

%-- Set X and Y axes to optimal position
DE2_MDTVol(MDT, X(1), 'x', 0); 
DE2_MDTVol(MDT, X(2), 'y', 0);

fclose(MDT);
delete(MDT);
clear MDT;

%-- Set focZ to optimal position
VFN_Zab_move(fibZ, X(3)/nrmFac);

%-- Clean up
VFN_cleanUpZabers;
if exist('vortX', 'var')
    clear vortX
end
if exist('vortY', 'var')
    clear vortY
end
if exist('fibZ', 'var')
    clear fibZ
end