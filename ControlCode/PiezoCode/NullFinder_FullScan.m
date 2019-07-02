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

global s itr pk_min_hists nl_min_hists FMTO_scale RBnd 

%% Run Parameters

%*** NOTE: make sure you've turned on and set the power level on the superK
%   before running this script. The script will take care of wavelength stuff
%   but it assumes the superK is already on and at the desired power.
% -- Use "VFN_NKT_quickSet" if needed

isvaria = false;         % Flag to use varia/superK
lams    = [633, 780];   % Wavelengths to use
bwd     = 3;            % Bandwidth to use (nm)

%% Varia/SuperK setup
if isvaria
    fprintf('\n-- Initializing Varia')
    VFN_setUpNKT;   % Instantiate the NKT (will be "NKT" variable)
end
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
sTol = .2;                     %Smallest step size

%-- Maximum number of function evals in minimizer
nFunEvals = 250;

%-- Set Options
%opts = optimoptions('patternsearch', 'StepTolerance', sTol, 'Display', 'iter');
%opts = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter', 'StepTolerance', sTol, 'FunctionTolerance', 0.01);
% ** Options are defined individually per minimizer below

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

%% Prepare for main loop (iterating through wavelengths)
if ~isvaria
    % Set lams to nan since we won't be using the varia but need single loop itr
    lams = [nan];
end

%-- Pre-allocate data matrices
nlam = length(lams);
% Matrix of raster scan images (these will already be transposed)
pk_rasters = nan(nlam, length(distY), length(distX));
nl_rasters = pk_rasters;
% Struct of iteration history values
pk_min_hists.X    = nan(nlam, nFunEvals, 3);        % State vector (peak scan)
pk_min_hists.PWR  = nan(nlam, nFunEvals);           % Measured power (peak scan)
pk_min_hists.SCLS = nan(nlam, nFunEvals);           % Scaling factor used (peak scan)
nl_min_hists.X    = nan(nlam, nFunEvals, 4);        % State vector (null scan)
nl_min_hists.PWR  = pk_min_hists.PWR;               % Measured power (null scan)
nl_min_hists.SCLS = pk_min_hists.SCLS;              % Scaling factor used (null scan)
% Minimization solutions
pk_fval = nan(nlam, 1);        % Peak found
pk_X    = nan(nlam, 3);     % Ideal state vector    [X, Y, Z]
% Minimization solutions
nl_fval = nan(nlam, 1);        % Null found
nl_X    = nan(nlam, 2);     % Ideal state vector    [X, Y]

%% Main loop
for iii = 1:nlam
if isvaria
    % Set varia wavelength
    VFN_NKT_setWvlRange(NKT, lams(iii)-bwd/2, lams(iii)+bwd/2);
    fprintf('\n-- Wavelength set to %f nm\n', lams(iii));
end
%% Course raster to see donut in general
measScl = CourseRaster_FullScan(MDT, dist);
pk_rasters(iii,:,:) = transpose(measScl);
figure; imagesc(distX, distY, transpose(measScl));    % Transpose so that axes are physically oriented
title('Initial Raster - for Peak');
axis image; axis xy;

%% Minimize to Find Peak 
%%--- Define Search Parameters for finding peak
%-- Use peak in course scan as starting point
[pkVal, pkInd] = max(measScl(:));
[pkX, pkY] = ind2sub(size(measScl), pkInd);

%-- Set initial params: 
fibX0 = distX(pkX);             % Fiber X position [in V]
fibY0 = distY(pkY);             % Fiber Y position [in V]
    %NOTE!!: The zaber functions use mm but I rescale here so that the a scalar 
    % step tolerance can be used. (We want ~0.3V step tol for the piezos and ~
    % 3 microns for zabers. By scaling x100, 3 microns = 0.3 and thus we can
    % sharing same step tolerance for both actuator types.
focZ0 = VFN_Zab_getPos(fibZ)*nrmFac;      % Fiber Z position [in mm converted 10s of microns]

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
UB = min(UB, [150;150;7.8*nrmFac]);

%%--- Perform Peak Search
fprintf('\n-- Searching for Peak\n')
fprintf('   Initial guess: (%6.2f V, %6.2f V, %7.4f mm)\n', X0(1), X0(2), X0(3)/nrmFac);

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFuncPeak_FullScan(X, fibZ, MDT, nrmFac);

%-- Set iteration counter for minimizer correctly
itr = [iii, 0];        % Iteration marker [itr(1) = nlam, itr(2) = nFunEval]

%[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', 0.075, 'TolFun', 0.00001, 'MaxFunEvals', nFunEvals);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%%--- Print results
fprintf('\n      - Peak found =  %f', fval);
fprintf('\n      - Ideal pos  =  (%5.2f V, %5.2f, %f mm)\n', X(1), X(2), X(3)/nrmFac);

%-- Save results
pk_fval(iii) = fval;
pk_X(iii,:) = X./[1,1,nrmFac];  %nrmFac to account for scaling

%-- Display results
figure; plot(pk_min_hists.X(iii,:,1)); title('X(1)');
figure; plot(pk_min_hists.X(iii,:,2)); title('X(2)');
figure; plot(pk_min_hists.X(iii,:,3)); title('X(3)');
figure; plot(pk_min_hists.PWR(iii,:)); title('PWR');

%% Course raster at ideal focus to get fresh donut
%-- Set focZ to optimal position
VFN_Zab_move(fibZ, X(3)/nrmFac);
%-- Do raster
measScl = CourseRaster_FullScan(MDT, dist);
nl_rasters(iii,:,:) = transpose(measScl);

%% Find donut in image
%fl = fitsread('5Don_780NullNew1_20x20_SEMIREDU.fits');
%[pkVal, pkInd] = max(fl(:));
%[pkX, pkY] = ind2sub(size(fl), pkInd);
%flThs = zeros(size(fl));
%flThs(fl>pkVal/3) = fl(fl>pkVal/3);
%[cent, radi] = imfindcircles(flThs,[4,15], 'Sensitivity',0.95)
%-- Search for circles in the scan
[cent, radi] = imfindcircles(transpose(measScl),[3,13], 'Sensitivity',0.92);
%-- Translate center coords from pixels to voltages
cent = [distX(round(cent(1))) distY(round(cent(2)))];
%-- Convert to polar for circular bounds
[pX, pY] = meshgrid(distX-cent(1), distY-cent(2));
[TH, R]  = cart2pol(pX, pY);
%-- Find null within circular bounds to use as starting point for scan
measIN = transpose(measScl);
measIN(R>radi*rstStp) = nan;
[nlVal, nlInd] = min(measIN(:));
[nlY, nlX] = ind2sub(size(measIN),nlInd);   % reversed order since transposed measIN
nlX = distX(nlX); nlY = distY(nlY);     % Translate to voltages
nlX = nlX - cent(1); nlY = nlY-cent(2);
%-- Plot results
figure; imagesc(distX, distY, transpose(measScl));    % Transpose so that axes are physically oriented
viscircles(cent,radi*rstStp);
axis image; axis xy;
title('Focused Raster - for Null');

%% Minimize to Find null 
%%--- Define Search Parameters for finding null
%-- Set initial params (based on polar coords): 
fibR0 = sqrt(nlX^2 + nlY^2);             % Fiber R position [in V]
if fibR0 == 0
    fibR0 = 0.00001;
end
fibT0 = acos(nlX/fibR0);             % Fiber Theta position [in rad]

%-- Vectorize initial params
X0 = [fibR0, fibT0];

%-- Set Bounds and Constraints
% Bounds
RBnd = 19*rstStp*radi/20;           %Max distance from fibR0 [in V]
%RBnd = min(cat(2, RBnd, cent, [150 150]-cent));  % Keep search region within limits 
LB = [0;                          %Lower bounds
      0];
UB = [RBnd;                 %Upper bounds
      2*pi];

%%--- Perform Null Search
fprintf('\n-- Searching for Null\n')
fprintf('   Initial guess: (R=%6.2f V, TH=%6.2f rad, X=%6.2f V, Y=%6.2f V)\n', X0(1), X0(2), X0(1)*cos(X0(2))+cent(1), X0(1)*sin(X0(2))+cent(2));

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFuncNullPolar_FullScan(X, MDT, cent);

%-- (re)Set iteration counter for minimizer correctly
itr = [iii, 0];        % Iteration marker [itr(1) = nlam, itr(2) = nFunEval]

%[X, fval, exitflag, output] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, output] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', sTol, 'TolFun', 0.000001, 'MaxFunEvals', nFunEvals);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%%--- Print results
fprintf('\n      - Null found =  %f', fval);
fprintf('\n      - Ideal pos  =  (R=%5.2f V, Th=%5.2f rad)', X(1), X(2));
fprintf('\n      - Ideal pos  =  (X=%5.2f V, Y=%5.2f V)\n', X(1)*cos(X(2))+cent(1), X(1)*sin(X(2))+cent(2));

%-- Save minimization results
nl_fval(iii) = fval;
nl_X(iii,:) =  [X(1)*cos(X(2))+cent(1), X(1)*sin(X(2))+cent(2)];

figure; plot(nl_min_hists.X(iii,:,1)); title('fibX');
figure; plot(nl_min_hists.X(iii,:,2)); title('fibY');
figure; plot(nl_min_hists.X(iii,:,3)); title('fibR');
figure; plot(nl_min_hists.X(iii,:,4)); title('fibT');
figure; plot(nl_min_hists.PWR(iii,:)); title('PWR');
%ylim([0,1])
end

%% Set final params on devices
%-- Set Femto back to 10^6 
FMTO_scale = 6;
s = VFN_FMTO_setGain(s, FMTO_scale);
fprintf('\n Femto Gain set to %i\n',FMTO_scale);

%-- Set X and Y axes to last optimal position
DE2_MDTVol(MDT, nl_X(end,1), 'x', 0); 
DE2_MDTVol(MDT, nl_X(end,2), 'y', 0);

%% Close connection to all devices
fclose(MDT);
delete(MDT);
clear MDT;

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