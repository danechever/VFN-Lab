%{
This version is used for finding the optimal vortex position for a null.
  It builds off the tests done with the _FullScan.m code.

%}

%% add the necessary functions
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

close all
clear all

global s itr nl_min_hists FMTO_scale RBnd 

%% Scan parameters
%-- Initial conditions for minimization
fibX0 = 75;             % Fiber X position [V]
fibY0 = 75;             % Fiber Y position [V]
focZ0 = 3.7200;         % Fiber Z position [mm]
vorX0 = 12.1450;        % Vortex X position [mm]
vorY0 = 12.4450;        % Vortex Y position [mm]

%-- Bound for fiber position in minimization
RBnd  = 10;             % Max fiber distance from center of scan [V]

%-- Bound for raster scan (how far to scan from null found)
rstBnd = 40;            % Max distance from null to scan in raster [V]
rstPts = 15;            % Number of data points in raster

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

zabs.vortX = vortX;
zabs.vortY = vortY;
zabs.fibZ  = fibZ;

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

%% Minimize to Find null 
%-- Pre-allocate data matrices
% Matrix of raster scan images (these will already be transposed)
nl_rasters = nan(length(distY), length(distX));
% Struct of iteration history values
nl_min_hists.X    = nan(nFunEvals, 5);        % State vector 
nl_min_hists.PWR  = nan(nFunEvals,1);           % Measured power 
nl_min_hists.SCLS = nan(nFunEvals,1);           % Scaling factor used 

%%--- Define Search Parameters for finding null
%-- Set initial params (in real axis positions): 
%fibX0 = 75; fibY0 = 75;             % Fiber Center position in X and Y [V]
     % above values are defined at top
cent  = [fibX0, fibY0];             % Save in cent vector
focZ0 = focZ0*nrmFac;               % Fiber focus position [10s of microns]
vorX0 = vorX0*nrmFac;               % Vortex initial X position [10s of microns]
vorY0 = vorY0*nrmFac;               % Vortex initial Y position [10s of microns]
%-- Set starting conditions to center here
fibR0 = 0;                          % Fiber R position [in V]
fibT0 = 0;                          % Fiber Theta position [in rad]

%-- Vectorize initial params
X0 = [fibR0, fibT0, focZ0, vorX0, vorY0];

%-- Set Bounds and Constraints
%RBnd = 10;                  %Max distance from cent [in V]
    % RBnd is defined at top 

%%--- Perform Null Search
fprintf('\n-- Searching for Null\n')
fprintf('   Initial guess: (\n\t\tfibX=%5.2f V, fibY=%5.2f V, \n\t\tfibZ =%7.4f mm, \n\t\tvortX=%7.4f mm, vortY=%7.4f mm)\n', fibX0, fibY0, X0(3)/nrmFac, X0(4)/nrmFac, X0(5)/nrmFac);

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFuncNullPolar_VortScan(X, MDT, zabs, cent, nrmFac);

%-- Set iteration counter for minimizer
itr = 1; 

%[X, fval, exitflag, output] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, output] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', sTol, 'TolFun', 0.001, 'MaxFunEvals', nFunEvals);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%-- Rescale X by nrmFac
X = X./[1, 1, nrmFc, nrmFac, nrmFac];
nlX = X(1)*cos(X(2))+cent(1); nlY = X(1)*sin(X(2))+cent(2);
%%--- Print results
fprintf('\n      - Null found =  %f', fval);
fprintf('\n      - Ideal pos  =  (fibX=%5.2f,, fibY=%5.2f, fibZ=%7.4f, vorX=%7.4f, vorY=%7.4f)\n', nlX, nlY, X(3),X(4),X(5));

figure; plot(nl_min_hists.X(:,1)); title('fibX');
figure; plot(nl_min_hists.X(:,2)); title('fibY');
figure; plot(nl_min_hists.X(:,3)); title('fibZ');
figure; plot(nl_min_hists.X(:,4)); title('vorX');
figure; plot(nl_min_hists.X(:,5)); title('vorY');
figure; plot(nl_min_hists.PWR); title('PWR');
%ylim([0,1])

%% Raster  to see resulting donut
%--- Define raster scan properties
distX  = linspace(nlX-rstBnd, nlX+rstBnd, rstPts); 
distY  = linspace(nlY-rstBnd, nlY+rstBnd, rstPts);
dist   = [distX' distY'];

%-- Move to ideal focus and vortex positions
% Remove backlash
VFN_Zab_move(fibZ, X(3)-0.05);
VFN_Zab_move(vortX, X(4)-0.05);
VFN_Zab_move(vortY, X(5)-0.05);
% Now move
VFN_Zab_move(fibZ, X(3));
VFN_Zab_move(vortX, X(4));
VFN_Zab_move(vortY, X(5));

%--- Do raster
measScl = CourseRaster_FullScan(MDT, dist);
measScl = transpose(measScl);   % Transpose so that axes are physically oriented
figure; imagesc(distX, distY, measScl);   
title('Raster of Ideal Position');
axis image; axis xy;

%% Set final params on devices
%-- Set Femto back to 10^6 
FMTO_scale = 6;
s = VFN_FMTO_setGain(s, FMTO_scale);
fprintf('\n Femto Gain set to %i\n',FMTO_scale);

%-- Set X and Y axes to last optimal position
DE2_MDTVol(MDT, nlX, 'x', 0); 
DE2_MDTVol(MDT, nlY, 'y', 0);

%% Close connection to all devices
fclose(MDT);
delete(MDT);
clear MDT;

%-- Clean up
VFN_cleanUpZabers;
if exist('vortX', 'var')
    clear vortX
    zabs=rmfield(zabs, 'vortX');
end
if exist('vortY', 'var')
    clear vortY
    zabs=rmfield(zabs, 'vortY');
end
if exist('fibZ', 'var')
    clear fibZ
    zabs=rmfield(zabs, 'fibZ');
end
clear zabs