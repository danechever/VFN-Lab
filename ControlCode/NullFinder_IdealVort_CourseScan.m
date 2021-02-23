%{
This version is used for finding the optimal vortex position for a null.
  It builds off the tests done with the _FullScan.m code.

%}

%% add the necessary functions
addpath(genpath(['..' filesep 'ControlCode']));
addpath(genpath(['..' filesep 'AnalysisCode' filesep 'AnalysisLib']));
addpath(genpath(['..' filesep '..' filesep 'VFN-Simulations' filesep 'VFNlib']));

close all
clear all

global itr FMTO nl_min_hists 

%%
%~~ ZABER STUFF 
% Distance to move zabers for backlash removal [in mm]
zabBackPos = [1, 1];        

%~~ END ZABER STUFF 

%~~ PI STUFF 
% Fiber X/Y center of scan in mm
Xcenter =  7.375299; % [mm]
Ycenter =  7.097086; % [mm]

% Fiber scan properties
%To include center values as a point in scan; use EVEN number of points
Xpoints = 10;% number of X points in scan (actual scan will have +1)
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle

% Fiber step sizes in Volts
refStep   = 10; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/1e-3);
%~~ END PIEZO STUFF 

%~~ FEMTO POWER METER STUFF 
FMTO.isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO.FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
FMTO.Nread   = 100;      % power samples at a given locaiton
%Nrate   = 1000;     % rate of scan [samples/second]
% Add delay for settling time of detector
% NOTE::: Chose arbitrary delay values for now
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END FEMTO POWER METER STUFF 

%~~ RED (NORM) POWER METER STUFF 
pmNread = 100;      % Number of samples to take
isPMNorm = true;   % Flag to mark whether Thorlabs PM100D should be read
                        % This flag is useful in case the PM is not in the
                        % system. -9999 will replace the pmRead values.
pmCalWvl = 635; % Wavelength for redPM for normalization
%~~ END RED (NORM) POWER METER STUFF

%% Scan parameters
%-- Initial conditions for minimization
fibX0 = Xcenter;             % Fiber X position [V]
fibY0 = Ycenter;             % Fiber Y position [V]

%% Zaber Setup
fprintf('---Performing Zaber Setup\n')
isZab = true;

%~~ Connect to Zabers to query position
VFN_setUpBenchZabers; % Instantiate and rename the zabers

%- Make sure pmX is out of the way
if isfield(Zabs, 'pmX')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
elseif isPMNorm
    error('PM Norm is requested but no pmX Zaber exists.')
end

%% PI stage (nulling fiber stage) Setup
fprintf('---Performing PI Stage Setup\n')
% Connect to the MDT Piezos
VFN_setUpBenchPIStages;

%% Femto setup
fprintf('---Performing Femto Setup\n')
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
FMTO.gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

% Set FMTO gain to user-provided value (known start value for code)
VFN_FMTO_setGain(FMTO);

fprintf('Current Gain setting: %i\n', FMTO.FMTO_scale)

%% Red Thorlabs PM setup
fprintf('---Performing redPM Setup (if needed)\n')
%-- Connect to red PM only if it will be used
if isPMNorm
    VFN_setUpPM100D
    
    % Set operating wavelength for calibration
    VFN_PM100D_setWavelength(PM, pmCalWvl);
end

%% Prepare scan
%-- Define Fiber scan locations
% Colon indexing for fib X,Y so that even xPoi has xCent as middle point in scan
distX=(Xcenter-StepSize*Xpoints/2):StepSize:(Xcenter+StepSize*Xpoints/2);
distY=(Ycenter-StepSize*Ypoints/2):StepSize:(Ycenter+StepSize*Ypoints/2);

%% Define General (Constant) Search Parameters
%-- Constraints
A   = [];
b   = [];
Aeq = [];
beq = [];
% Set nonlcon blank so that we can provide 'options'
nonlcon = [];

%-- Set tolerance and step size values
TolX = .0005;                     %Smallest step size
TolFun = 0.01;

%-- Maximum number of function evals in minimizer
nFunEvals = 10;

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
dist = [distX' distY'];

%% Minimize to Find null 
%-- Pre-allocate data matrices
% Struct of iteration history values
nl_min_hists.X    = nan(nFunEvals, 3);        % State vector 
nl_min_hists.PWR  = nan(nFunEvals,1);           % Measured power 
nl_min_hists.nullPt = nan(nFunEvals, 2);        % coords of null poit
nl_min_hists.nullVl = nan(nFunEvals,1);     % normalize null depth
nl_min_hists.peakVl = nan(nFunEvals,1);     % normalized radial avg. peak
nl_min_hists.normVl = nan(nFunEvals,1);     % normalization value

%%--- Define Search Parameters for finding null
focZ0 = VFN_PIStage_getPos(PIdevs.fibZ);               % Fiber focus position 
vorX0 = VFN_Zab_getPos(Zabs.vortX);               % Vortex initial X position 
vorY0 = VFN_Zab_getPos(Zabs.vortY);               % Vortex initial Y position 

% Define zabBackPos to always start zabers at same point
zabBackPos = [vorX0 - zabBackPos(1), vorY0 - zabBackPos(2)];

%-- Vectorize initial params
X0 = [focZ0, vorX0, vorY0];

%-- Set Bounds and Constraints
% Bounds
focBnd = 0.5;                                %Max distance from focZ0 [in mm]
LB = X0 - focBnd;
UB = X0 + focBnd;
% Make sure bounds don't go beyond range
LB = max(LB, [PIdevs.fibZ_lower,Zabs.vortX_lower,Zabs.vortY_lower]);
UB = min(UB, [PIdevs.fibZ_upper,Zabs.vortX_upper,Zabs.vortY_upper]);

%%--- Perform Null Search
fprintf('\n-- Searching for Null\n')
%fprintf('   Initial guess: (\n\t\tfibX=%5.2f V, fibY=%5.2f V, \n\t\tfibZ =%7.4f mm, \n\t\tvortX=%7.4f mm, vortY=%7.4f mm)\n', fibX0, fibY0, X0(3)/nrmFac, X0(4)/nrmFac, X0(5)/nrmFac);
fprintf('   Initial guess: (\n\t\tfibX=%f mm, fibY=%f mm, fibZ=%f mm, \n\t\tvortX=%7.4f mm, vortY=%7.4f mm)\n', fibX0, fibY0, X0(1), X0(2), X0(3));

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
%func = @(X) MinFuncNullPolar_VortScan(X, MDT, zabs, cent, nrmFac, RBnd, zabBnd);
func = @(X) MinFunc_Transmissive_VortFocScanner(X, PIdevs, Zabs, PM, dist, pmNread, isPMNorm, pmCalWvl, zabBackPos);
%-- Set iteration counter for minimizer
itr = 1; 

%[X, fval, exitflag, output] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%opts = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter', 'StepTolerance', sTol, 'FunctionTolerance', 0.000001, 'MaxFunctionEvaluations', nFunEvals);
%[X, fval, exitflag, output] = fmincon(func, X0, A, b, Aeq, beq, [], [], nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', TolX, 'TolFun', TolFun, 'MaxFunEvals', nFunEvals, 'PlotFcns', @optimplotfval);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%%--- Print results
% Find iteration with best results
bstind = find(nl_min_hists.PWR == fval);
nlpt = nl_min_hists.nullPt(bstind);
fprintf('\n      - Rel Tint found =  %f', fval);
%fprintf('\n      - Ideal pos  =  (fibX=%5.2f, fibY=%5.2f, fibZ=%7.4f, vorX=%7.4f, vorY=%7.4f)\n', nlX, nlY, X(3),X(4),X(5));
fprintf('\n      - Ideal pos  =  (fibX=%f, fibY=%f, fibZ=%f, vorX=%7.4f, vorY=%7.4f)\n', nlpt(1), nlpt(2), X(1), X(2), X(3));

figure; plot(nl_min_hists.X(:,1)); title('fibZ');
figure; plot(nl_min_hists.X(:,2)); title('vortZ');
figure; plot(nl_min_hists.X(:,3)); title('vortY');
%figure; plot(nl_min_hists.X(:,5)); title('vorY');
figure; plot(nl_min_hists.PWR); title('relTint');

%% close the connections to all devices
% Clean up FMTO
VFN_cleanUpFMTO;
fprintf('\n Femto Gain set to %i',FMTO.FMTO_scale);

if isZab
    %-- Return vortex to center position if needed
    if X(2) ~= VFN_Zab_getPos(Zabs.vortX)
        % Remove Zab (X) backlash
        VFN_Zab_move(Zabs.vortX, zabBackPos(1));
        % Center Zab (X)
        VFN_Zab_move(Zabs.vortX, X(2));
    end
    if X(3) ~= VFN_Zab_getPos(Zabs.vortY)
        % Remove Zab (Y) backlash
        VFN_Zab_move(Zabs.vortY, zabBackPos(2));
        % Center Zab (Y)
    	VFN_Zab_move(Zabs.vortY, X(3));
    end
    
    % Print position of vortex
    fprintf('\nVortX pos: %f  VortY pos: %f',...
        VFN_Zab_getPos(Zabs.vortX), ...
        VFN_Zab_getPos(Zabs.vortY));
    
    
    %-- Move Calibration PM out of the beam (again; just to be safe)
    if isPMNorm
        fprintf('  |  PM Zab pos: %f\n', VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower));
    else
        fprintf('\n');  % terminate line from zaber recentering
    end
    
    %-- Clean up
    VFN_cleanUpZabers;
end

%Set PI stages to central position
VFN_PIStage_move(PIdevs.fibX, nlpt(1));
VFN_PIStage_move(PIdevs.fibY, nlpt(2));
fprintf('\nfibZ pos: %f', VFN_PIStage_move(PIdevs.fibZ, X(1)));
% Close connection to the PI Stages
VFN_cleanUpPIStages

%-- Disconnect from calibration PM
if isPMNorm
    %-- Only disconnect if we connected earlier
    VFN_cleanUpPM100D;
end