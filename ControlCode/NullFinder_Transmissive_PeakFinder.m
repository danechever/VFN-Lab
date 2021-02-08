

%% add the necessary functions
addpath(genpath(['..' filesep 'ControlCode']));
addpath(genpath(['..' filesep 'AnalysisCode' filesep 'AnalysisLib']));
addpath(genpath(['..' filesep '..' filesep 'VFN-Simulations' filesep 'VFNlib']));

close all
clear all

global itr pk_min_hists FMTO


%% PI stage (nulling fiber stage) Setup
fprintf('---Performing PI Stage Setup\n')
% Connect to the MDT Piezos
VFN_setUpBenchPIStages;

%% Femto setup

%~~ FEMTO POWER METER STUFF 
FMTO.isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO.FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
FMTO.Nread   = 100;      % power samples at a given locaiton
% Add delay for settling time of detector
% NOTE::: Chose arbitrary delay values for now
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END FEMTO POWER METER STUFF 

fprintf('---Performing Femto Setup\n')
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
FMTO.gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

% Set FMTO gain to user-provided value (known start value for code)
VFN_FMTO_setGain(FMTO);

fprintf('Current Gain setting: %i\n', FMTO.FMTO_scale)

%% Define General (Constant) Search Parameters

%-- Constraints
A   = [];
b   = [];
Aeq = [];
beq = [];
% Set nonlcon blank so that we can provide 'options'
nonlcon = [];

%-- Set tolerance and step size values
TolFun = 0.001;                  % Smallest gain required
TolX = .005;                     %Smallest step size

%-- Maximum number of function evals in minimizer
nFunEvals = 150;

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

%% Pre-allocate data matrices
% Struct of iteration history values
pk_min_hists.X    = nan(nFunEvals, 3);        % State vector 
pk_min_hists.PWR  = nan(nFunEvals);           % Measured power 
pk_min_hists.SCLS = nan(nFunEvals);           % Scaling factor used 


%% Minimize to Find Peak 
%%--- Define Search Parameters for finding peak

%-- Set initial params: 
fibX0 = 7.375895;%7.343111;%    % Fiber X position [in mm]
fibY0 = 7.096191;%7.141939;%    % Fiber Y position [in mm]
fibZ0 = 9.510887;%9.700686;%     % Fiber Z position [in mm]

%      - Peak found =  5.727132
%      - Ideal pos  =  (7.375895 mm, 7.096191 mm, 9.510887 mm)

%      - Peak found =  5.745480
%      - Ideal pos  =  (7.375983 mm, 7.096441 mm, 9.511552 mm)

%-- Vectorize initial params
X0 = [fibX0, fibY0, fibZ0];

%-- Set Bounds and Constraints
% Bounds
fibBnd = 0.5;                                %Max distance from fibX0/Y0 [in mm]
focBnd = 0.5;                                %Max distance from focZ0 [in mm]
LB = [fibX0-fibBnd;                          %Lower bounds
      fibY0-fibBnd; 
      fibZ0-focBnd];
UB = [fibX0+fibBnd;                          %Upper bounds
      fibY0+fibBnd;
      fibZ0+focBnd];
% Make sure bounds don't go beyond range
LB = max(LB, [PIdevs.fibX_lower;PIdevs.fibY_lower;PIdevs.fibZ_lower]);
UB = min(UB, [PIdevs.fibX_upper;PIdevs.fibY_upper;PIdevs.fibZ_upper]);

%%--- Perform Peak Search
fprintf('\n-- Searching for Peak\n')
fprintf('   Initial guess: (%f mm, %f mm, %f mm)\n', X0(1), X0(2), X0(3));

%-- Define iteration function
    % This moves the actuators to the new position and measures the power
func = @(X) MinFunc_Transmissive_PeakFinder(X, PIdevs);

%-- Set iteration counter for minimizer correctly
itr = 0;        % Iteration marker [nFunEval]

%[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', TolX, 'TolFun', TolFun, 'MaxFunEvals', nFunEvals, 'PlotFcns', @optimplotfval);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%%--- Print results
fprintf('\n      - Peak found =  %f', -1*fval);
fprintf('\n      - Ideal pos  =  (%f mm, %f mm, %f mm)\n', X(1), X(2), X(3));

%-- Display results
figure; plot(pk_min_hists.X(:,1)); title('X(1)');
figure; plot(pk_min_hists.X(:,2)); title('X(2)');
figure; plot(pk_min_hists.X(:,3)); title('X(3)');
figure; plot(pk_min_hists.PWR(:)); title('PWR');
figure; semilogy(pk_min_hists.PWR(:)); title('LOG PWR');


%% Close connection to all devices
% Clean up FMTO
VFN_cleanUpFMTO;
fprintf('\n Femto Gain set to %i\n',FMTO.FMTO_scale);

%Set PI stages to central position
VFN_PIStage_move(PIdevs.fibX, X(1));
VFN_PIStage_move(PIdevs.fibY, X(2));
VFN_PIStage_move(PIdevs.fibZ, X(3));
% Close connection to the PI Stages
VFN_cleanUpPIStages