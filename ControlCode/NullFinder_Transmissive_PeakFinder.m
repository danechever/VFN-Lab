

%% add the necessary functions
addpath(genpath(['..' filesep 'ControlCode']));
addpath(genpath(['..' filesep 'AnalysisCode' filesep 'AnalysisLib']));
addpath(genpath(['..' filesep '..' filesep 'VFN-Simulations' filesep 'VFNlib']));

close all
clear all

global itr pk_min_hists FMTO


%% Zaber Setup
fprintf('---Performing Zaber Setup\n')

%~~ Connect to Zabers to query position
VFN_setUpBenchZabers; % Instantiate and rename the zabers

%- Make sure pmX is out of the way
if isfield(Zabs, 'pmX')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
end

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

%% Red Thorlabs PM setup
fprintf('---Performing redPM Setup (if needed)\n')
VFN_setUpPM100D

pmCalWvl = 635;
pmNread  = 100;

% Set operating wavelength for calibration
VFN_PM100D_setWavelength(PM, pmCalWvl);

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
TolX = .001;                     %Smallest step size

%-- Maximum number of function evals in minimizer
nFunEvals = 150;

%-- Set Options
%opts = optimoptions('patternsearch','StepTolerance',TolX,'FunctionTolerance', TolFun, 'MaxTime',60*15,'MaxFunctionEvaluations',nFunEvals,'Display','iter','PlotFcn','psplotbestf');
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
fibX0 = 7.374965;% MMF: 7.408107;%    % Fiber X position [in mm]
fibY0 = 7.096703;% MMF: -9.862401;%   % Fiber Y position [in mm]
fibZ0 = 9.508524;% MMF: 9.209940;%    % Fiber Z position [in mm]

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

%opts = optimoptions('patternsearch','StepTolerance',TolX,'FunctionTolerance', TolFun, 'MaxFunctionEvaluations',nFunEvals,'Display','iter','PlotFcn','psplotbestf');%'MaxTime',60*15,
%[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
opts = optimset('Display', 'iter', 'TolX', TolX, 'TolFun', TolFun, 'MaxFunEvals', nFunEvals, 'PlotFcns', @optimplotfval);
[X, fval, exitflag, output] = fminsearch(func, X0, opts);

%% Read norm value on redPM
%-- Move redPM into the beam
%Move fiber stage out of the way to avoid collision with pmX
VFN_PIStage_move(PIdevs.fibZ, PIdevs.fibZ_lower);
fprintf('  fibZ at %0.2f',VFN_PIStage_getPos(PIdevs.fibZ))

% Brief pause just to be safe
pause(2);

%Move if no collision
fprintf('\nMoving pmX...')
VFN_Zab_move(Zabs.pmX, Zabs.pmX_upper);
fprintf(' Doing measurments...')

%-- Read from red PM 
% Wait for PM value to settle
pause(2);

% Take measurements of the power (use try catch to avoid errors)
NTries = 5;
for i = 1:NTries
    try 
        % This is where the read often fails 
        pmRead = VFN_PM100D_read(PM, pmNread);
    catch ME
        % An error occurred, let the user know
        warning('PM read failed')
        if i < 5
            % We have tried less than NTries times, try again

            % Wait a moment to give the redPM time to figure itself out
            pause(1);

            % First, try playing with the wavelength settings
                % This seems to help for some reason
            VFN_PM100D_setWavelength(PM, pmCalWvl);
            pause(0.1);
            VFN_PM100D_getWavelength(PM);

            % Now try again (continue will jump straight to next iter.)
            fprintf('Trying again\n')
            continue
        else
            % We've tried 5 times, it's time to give up
            error('Multiple redPM reads failed')
        end
    end
    % If we get here, succeeded in reading so exit for-loop and proceed
    break
end

% Move redPM zaber, pmX, out of beam again
VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
fprintf('\nDone. pmX moved to %0.2f',VFN_Zab_getPos(Zabs.pmX))

%-- Get norm value
FL = 0.9654;  % (1-loss) Fraction per surface
PL = 0.9966;  % (1-loss) Fraction per meter of fiber

%-- Get FMTO sensitivity at given wavelength
%fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
fmtoSnsScl = 0.62; % value from manual comparison

%-- When actual PM reading was done, normalize as usual
% Account for loss from fiber lens
    % 3.4% reflection spec at 780nm, 0.3% reflection measured at 635nm
%pmRead1 = mean(pmRead)*(1-0.034);%0.997;    % Updated on 2/8/19 based on 10/26/18
% Convert the calibration PM value from W to uW
pmRead1 = mean(pmRead)*10^6;
% Translate the calibration PM value from uW to V (at 10^6)
pmRead1 = pmRead1*fmtoSnsScl;     % V/uW conversion is ~ 0.63 @633; ~0.9 @ 780
% Account for Fresnel and Propagation losses
nrmVal = (1/pmRead1)*(1/(FL^2))*(1/(PL^2));

%% Display results

%%--- Print results
fprintf('\n      - Peak found =  %f', -1*fval);
fprintf('\n      -   (normed) =  %f', -1*fval*nrmVal);
fprintf('\n      - Ideal pos  =  (%f mm, %f mm, %f mm)\n', X(1), X(2), X(3));

%-- Display results
figure; plot(pk_min_hists.X(:,1)); title('X(1)');
figure; plot(pk_min_hists.X(:,2)); title('X(2)');
figure; plot(pk_min_hists.X(:,3)); title('X(3)');
figure; plot(pk_min_hists.PWR(:)*nrmVal*100); title('PWR (Normed)'); ylabel('Normed Coup [%]');
figure; semilogy(pk_min_hists.PWR(:)*nrmVal*100); title('LOG PWR (Normed)');ylabel('logNormed Coup [%]');


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

%-- Move Calibration PM out of the beam (again; just to be safe)
fprintf('PM Zab pos: %f\n', VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower));

%-- Clean up Zabers
VFN_cleanUpZabers;

%-- Clean up redPM
VFN_cleanUpPM100D;