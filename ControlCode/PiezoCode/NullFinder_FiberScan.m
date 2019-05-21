%% add the necessary functions
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

close all
clear all

global s min_hist

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

%% Define Search Parameters
%-- Set initial params: 
fibX0 = 75;             % Fiber X position
fibY0 = 75;             % Fiber Y position
    %NOTE!!: The zaber functions use mm but I rescale here so that the a scalar 
    % step tolerance can be used. (We want ~0.3V step tol for the piezos and ~
    % 3 microns for zabers. By scaling x100, 3 microns = 0.3 and thus we can
    % sharing same step tolerance for both actuator types.
nrmFac = 100;
focZ0 = 3.64450*nrmFac;      % Fiber Z position [in microns]

%-- Vectorize initial params
X0 = [fibX0, fibY0, focZ0];

%-- Set Bounds and Constraints
% Bounds
fibBnd = 15;                                %Max distance from fibX0/Y0
focBnd = 2;                                %Max distance from focZ0
LB = [fibX0-fibBnd;                          %Lower bounds
      fibY0-fibBnd; 
      focZ0-focBnd];
UB = [fibX0+fibBnd;                          %Upper bounds
      fibY0+fibBnd;
      focZ0+focBnd];
% Constraints
A   = [];
b   = [];
Aeq = [];
beq = [];
% Set nonlcon blank so that we can provide 'options'
nonlcon = [];


%-- Set tolerance and step size values
sTol = .3;                     %Smallest step size

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

%% Perform minimization
fprintf('\n-- Running Minimizer\n')
%-- Define history vector to contain iteration information
min_hist.X      = [];
min_hist.PWR    = [];

%-- Define iteration function
    % This moves the actuators to the new position and measures the donut power
func = @(X) MinFunc_FiberScan(X, fibZ, MDT, FMTO_scale, nrmFac);

[X, fval, exitflag, ouput] = patternsearch(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
%[X, fval, exitflag, ouput] = fmincon(func, X0, A, b, Aeq, beq, LB, UB, nonlcon, opts);
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

%% Print results
fprintf('\n      - Null found =  %f', fval);
fprintf('\n      - Ideal pos  =  (%5.2f V, %5.2f, %f mm)\n', X(1), X(2), X(3)/nrmFac);
