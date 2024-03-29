%Script to prepare the PI Stages specifically for the VFN bench
%
% --> Instantiates the stages
% --> Turns on the servos (only on stages in "names" below)
% --> Homes if needed (only for stages in "names" below)
% --> Renames the stages in the workspace
% --> Adds axis limits defined by user to PIdevs struct
%
% Any unrecognized stages (not in "names" below) are left unchanged

%% Define name list
% Populate a struct that translates an old name (w/ SN) to the simple name
names = struct('AX_119063717', 'fibX', ...
               'AX_119063721', 'fibY', ...
               'AX_119069544', 'fibZ');
           
%% Define Axis limits
% Add limits as struct fields. 1st element is struct name, 2nd is limit
    % Name should have the same prefix as above with a descriptive suffix
    % NOTE: current alignment places no hard limits on stages
    % Just need to make sure pmX is out of way. 
    % Still, avoid far edges since they can be unstable in closed loop
% Limits are inclusive (ie. stage can move up to but not beyond the limit)

% LIMITS FROM TRANSMISSIVE BENCH:
%limits = struct('fibX_lower', -10.0000, ...
%                'fibX_upper',  12.0000, ...
%                'fibY_lower', -12.0000, ...
%                'fibY_upper',  10.0000, ...
%                'fibZ_lower', -10.0000, ... % fibZ moves here to avoid pmX
%                'fibZ_upper',  12.0000);

% LIMITS FROM PoRT:
limits = struct('fibX_lower', -12.0000, ...
                'fibX_upper',  10.0000, ... % Could be set to 12.0 but gets very close to iris in some Z positions
                'fibY_lower', -12.0000, ...
                'fibY_upper',  12.0000, ...
                'fibZ_lower', -12.0000, ... % fibZ moves here to avoid pmX
                'fibZ_upper',  12.0000);


%% Connect to PI stages
VFN_setUpPIStages; 
           
%% Rename, turn on servo, and home the given stages
% Get all elements in the struct so that we can iterate through them
axs = fieldnames(names);

% Iterate through stages
for i = 1:numel(axs)
    ax = axs{i};
    % Check if desired stage is present in connected PI stages
    if isfield(PIdevs, ax)
        % Turn on the servo
        VFN_PIStage_setServo(PIdevs.(ax), 1);
        
        % Home the stage if needed
        if ~VFN_PIStage_isHomed(PIdevs.(ax))
            % Perform reference move
            fprintf('\n%s (%s) needs to be homed\n', ax, names.(ax))
            VFN_PIStage_home(PIdevs.(ax));
        end
        
        % Rename the axis
        PIdevs.(names.(ax)) = PIdevs.(ax);
        % Remove old field name
        PIdevs = rmfield(PIdevs, ax);
    end
end

%% Add axis limits if the given axis exists
% Get all elements of limits struct to iterate through them
lms = fieldnames(limits);

% Get all elements of PIdevs struct for reference
axs = fieldnames(PIdevs);

% Iterate through limits
for i = 1:numel(lms)
    lm = lms{i};
    
    if contains(lm,axs)
        % An axis present in PIdevs matches the prefix for this limit
            % Thus, add the limit to the PIdevs struct
        PIdevs.(lm) = limits.(lm);
    end
end    

clear('axs', 'ax', 'names', 'lms', 'lm', 'limits', 'i')