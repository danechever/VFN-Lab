%Script to prepare the Zabers specifically for the VFN bench
%
% This will rename the various axes to match the nomencalture used in 
% our main control scripts. It will also add upper and lower limits to the
% struct. 
% --> Instantiates the axes with the appropriate names in the workspace.
% --> Adds limit values to struct. These limits are NOT automatically
%       handled. User can reference them in their code though.
%
% Any unrecognized axes (not in "names" below) are left unchanged

%% Define name list
% Populate a struct that translates an old name (w/ SN) to the simple name
names = struct('AX_52714', 'vortX', ...
               'AX_51463', 'vortY', ...
               'AX_59214', 'pmX'); 
               %'AX_52893', 'fibZ', ...

%% Define Axis limits
% Add limits as struct fields. 1st element is struct name, 2nd is limit
    % Name should have the same prefix as above with a descriptive suffix
limits = struct('vortX_lower', 0, ...
                'vortX_upper', 25, ...
                'vortY_lower', 0, ...
                'vortY_upper', 25, ...
                'pmX_lower', 0, ...
                'pmX_upper', 25);

%% Connect to zabers
VFN_setUpZabers; 
           
%% Rename struct fields to match script nomenclature

% Get all elements in the struct so that we can iterate through them
axs = fieldnames(names);

% Iterate through axes
for i = 1:numel(axs)
    % Check if desired axis to rename is present in connected Zabers
    if isfield(Zabs, axs{i})
        % Rename the axis
        Zabs.(names.(axs{i})) = Zabs.(axs{i});
        % Remove old field name
        Zabs = rmfield(Zabs, axs{i});
    end
end

%% Add axis limits if the given axis exists
% Get all elements of limits struct to iterate through them
lms = fieldnames(limits);

% Get all elements of PIdevs struct for reference
axs = fieldnames(Zabs);

% Iterate through limits
for i = 1:numel(lms)
    lm = lms{i};
    
    if contains(lm,axs)
        % An axis present in PIdevs matches the prefix for this limit
            % Thus, add the limit to the PIdevs struct
        Zabs.(lm) = limits.(lm);
    end
end    

clear('axs', 'names', 'limits', 'lms', 'i', 'lm')
