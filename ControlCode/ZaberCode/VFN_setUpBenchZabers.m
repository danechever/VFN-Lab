%Script to prepare the Zabers specifically for the VFN bench
%
% This will rename the various axes to match the nomencalture used in 
% our main control scripts. 
% --> Instantiates the axes with the appropriate names in the workspace.
%
% Any unrecognized axes (not in "names" below) are left unchanged

%% Define name list
% Populate a struct that translates an old name (w/ SN) to the simple name
names = struct('AX_52714', 'vortX', ...
               'AX_51463', 'vortY', ...
               'AX_52893', 'fibZ', ...
               'AX_59214', 'pmX'); 


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

clear('axs', 'names')