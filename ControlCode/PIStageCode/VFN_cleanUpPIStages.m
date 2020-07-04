% This will close the connection to the E-873 controllers
%
% It assumes the following variables were already created by _setUpPIStages
%   'Controller' = instance of controller object
%   'PIdevs'     = struct with instances of device (USB Connection) objects
%
% NOTE: This code no longer sets the controllers to open loop
    
    
%% Close connection to all devices in the struct
axs = fieldnames(PIdevs);
for i = 1:numel(axs)
    PIdevs.(axs{i}).CloseConnection();
end

%% Unnload the dll and "destroy" the Controller object
Controller.Destroy ();
clear('axs', 'Controller', 'PIdevs', 'i')