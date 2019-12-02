%% Comment:
% This will close the connection to the E-873 controller
%
% NOTE: This will turn off the servo, thereby setting the axis to open-loop
%
% It assumes the following variables have already been created by _init.m
    % 1) 'axis'         = multidimensional cell containing the axis names for use by other scripts
    % 2) 'Controller'   = instance of controller object
    % 3) 'PIdevice'     = instance of device (USB connection) object
    
%% Open servo
PIdevice.SVO ( axis, zeros(1, length(axis)) );
    
%% If you want to close the connection
PIdevice.CloseConnection ();

%% If you want to unload the dll and destroy the class object
Controller.Destroy ();
clear Controller PIdevice axis