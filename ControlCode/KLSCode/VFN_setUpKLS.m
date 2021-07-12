%% Set python path correctly
%-- Define the paths
% Get the current matlab path (since that'll put us relative to the calling
% script
matpath = pwd();
% Define the KLS path relative to the matlab path
klspathdir = 'KLSCode';
klspathfull = [matpath, filesep, klspathdir];

%-- Set the path
% Get current python path
pypath = py.sys.path;
% Iterate over elements of the path list to see if desired path is present
klsFlag = false;
for ind = 1:length(pypath)
    if strcmp(string(pypath{ind}), klspathfull)
        klsFlag = true;
    end
end

% Add KLS code directory to python path if needed
if ~klsFlag 
    insert(py.sys.path, int32(length(py.sys.path)), klspathfull)
end

%% Connect to the laser
% NOTE: set the path here based on the actual connected path
KLS = py.TLS.TLS_Device('/dev/ttyUSB6');

% Open the port
KLS.open()

% Check communications
disp('Laser settings')
KLS.printDev()

%% Cleanup after ourselves
clear matpath klspathdir klspathfull pypath klsFlag