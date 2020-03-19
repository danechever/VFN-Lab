% Load dotNET assembly
asm = NET.addAssembly('C:\Program Files\Thorlabs\Thorlabs OSA\ThorlabsOSAWrapper.dll');
import ThorlabsOSAWrapper.*

% Define common interface as dedicated object
%lib = ThorlabsOSAWrapper.LibDeviceInterface;

% Find/Initialize OSA on computer
devs = DeviceLocator.InitializeSpectrometers();

% Connect to first device found
osa = LibDeviceInterface(0)
disp('Connected to OSA')

% Set measurement parameters
disp('Setting Sensitivity Mode')
osa.SetSensitivityMode(1)               % 1 = low sens
disp('Setting Resolution Mode')
osa.SetResolutionMode(1)                % 1 = low res
disp('Setting Setting Auto Gain')
osa.AcquisitionSettings.AutomaticGain = true;    % Set to auto gain
disp(['Is auto gain on?: ', mat2str(osa.AcquisitionSettings.AutomaticGain)]);

% Do a single spectrum scan
disp('Acquiring Spectrum.')
osa.AcquireSingleSpectrum()
disp('Specturm Acquired')

% Wait for spectrum to finish, just in case
pause(2)

% Close connection
osa.CloseSpectrometer()
disp('Connection to OSA closed')

% Check that spectrum is not blank
fprintf('\nLength of measuered data: %i', osa.GetLastSpectrumLength());
fprintf('\n  - If non-zero, data is probably present\n')

% Create spectrum object to hold read data
%%% THIS IS WHERE THE FAILURE OCCURRS %%%
disp('Attempt to create SpectrumStruct ')
spec = SpectrumStruct(osa.GetLastSpectrumLength());

% The following are things I've tried that did not work
disp('Try Other stuff')
osa.GetLastSpectrum()   % call without provided struct
osa.GetLastSpectrum(ThorlabsOSAWrapper.SpectrumStruct(osa.GetLastSpectrumLength())) % creat struct in-line

% It seems 'spectrum_t*' is defined in FTSlib.dll. As such, try adding this dll as well
asm2 = NET.addAssembly('C\Program Files\Thorlabs\Thorlabs OSA\FTSLib.dll');