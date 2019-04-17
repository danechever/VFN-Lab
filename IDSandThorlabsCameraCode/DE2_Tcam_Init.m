%{
Script to initialize Thorlabs (DCC 1545m) camera.
- Creates instance of camera as well as image data structure for accessing
    image data later.
- Sets the Exposure time per user input as well as other camera settings
    Pixel Clock Range: 5 - 43 MHz (IDS Cam: 5 - 36MHz)
    Exposure Time Range: .037 - 983ms (IDS Cam: .340 - 14582 ms)

** DOES NOT CLOSE CAMERA INSTANCE. This script is meant to initialize the camera
    for use in other code. Therefore, it creates a camera instance which
    must then be closed to free the camera for later use. This must be done
    by CALLING DE2_Tcam_Close at the end of the main function. FAILURE TO
    DO SO WILL LEAD TO CAMERA AND CODE PROBLEMS.

- Useful PDF for understanding function calls within .NET profile
    https://www.thorlabs.com/drawings/46f3af850ccfa06a-9BB3D44A-0ED1-7FB0-3FB956E232A2860C/DCC1545M-.NETProgrammingInterface.pdf
******************************************************
- Arguments:
    expTime         exposure time for captures, in ms

- Returns:
    cam             Data structure containing camera instance; must be
                        passed between uc480_ functions
    img             Data structure containing image data; needed for
                        accessing data from RAM and buffer    
                
- Dependencies:
    uc480DotNet.dll .NET library for camera. Available from SDK.
******************************************************

Initial Version:    Daniel Echeverri, 10/24/2017
Last Edit:          Daniel Echeverri
Last Modified:      11/05/2017
%}
function [cam, img] = DE2_Tcam_Init(expTime)
    %Check if assembly is present. Import Otherwise
    %May need to change search location if library is moved
    asm = System.AppDomain.CurrentDomain.GetAssemblies;
    if ~any(arrayfun(@(n) strncmpi(char(asm.Get(n-1).FullName), ...
            'uc480DotNet', length('uc480DotNet')), 1:asm.Length))
        NET.addAssembly(...
            'C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
    end

    %% Create and initialize camera instance
    %Create CCD object handle
    cam = uc480.Camera();

    %Open and connect to first camera: 0 = 1st cam
    char(cam.Init(0));

    %Set display mode to bitmap (DiB); Captured to RAM
    if ~strcmp(char(cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB)), ...
            'SUCCESS')
        error('Could not set display mode');
    end

    %Set color mode to 8-bit RAW
    if ~strcmp(char(cam.PixelFormat.Set(uc480.Defines.ColorMode.SensorRaw8)), ...
            'SUCCESS')
        error('Could not set pixel format');
    end

    %Set trigger mode to software (single image acquisition when Freeze() is called)
    if ~strcmp(char(cam.Trigger.Set(uc480.Defines.TriggerMode.Software)), ...
            'SUCCESS')
        error('Could not set trigger format');
    end
    %Set pixel clock which manages fundamental timing unit of camera. 
    cam.Timing.PixelClock.Set(30);              %MHz
    cam.Timing.Exposure.Set(expTime);

    %Allocate image memory, define id as img.ID
    [ErrChk, img.ID] = cam.Memory.Allocate(true);
    if ~strcmp(char(ErrChk), 'SUCCESS')
        error('Could not allocate memory');
    end

    %% Take single image to instantiate img.() structure
    %Obtain image information
    [ErrChk, img.Width, img.Height, img.Bits, img.Pitch] ...
        = cam.Memory.Inquire(img.ID);
    if ~strcmp(char(ErrChk), 'SUCCESS')
        error('Could not get image information');
    end

    %Acquire image
    if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
        error('Could not acquire image');
    end

    %Extract image, store in tmp as vector array
    [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    if ~strcmp(char(ErrChk), 'SUCCESS')
        error('Could not obtain image data');
    end

    %Reshape image into WxH matrix; convert tmp to uint8 vector
    img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
end