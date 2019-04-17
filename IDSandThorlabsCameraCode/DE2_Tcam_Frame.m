%{
Script to capture single image data from Thorlabs (DCC 1545m) camera.
- Captures and converts image data from the camera into useable image
    matrix. Returns this matrix.
** DOES NOT CREATE CAMERA INSTANCE. This script is meant to take a single
    image on an existing camera object. CALL DE2_Tcam_Init to create
    instance.
** DOES NOT CLOSE CAMERA INSTANCE. Camera must be closed by CALLING DE2_Tcam_Close 
    at the end of the main function. FAILURE TO DO SO WILL LEAD TO CAMERA 
    AND CODE PROBLEMS.
******************************************************
- Arguments:
    cam             Data structure containing camera instance; must be
                        passed between uc480_ functions
    img             Data structure containing image data; needed for
                        accessing data from RAM and buffer

- Returns:
    Im              img.Height x img.Width matrix containing image data
                        that is directly usable within matlab image scripts
    img             Updated image data
                
- Dependencies:
    uc480DotNet.dll .NET library for camera. Available from SDK.
******************************************************

Initial Version:    Daniel Echeverri, 10/24/2017
Last Edit:          Daniel Echeverri
Last Modified:      10/24/2017
%}
function [Im, img] = DE2_Tcam_Frame(cam, img)
    %Acquire image
    if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
        error('Could not acquire image');
    end
    %Extract image
    [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    if ~strcmp(char(ErrChk), 'SUCCESS')
        error('Could not obtain image data');
    end
    %Reshape image
    img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);

    %Provide image data as usable double matrix
    Im = double(img.Data);
end