%{
Script to close Thorlabs (DCC 1545m) camera instance.
- Closes camera instance to free camera for other use.
** MUST BE CALLED AT THE END OF MAIN SCRIPT TO ENSURE CAMERA IS RELEASED.
    Failure to do so will lead to camera and code problems.
** Calling this function does not cause code or communication errors. As
    such if there is any concern with the status of a camera instance, call
    this function to ensure instance is closed. At worst, an error will be
    thrown if camera instance did not exist.
******************************************************
- Arguments:
    cam             Data structure containing camera instance; must be
                        passed between uc480_ functions

- Returns:
    NONE
                
- Dependencies:
    uc480DotNet.dll .NET library for camera. Available from SDK.
******************************************************

Initial Version:    Daniel Echeverri, 10/24/2017
Last Edit:          Daniel Echeverri
Last Modified:      10/24/2017
%}
function DE2_Tcam_Close(cam)
if ~strcmp(char(cam.Exit), 'SUCCESS')
    error('Could not close camera');
end