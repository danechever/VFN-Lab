%{
Script to modify exposure time of Thorlabs (DCC 1545m) camera.
- Modifies camera exposure time based on input argument
- Attempts to modify pixelClock accordingly to ensure appropriate clocking
- Pixel Clock Range: 5 to 43 MHz (IDS Cam: 5 - 36MHz)
       Pixel clock defines the rate at which pixels are read. It is the
       fundamental unit for camera timing. Higher rates means shorter
       exposure times are possible but also higher camera temperatures.
- Exposure time range: .037ms to 983ms  (IDS Cam: .340 - 14582 ms)
       The exposure time is limited by the framerate and pixel clock.
       Minimum exposure time requires maximum clock rate and vice versa.
- Could expirement with other functions to ensure that exposure time value
        is valid and to optimize pixelClock relative to exposure time.
** DOES NOT CREATE CAMERA INSTANCE. This script is meant to change exposure
    time of an existing camera object. CALL DE2_Tcam_Init to create
    instance.
** DOES NOT CLOSE CAMERA INSTANCE. Camera must be closed by CALLING DE2_Tcam_Close 
    at the end of the main function. FAILURE TO DO SO WILL LEAD TO CAMERA 
    AND CODE PROBLEMS.
******************************************************
- Arguments:
    cam             Data structure containing camera instance; must be
                        passed between uc480_ functions
    expTime         Exposure time value in ms. Range: .037-983 ms

- Returns:
    cam             Modified cam structure containing updated exposure time
                
- Dependencies:
    uc480DotNet.dll .NET library for camera. Available from SDK.
******************************************************
Initial Version:    Daniel Echeverri, 10/24/2017
Last Edit:          Daniel Echeverri
Last Modified:      10/24/2017   
%}

function [cam] = DE2_Tcam_ExpTime(cam, expTime)
    cam.Timing.Exposure.Set(expTime);    
end