%{
Script to run the Thorlabs Camera in live view mode
- V2: implements object abstraction to simplify commands
******************************************************
- Arguments:
    NONE            USER CAN ADJUST EXPOSURE TIME on CAM

- Returns:
    NONE            
                
- Dependencies:
    uc480DotNet.dll .NET library for camera. Available from SDK.
    DE2_Tcam_Init   Function to create and set camera object
    DE2_Tcam_Frame  Function to capture single camera frame
    DE2_Tcam_Close  Function to close and release camera
******************************************************

Compiled By:    Daniel Echeverri and JR Delorme
Last Modified:  10/24/2017
%}
clear all

%% Programmable Parameters
%****************************************************
% Define your variables and other settings for the analysis
%****************************************************

%Variables
nFrames     = 25    ;   % Number of frames to capture.
frameRate   = 5    ;   % Frames per second, limited by exposure time
expTime     = 0.1   ;   % Exopsure time
xCrop       = 450:550;  % Columns to include when cropped (isCrop)
yCrop       = 500:600;  % Rows to include when cropped (isCrop)

%Save folder (Fl) and name (Nm)
fileFl      = 'C:\Users\Daniel\Desktop\Mawet201718';
fileNm      = 'Test';   %Image number automatically appended to end of name

%Settings
isLog       = 0     ;   % 0 for linear scale, 1 for log
isCrop      = 1     ;   % 0 for full, 1 for cropped
                        %       Gets cropped to (xCrop, yCrop)
isDispFull  = 1     ;   % 0 for just cropped or full, 1 for both cropped and full
                        %       0 displays image defined by isCrop
                        %       If isCrop == 0, only Full image will
                        %       display regardless of isDispFull to avoid
                        %       plotting two Full images
isSaveData  = 0     ;   % 0 for no save, 1 to save into fileFl as fileNm
                        %       Saves the original image, regardless of
                        %       isLog or isCrop.
                        %       Will save as many frames as nFrames

                        
%% Camera Intialization
[cam, img] = DE2_Tcam_Init(expTime);

%% Main Loop
%****************************************************
% Iterates through image capture to emulate live video
%****************************************************
for n = 1:nFrames
    %Capture frame
    Im = DE2_Tcam_Frame(cam, img);
    
    
    dispCrop = Im;        %dummy variable for crop manipulation
    dispFull = Im;        %dummy variable for Full manipulation
    %Logic to change scale
    if isLog
        dispCrop = log(dispCrop+1);
        dispFull = log(dispFull+1);
    end
    %Logic to crop image
    if isCrop
        dispCrop = dispCrop(yCrop, xCrop);
    end
    
    %Display first plot, disp
    figure(1)
      imagesc(dispCrop);
      axis image
  
    if (isDispFull && isCrop)  
      figure(2)  
        imagesc(dispFull)
        colorbar
        axis image 
        caxis([0.1 5])
    end
  
    if isSaveData
        fnm = [fileFl '\' fileNm '_' num2str(i,'%03.0f')];
        fitswrite(Im, fnm);
    end
    max(Im(:))
    pause(1/frameRate)
end

%% Close Camera: 
%THIS SHOULD RUN AT THE END EVERY TIME. IF CODE FAILS FOR WHATEVER REASON, 
%  RUN THIS CHUNK TO SUCCESSFULLY SHUT DOWN AND RESET CAMERA. 
DE2_Tcam_Close(cam);
