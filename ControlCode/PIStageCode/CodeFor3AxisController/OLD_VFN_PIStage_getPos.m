function currentPos = VFN_PIStage_getPos(axis, PIdevice)
% VFN_PIStage_getPos Function for quering the position of the PI stages
%   
%   EXAMPLE:______
%   currentPos = VFN_PIStage_getPos(axis, PIdevice)
%         1) 'axis'         = multidimensional cell containing the axis names for use by other scripts
%         2) 'PIdevice'     = instance of device (USB connection) object

%% Get Positions axis by axis
currentPos = zeros(1, length(axis));

for axisID=1:length(axis)
    try
        currentPos = PIdevice.qPOS(char(axis(axisID))); % Returns mm already
    catch exception
        % Close port if a MATLAB error occurs, otherwise it remains locked
        VFN_PIStage_cleanUp;
        %fclose(axis.Protocol.Port); %From Zab_getpos
        rethrow(exception);
    end
end

clear axisID
end