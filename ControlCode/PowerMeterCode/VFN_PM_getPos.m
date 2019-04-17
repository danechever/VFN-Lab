function resPos = VFN_PM_getPos(pmObj, pm_scale)
% VFN_PM_resPos Function for quering the power on the femto power meter
%   
%   EXAMPLE:______
%   resPos = VFN_PM_getPos(pmObj)
%       pmObj:      Instance of power meter object ('PM' by default)
%       pm_scale:   Scaling value for power meter. 
%           *** THIS IS THE VALUE THAT KNOB ON THE POWER METER IS POINTING
%               TO. FOR EXAMPLE, '10^8'.
%       resPos:     resulting power (in nW)

if ~isa(pmObj, )
    error('The first argument must be a daq object')
end


%% get power level
resPos = pmObj.inputSingleScan()/pm_scale;
end