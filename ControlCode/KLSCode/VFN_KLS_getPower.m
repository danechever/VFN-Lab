function [actpwr, setpwr] = VFN_KLS_getPower(KLS)
% VFN_KLS_getPower Function to get the KLS power:
%
%   EXAMPLE:______
%   [actpwr, setpwr] = VFN_KLS_getPower(KLS)
%       KLS:    the python KLS object returned by VFN_setUpKLS
%       actpwr: The actual power state, returned by the KLS laser
%       setpwr: The power setpoint (goal value set by user's previos comm)
%       

actpwr = KLS.reqPowerAct();
setpwr = KLS.reqPowerSet();

end