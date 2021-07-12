function VFN_KLS_setPower(KLS, power)
% VFN_KLS_setPower Function to set the KLS power:
%
%   EXAMPLE:______
%   VFN_KLS_setPower(KLS, 0.5)
%       KLS:   the python KLS object returned by VFN_setUpKLS
%       power: Power to which the laser should be set (fractional)
%       

KLS.setPowerPrct(power);

end