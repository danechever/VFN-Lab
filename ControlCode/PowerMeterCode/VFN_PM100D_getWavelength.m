function lam = VFN_PM100D_getWavelength(PM)
% VFN_PM100D_GETWAVELENGTH Function to get the operating wavelength of the PM100D
%
%   - The PM100D should have already been set up using VFN_setUpPM100D
%
%   Arguments: 
%   - PM:       struc containing objects related to PM100D
%                 Created by VFN_setUpPM100D
%
%   Returns: 
%   - lam:      Current wavelength as read from device
%
%   EXAMPLE:______
%   lam = VFN_PM100D_GETWAVELENGTH(PM)
%       Returns the current operating wavelength of the device
%       * lam: (775) current wavelength
%
%   See also VFN_setUpPM100D, VFN_PM100D_read

if ~isa(PM.PM, 'py.ThorlabsPM100.ThorlabsPM100.ThorlabsPM100')
    error('PM.PM must be an instance of the ThorlabsPM100 class')
end

lam = PM.PM.sense.correction.wavelength;

end