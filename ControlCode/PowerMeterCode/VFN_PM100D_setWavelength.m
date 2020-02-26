function res = VFN_PM100D_setWavelength(PM, lam)
% VFN_PM100D_SETWAVELENGTH Function to set the operating wavelength of the PM100D
%
%   - The PM100D should have already been set up using VFN_setUpPM100D
%
%   Arguments: 
%   - PM:       struc containing objects related to PM100D
%                 Created by VFN_setUpPM100D
%   - lam:    	wavelength to set for PM100D
%
%   Returns: 
%   - res:      Resulting wavelength as read from device
%
%   EXAMPLE:______
%   res = VFN_PM100D_SETWAVELENGTH(PM, 775)
%       Sets the operating wavelength to 775nm
%       * res: (775) resulting wavelength
%
%   See also VFN_setUpPM100D, VFN_PM100D_read

if ~isa(PM.PM, 'py.ThorlabsPM100.ThorlabsPM100.ThorlabsPM100')
    error('PM.PM must be an instance of the ThorlabsPM100 class')
end

PM.PM.sense.correction.wavelength = lam;

res = PM.PM.sense.correction.wavelength;

end