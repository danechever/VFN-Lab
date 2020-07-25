function read = VFN_PM100D_read(PM, varargin)
% VFN_PM100D_READ Function to read the power from the PM100D
%
%   - The PM100D should have already been set up using VFN_setUpPM100D
%
%   - Optionally takes in number of consecutive reads to perform
%
%   Arguments: 
%   - PM:       struc containing objects related to PM100D
%                 Created by VFN_setUpPM100D
%   - nRead:    (optional) number of consecutive reads to perform
%
%   Returns: 
%   - read:     values read from PM. If nRead>1, this will be a vector
%
%   EXAMPLE:______
%   read = VFN_PM100D_READ(PM)
%       Reads a single power value from the PM100D connected to PM 
%       * read:   (scalar)power value read from PM100D
% 
%   read = VFN_PM100D_READ(PM, 100)
%       Reads 100 power values from the PM100D connected to PM 
%       * read:   (vector)power values read from PM100D
%
%   See also VFN_setUpPM100D, VFN_PM100D_setWavelength

if ~isa(PM.PM, 'py.ThorlabsPM100.ThorlabsPM100.ThorlabsPM100')
    error('PM.PM must be an instance of the ThorlabsPM100 class')
end

if nargin > 1
    % Argument provided; perform number of requested reads
    nRead = varargin{1};
    % Preallocate vector
    read = nan(nRead,1);
    for i = 1:nRead
        read(i) = PM.PM.read;
    end
else
    % No argument provided; perform single read
    read = PM.PM.read;
end

end