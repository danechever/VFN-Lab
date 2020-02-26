% Close the connection to the Thorlabs MDT
%
%   - Also deletes MDT variable from workspace
%
%   See also VFN_setUpMDT, VFN_MDT_move


% Close connection to the device
fclose(MDT);

% Delete MDT variable
clear('MDT')