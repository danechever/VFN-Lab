% Clearn up the NKT superK and varia. Updates the NKT sub-struct.
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019


% Close the Nkt object
NKT.nktobj.close();

% remove the connection object from the bench.NKT struct
NKT = rmfield(NKT,'nktobj');

NKT.CONNECTED = false;

disp('*** NKT disconnected. ***');