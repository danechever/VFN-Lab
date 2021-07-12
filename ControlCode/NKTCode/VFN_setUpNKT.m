%Set up the NKT super K and varia. Updates the NKT sub-struct
% which contains pertient information.
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

%-- Check if NKT exists, has the CONNECTED field, and is not connected
    % This takes advantage of the MATLAB logical short-circuiting 
if ~exist('NKT', 'var') || ~isfield(NKT, 'CONNECTED') || ~NKT.CONNECTED
    
    %NKT_lib_PATH = 'C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode\NKTCode\';
    NKT_lib_PATH = '/home/vfndev/Documents/MATLAB/VFN-Lab/ControlCode/NKTCode/';
    if count(py.sys.path, NKT_lib_PATH) == 0
        insert(py.sys.path, int32(0), NKT_lib_PATH);
    end

    NKT.numTries = 10;
    NKT.delay = 1;
    
    disp('*** Connecting to NKT ... ***');
    
    % Create Nkt object using nkt_mod.py
      %NOTE: if one doesnt work, try the other
    %NKT.nktobj = py.nkt_mod_falco_py3.Nkt('/dev/ttyUSB5',115200);%115200%
     NKT.nktobj = py.nkt_mod_falco.Nkt('/dev/ttyUSB5',115200);%115200%
    
    NKT.CONNECTED = true;
    
    disp('*** NKT connected. ***');
else
    
    disp('*** NKT already connected. ***');
    
end