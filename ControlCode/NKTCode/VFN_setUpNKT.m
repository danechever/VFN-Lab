%Set up the NKT super K and varia. Updates the NKT sub-struct
% which contains pertient information.
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

if(~NKT.CONNECTED)
    
    NKT_lib_PATH = 'C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode\NKTCode\';
    if count(py.sys.path, NKT_lib_PATH) == 0
        insert(py.sys.path, int32(0), NKT_lib_PATH);
    end

    NKT.numTries = 10;
    NKT.delay = 1;
    
    disp('*** Connecting to NKT ... ***');
    
    % Create Nkt object using nkt_mod.py
    NKT.nktobj = py.nkt_mod_falco.Nkt('/dev/ttyNKT',115200);%115200%
    
    NKT.CONNECTED = true;
    
    disp('*** NKT connected. ***');
else
    
    disp('*** NKT already connected. ***');
    
end