function VFN_KLS_disableOut(KLS)
% VFN_KLS_disableOut Function to set the laser emission off
%
%   EXAMPLE:______
%   VFN_KLS_disableOut(KLS)
%       KLS:    the python KLS object returned by VFN_setUpKLS
%       

KLS.disableOut();

end