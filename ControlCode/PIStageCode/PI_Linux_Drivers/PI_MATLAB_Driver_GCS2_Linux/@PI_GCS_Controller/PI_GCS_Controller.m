function c = PI_GCS_Controller()
% PI_GCS_Controller controller class constructor
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.


if nargin==0
    c.libalias = 'PI';
    
    if ( libisloaded ( c.libalias ) )
        disp('PI_MATLAB_Driver_GCS2 already loaded.');
    else
        disp ( 'Loading PI_MATLAB_Driver_GCS2 ...' );
        LoadGCSDLL ( c );
        disp('PI_MATLAB_Driver_GCS2 loaded successfully.');
    end
        
    c = SetDefaults(c);
    c = class(c,'PI_GCS_Controller');
end