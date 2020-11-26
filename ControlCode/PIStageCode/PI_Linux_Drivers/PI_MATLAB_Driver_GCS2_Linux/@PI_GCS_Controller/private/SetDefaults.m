function c = SetDefaults(c)

% ID to use C API
c.ID = -1;

% ID to use C API with DaisyChain
c.DC_ID = -1;

% This is needed to retrieve the parameters for all axes without specifing them (create array for C-interface)
c.NumberOfAxes = 0;     

% Used for daisy chain
c.ConnectedDaisyChainDevices = [];

% load available dll functions
c.dllfunctions = libfunctions ( c.libalias );