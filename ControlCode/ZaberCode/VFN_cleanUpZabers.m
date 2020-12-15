%% Close the zaber port object and delete the Zabs struct

%-- Close the port if it exists (repeated .close() calls cause no issues)
if isfield(Zabs, 'port')
    Zabs.port.close();
    Zabs = rmfield(Zabs, 'port');
    fprintf('Zaber port closed\n')
else 
    fprintf('Zaber port did not exist\n')
end

% %-- Delete the Zabs object once port is closed
clear Zabs