function versionNumber = GetVersionNumber(c)

driverDirectoryOnWindows = 'C:\Users\Public\PI\PI_MATLAB_Driver_GCS2';
fileIdRead  = fopen ( [c.driverDirectoryOnWindows, '\', 'Version.txt'],  'r' );
currentLine = fgets(fileIdRead);
versionNumber = '';

while ( currentLine > 0)
    if ( ~isempty ( regexp ( currentLine, 'Version:\s*', 'match' ) ) )
        versionNumber = regexp ( currentLine, 'Version:\s*', 'split' );
        versionNumber = versionNumber{end};
        break;
    end    
    
    currentLine = fgets(fileIdRead);
end

versionNumber = ['PI MATLAB Driver Version Number: ', versionNumber]; 