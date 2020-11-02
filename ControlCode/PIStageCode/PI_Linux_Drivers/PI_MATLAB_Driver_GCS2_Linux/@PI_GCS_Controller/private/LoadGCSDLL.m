function c = LoadGCSDLL(c)

isWindows   = any (strcmp (mexext, {'mexw32', 'mexw64'}));
isLinux     = any (strcmp (mexext, {'mexglx', 'mexa64'}));



%% Windows
if (isWindows)
    % Retrieve whether MATLAB is used in its 32bit or 64bit version
    
    dllFileNameStub = 'PI_GCS2_DLL';
    
    % Get architecture information
    if(strcmp(mexext , 'mexw64'))
        matlabIs64bit = true;
    else
        matlabIs64bit = false;
    end
    
    % Set dll extension
    if (matlabIs64bit)
        extension = '_x64.dll';
    else
        extension = '.dll';
    end
    
    dllFileName = [dllFileNameStub, extension];
    
    %Retrieve library path
    
    path = '';
    
    % To use the default dll, uncomment the following section
    % Look in regedit entry
    try
        path = winqueryreg('HKEY_LOCAL_MACHINE','SOFTWARE\PI\GCSTranslator','Path');
    catch
        try
            path = winqueryreg('HKEY_LOCAL_MACHINE','SOFTWARE\Wow6432Node\PI\GCSTranslator','Path');
        catch
            % Do nothing
        end
    end
    
    
    if ~isempty ( path )
        % Check whether path name is: 'xxx\xxx' (instead of 'xxx\xxx\')
        if(~strcmp(path(end),'\'))
            path = [path,'\'];
            
            % Check whether path name is: 'xxx\xxx\\' (instead of 'xxx\xxx\')
        elseif((strcmp(path(end-1:end),'\\')))
            path = path(1:end-1);
        end
    end
    
    piMatlabDriverPath = strrep ( mfilename ( 'fullpath' ), '@PI_GCS_Controller\private\LoadGCSDLL', '' );
   
end


%% Linux
if (isLinux)
    
     warning ('Linux is still a beta version! If you encounter any problems please contact support@pi.ws.');
     dllFileName = 'libpi_pi_gcs2.so';
     
     % Get architecture information
     if(strcmp(mexext , 'mexa64'))
         matlabIs64bit = true;
     else
         matlabIs64bit = false;
     end
     
     % Set path
     if (matlabIs64bit)
         path = '/usr/local/PI/lib64/';
     else
         path = '/usr/local/PI/lib32/';
     end
     
     piMatlabDriverPath = strrep ( mfilename ( 'fullpath' ), '@PI_GCS_Controller/private/LoadGCSDLL', '' );
end


%% Build file names

dllFullFileName = [path, dllFileName];
headerFileName  = 'PI_GCS2_DLL.h';      % No path because header file is in the same folder as "@PI_GCS_Controller"

if ~(exist ( dllFullFileName, 'file' ) == 2)
    error ( ['PI MATLAB Driver requires a software component which is not installed on the PC. ' ...
        'Please run "PI_GCS_Library_PI_GCS2_DLL_Setup.exe" ' ...
        'from the product CD to install the required component on the PC.'] );
end


%% Load DLL

% only load dll if it wasn't loaded before
if(~libisloaded(c.libalias))
    disp('Loading PI_MATLAB_Driver_GCS2 ...');
    
    directoryOfCallingScript = cd;
    
    try
        cd (piMatlabDriverPath);
            
        warning('off', 'all'); % If not disabled missing linux only functions will rise a warning
        
        % create prototype if it doesn't exist yet.
        if ( ~ (exist ( 'PI_GCS2_DLL_prototype', 'file' ) ||  exist ( 'PI_GCS2_DLL_prototype_x64', 'file' ) ) )
            if ( matlabIs64bit )
                [notfound,warnings] = loadlibrary (dllFullFileName, headerFileName, 'mfilename', 'PI_GCS2_DLL_prototype_x64', 'alias', c.libalias);
            else
        	    [notfound,warnings] = loadlibrary (dllFullFileName, headerFileName, 'mfilename', 'PI_GCS2_DLL_prototype', 'alias', c.libalias);
            end
            unloadlibrary('PI');
        end
        
        if ( matlabIs64bit )
            [notfound,warnings] = loadlibrary (dllFullFileName, @PI_GCS2_DLL_prototype_x64, 'alias', c.libalias);
        else
            [notfound,warnings] = loadlibrary (dllFullFileName, @PI_GCS2_DLL_prototype, 'alias', c.libalias);
        end
        warning('on', 'all');
        
    catch ME
        errorMessage = sprintf('Error while loading %s.', dllFileName);
        errorMessage = sprintf('%s\n\nThe %s was assumed to have the following path:\n%s', errorMessage, dllFileName, dllFullFileName);
        errorMessage = sprintf('%s\n\nThe PI MATLAB DRIVER GCS2 either could not find or could not use the "%s".', errorMessage, dllFileName);
        errorMessage = sprintf('%s\nMake sure you have installed the newest driver for your PI controller and look at the PI MATLAB DRIVER GCS2 Manual in chapter "Troubleshooting" headword "PI GCS2 DLL" for known problems and fixes for problems while loading the %s.', errorMessage, dllFileName);
        errorMessage = sprintf('%s\n\nAdditional Information: MATLAB itself raised the following error:', errorMessage);
        errorMessage = sprintf('%s\n%s', errorMessage, ME.message);
        error(errorMessage)
    end
    
    cd (directoryOfCallingScript);
end

disp('PI_MATLAB_Driver_GCS2 loaded successfully.');

c.dllfunctions = libfunctions(c.libalias);



