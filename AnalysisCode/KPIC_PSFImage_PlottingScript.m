% clean workspace (optional)
close all; clear all;

%% Analysis Setup

%-- Directories
% Raw data directory (where raw data is located for import)** End with filesep**
rwfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\PupilVFN\PSFImaging\KBand_Vortex_Apod_Good_072519\';
% Output data directory (where processed data will be saved)
svfld = [rwfld 'Processed\'];

%-- Control Flags
isSaveFigs = true;  % Flag to save figures
isPlotFull = false; % Flag to display full (un-cropped) frames
isPlotLogs = true;  % Flag to display the log of the cropped images

%-- Outlier tolerance (for removing frames with flux changes)
tol = 5;

%-- Crop Parameters
vort1_crp   = 60;           % Crop window (+/- this value)
vort1_cnt = [171 242];      % Crop center for IN frame. [row col]
vort2_crp   = vort1_crp;    % Crop window (+/- this value)
vort2_cnt = [81 232];      % Crop center for IN frame. [row col]
vort3_crp   = 51;    % Crop window (+/- this value)
vort3_cnt = [86 224];      % Crop center for IN frame. [row col]
apod_crp = 26;
% _Out is automatically centered at peak

%-- Choose colormap for saved images
cmap = 'gray';
%cmap = 'parula';

%{ 
NOTE: For some reason, the camera code used to create these fits files does not
save the image data into the "Primary" extension of the fits file. It saves the
images into the "image" extension so we need to keep this in mind when using
fitsread to load the images.
%}
%% Vortex 1: DS081111-B  (12.5 x 12.5 mm)
fprintf('\n______VORT 1 DS081111-B  (12.5 x 12.5 mm)______\n')

%-- Load the IN file
flnm = [rwfld 'KPIC_Vort1_VortIN_PSF2.fits'];
vort1In = fitsread(flnm, 'image');

%-- Load the Out file
flnm = [rwfld 'KPIC_Vort1_VortOut_PSF2.fits'];
vort1Out = fitsread(flnm, 'image');

%-- Load the Background file
flnm = [rwfld 'KPIC_Vort1_Background_PSF2.fits'];
vort1Bck = fitsread(flnm, 'image');

%-- Remove outliers (where flux changed significantly)
% Remove from IN file
mns = squeeze(mean(vort1In, [1 2]));    % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort1In = vort1In(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort1In cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort1In, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
% Remove from Out file
mns = squeeze(mean(vort1Out, [1 2]));   % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort1Out = vort1Out(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort1Out cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort1Out, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
%-- NO NEED TO REMOVE OUTLIERS FROM BACKGROUND SINCE LASER NOT PRESENT
mns = squeeze(mean(vort1Bck, 'all'));
fprintf('-- No Outliers Removed From Vort1 Background cube\n')
fprintf('   Mean of background:   %8.4f\n', mns);

%-- Average outlier-removed data
vort1In  = mean(vort1In, 3);
vort1Out = mean(vort1Out, 3); 
vort1Bck = mean(vort1Bck, 3);

%-- Subtract Background 
vort1In_bk  = vort1In - vort1Bck;
vort1Out_bk = vort1Out - vort1Bck;

if isPlotFull
    %-- Display data [FULL FRAME]
    % In data
    figure;
    imagesc(vort1In_bk);
    axis image
    colorbar;
    title('Vortex 1 (DS081111-B) IN')
    % Out data
    figure;
    imagesc(vort1Out_bk);
    axis image
    colorbar;
    title('Vortex 1 (DS081111-B) Out')
    % Background
    figure;
    imagesc(log10(vort1Bck));
    axis image
    colorbar;
    title('Vortex 1 (DS081111-B) LOG10(Background)')
end

%-- Display data [CROPPED FRAME]
% In data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
crp(1) = max(vort1_cnt(1)-vort1_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(vort1_cnt(1)+vort1_crp, size(vort1In_bk, 1));
crp(3) = max(vort1_cnt(2)-vort1_crp, 0);
crp(4) = min(vort1_cnt(2)+vort1_crp, size(vort1In_bk, 2));
vort1In_crp = vort1In_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort1In_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 1 (DS081111-B) IN - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort1In_Crp_BckSub.png'],'-r300', '-painters');
end
% Out data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
[~, cntI] = max(vort1Out_bk(:));
[cnt(1), cnt(2)] = ind2sub(size(vort1Out_bk), cntI);
crp(1) = max(cnt(1)-vort1_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(cnt(1)+vort1_crp, size(vort1Out_bk, 1));
crp(3) = max(cnt(2)-vort1_crp, 0);
crp(4) = min(cnt(2)+vort1_crp, size(vort1Out_bk, 2));
vort1Out_crp = vort1Out_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort1Out_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 1 (DS081111-B) Out - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort1Out_Crp_BckSub.png'],'-r300', '-painters');
end
% Background
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
vort1Bck_crp = vort1Bck(crp(1):crp(2),crp(3):crp(4),:); % Crop to "Out" dims
imagesc(log10(vort1Bck_crp));
axis image
colormap(cmap)
colorbar;
title('Vortex 1 (DS081111-B) LOG10(Background) - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort1Bck_Crp.png'],'-r300', '-painters');
end

if isPlotLogs
    % In data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort1In_crpLog = vort1In_crp;
    vort1In_crpLog(vort1In_crpLog<=0) = min(vort1In_crp(vort1In_crp>0));
    imagesc(log10(vort1In_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 1 (DS081111-B) LOG10(IN) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort1In_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
    % Out data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort1Out_crpLog = vort1Out_crp;
    vort1Out_crpLog(vort1Out_crpLog<=0) = min(vort1Out_crp(vort1Out_crp>0));
    imagesc(log10(vort1Out_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 1 (DS081111-B) LOG10(Out) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort1Out_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
end

%% Vortex 2: DS110414-10  (12.5 x 12.5 mm)
fprintf('\n______VORT 2 DS110414-10  (12.5 x 12.5 mm)______\n')

%-- Load the IN file
flnm = [rwfld 'KPIC_Vort2_VortIN_PSF2.fits'];
vort2In = fitsread(flnm, 'image');

%-- Load the Out file
flnm = [rwfld 'KPIC_Vort2_VortOut_PSF2.fits'];
vort2Out = fitsread(flnm, 'image');

%-- Load the Background file
flnm = [rwfld 'KPIC_Vort2_Background_PSF2.fits'];
vort2Bck = fitsread(flnm, 'image');

%-- Remove outliers (where flux changed significantly)
% Remove from IN file
mns = squeeze(mean(vort2In, [1 2]));    % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort2In = vort2In(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort2In cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort2In, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
% Remove from Out file
mns = squeeze(mean(vort2Out, [1 2]));   % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort2Out = vort2Out(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort2Out cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort2Out, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
%-- NO NEED TO REMOVE OUTLIERS FROM BACKGROUND SINCE LASER NOT PRESENT
mns = squeeze(mean(vort2Bck, 'all'));
fprintf('-- No Outliers Removed From Vort2 Background cube\n')
fprintf('   Mean of background:   %8.4f\n', mns);

%-- Average outlier-removed data
vort2In  = mean(vort2In, 3);
vort2Out = mean(vort2Out, 3); 
vort2Bck = mean(vort2Bck, 3);

%-- Subtract Background 
vort2In_bk  = vort2In - vort2Bck;
vort2Out_bk = vort2Out - vort2Bck;

if isPlotFull
    %-- Display data [FULL FRAME]
    % In data
    figure;
    imagesc(vort2In_bk);
    axis image
    colorbar;
    title('Vortex 2 (DS110414-10) IN')
    % Out data
    figure;
    imagesc(vort2Out_bk);
    axis image
    colorbar;
    title('Vortex 2 (DS110414-10) Out')
    % Background
    figure;
    imagesc(log10(vort2Bck));
    axis image
    colorbar;
    title('Vortex 2 (DS110414-10) LOG10(Background)')
end

%-- Display data [CROPPED FRAME]
% In data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
crp(1) = max(vort2_cnt(1)-vort2_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(vort2_cnt(1)+vort2_crp, size(vort2In_bk, 1));
crp(3) = max(vort2_cnt(2)-vort2_crp, 0);
crp(4) = min(vort2_cnt(2)+vort2_crp, size(vort2In_bk, 2));
vort2In_crp = vort2In_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort2In_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 2 (DS110414-10) IN - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort2In_Crp_BckSub.png'],'-r300', '-painters');
end
% Out data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
[~, cntI] = max(vort2Out_bk(:));
[cnt(1), cnt(2)] = ind2sub(size(vort2Out_bk), cntI);
crp(1) = max(cnt(1)-vort2_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(cnt(1)+vort2_crp, size(vort2Out_bk, 1));
crp(3) = max(cnt(2)-vort2_crp, 0);
crp(4) = min(cnt(2)+vort2_crp, size(vort2Out_bk, 2));
vort2Out_crp = vort2Out_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort2Out_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 2 (DS110414-10) Out - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort2Out_Crp_BckSub.png'],'-r300', '-painters');
end
% Background
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
vort2Bck_crp = vort2Bck(crp(1):crp(2),crp(3):crp(4),:); % Crop to "Out" dims
imagesc(log10(vort2Bck_crp));
axis image
colormap(cmap)
colorbar;
title('Vortex 2 (DS110414-10) LOG10(Background) - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort2Bck_Crp.png'],'-r300', '-painters');
end

if isPlotLogs
    % In data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort2In_crpLog = vort2In_crp;
    vort2In_crpLog(vort2In_crpLog<=0) = min(vort2In_crp(vort2In_crp>0));
    imagesc(log10(vort2In_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 2 (DS110414-10) LOG10(IN) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort2In_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
    % Out data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort2Out_crpLog = vort2Out_crp;
    vort2Out_crpLog(vort2Out_crpLog<=0) = min(vort2Out_crp(vort2Out_crp>0));
    imagesc(log10(vort2Out_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 2 (DS110414-10) LOG10(Out) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort2Out_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
end

%% Vortex 3: DS110419-20  (12.5 x 12.5 mm)
fprintf('\n______VORT 3 DS110419-20  (45 x 45 mm)______\n')

%-- Load the IN file
flnm = [rwfld 'KPIC_Vort3_VortIN_PSF2.fits'];
vort3In = fitsread(flnm, 'image');

%-- Load the Out file
flnm = [rwfld 'KPIC_Vort3_VortOut_PSF2.fits'];
vort3Out = fitsread(flnm, 'image');

%-- Load the Background file
flnm = [rwfld 'KPIC_Vort3_Background_PSF2.fits'];
vort3Bck = fitsread(flnm, 'image');

%-- Remove outliers (where flux changed significantly)
% Remove from IN file
mns = squeeze(mean(vort3In, [1 2]));    % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort3In = vort3In(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort3In cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort3In, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
% Remove from Out file
mns = squeeze(mean(vort3Out, [1 2]));   % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
vort3Out = vort3Out(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from Vort3Out cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(vort3Out, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
%-- NO NEED TO REMOVE OUTLIERS FROM BACKGROUND SINCE LASER NOT PRESENT
mns = squeeze(mean(vort3Bck, 'all'));
fprintf('-- No Outliers Removed From Vort3 Background cube\n')
fprintf('   Mean of background:   %8.4f\n', mns);

%-- Average outlier-removed data
vort3In  = mean(vort3In, 3);
vort3Out = mean(vort3Out, 3); 
vort3Bck = mean(vort3Bck, 3);

%-- Subtract Background 
vort3In_bk  = vort3In - vort3Bck;
vort3Out_bk = vort3Out - vort3Bck;

if isPlotFull
    %-- Display data [FULL FRAME]
    % In data
    figure;
    imagesc(vort3In_bk);
    axis image
    colorbar;
    title('Vortex 3 (DS110419-20) IN')
    % Out data
    figure;
    imagesc(vort3Out_bk);
    axis image
    colorbar;
    title('Vortex 3 (DS110419-20) Out')
    % Background
    figure;
    imagesc(log10(vort3Bck));
    axis image
    colorbar;
    title('Vortex 3 (DS110419-20) LOG10(Background)')
end

%-- Display data [CROPPED FRAME]
% In data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
crp(1) = max(vort3_cnt(1)-vort3_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(vort3_cnt(1)+vort3_crp, size(vort3In_bk, 1));
crp(3) = max(vort3_cnt(2)-vort3_crp, 0);
crp(4) = min(vort3_cnt(2)+vort3_crp, size(vort3In_bk, 2));
vort3In_crp = vort3In_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort3In_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 3 (DS110419-20) IN - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort3In_Crp_BckSub.png'],'-r300', '-painters');
end
% Out data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
[~, cntI] = max(vort3Out_bk(:));
[cnt(1), cnt(2)] = ind2sub(size(vort3Out_bk), cntI);
crp(1) = max(cnt(1)-vort3_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(cnt(1)+vort3_crp, size(vort3Out_bk, 1));
crp(3) = max(cnt(2)-vort3_crp, 0);
crp(4) = min(cnt(2)+vort3_crp, size(vort3Out_bk, 2));
vort3Out_crp = vort3Out_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(vort3Out_crp);
axis image
colormap(cmap)
colorbar;
title('Vortex 3 (DS110419-20) Out - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort3Out_Crp_BckSub.png'],'-r300', '-painters');
end
% Background
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
vort3Bck_crp = vort3Bck(crp(1):crp(2),crp(3):crp(4),:); % Crop to "Out" dims
imagesc(log10(vort3Bck_crp));
axis image
colormap(cmap)
colorbar;
title('Vortex 3 (DS110419-20) LOG10(Background) - Cropped')
if isSaveFigs
    export_fig([svfld 'Vort3Bck_Crp.png'],'-r300', '-painters');
end

if isPlotLogs
    % In data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort3In_crpLog = vort3In_crp;
    vort3In_crpLog(vort3In_crpLog<=0) = min(vort3In_crp(vort3In_crp>0));
    imagesc(log10(vort3In_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 3 (DS110419-20) LOG10(IN) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort3In_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
    % Out data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    vort3Out_crpLog = vort3Out_crp;
    vort3Out_crpLog(vort3Out_crpLog<=0) = min(vort3Out_crp(vort3Out_crp>0));
    imagesc(log10(vort3Out_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Vortex 3 (DS110419-20) LOG10(Out) - Cropped')
    if isSaveFigs
        export_fig([svfld 'Vort3Out_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
end


%% Apodizer: 
fprintf('\n______ Apoodizer ______\n')

%-- Load the IN file
flnm = [rwfld 'KPIC_Apod_ApodIn_PSF2.fits'];
apodIn = fitsread(flnm, 'image');

%-- Load the Out file
flnm = [rwfld 'KPIC_Apod_ApodOut_PSF2.fits'];
apodOut = fitsread(flnm, 'image');

%-- Load the Background file
flnm = [rwfld 'KPIC_Apod_Background_PSF2.fits'];
apodBck = fitsread(flnm, 'image');

%-- Remove outliers (where flux changed significantly)
% Remove from IN file
mns = squeeze(mean(apodIn, [1 2]));    % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
apodIn = apodIn(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from ApodIn cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(apodIn, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
% Remove from Out file
mns = squeeze(mean(apodOut, [1 2]));   % Average frames individually 
mn2 = mean(mns);                        % Average all frames to get outliers
excl = (mns > mn2 + tol)+(mns < mn2 - tol); % Find outliers (+ for logical or)
apodOut = apodOut(:,:,~logical(excl));  % Extract non-outlier frames
fprintf('-- Removed %d Outlier Frames from ApodOut cube\n', sum(excl))
fprintf('   Original STD of mean: %8.4f\n', std(mns));
mns = squeeze(mean(apodOut, [1 2]));    % Calculate new std
fprintf('   STD of new mean:      %8.4f\n', std(mns));
fprintf('   New mean:             %8.4f\n', mean(mns));
%-- NO NEED TO REMOVE OUTLIERS FROM BACKGROUND SINCE LASER NOT PRESENT
mns = squeeze(mean(apodBck, 'all'));
fprintf('-- No Outliers Removed From Apod Background cube\n')
fprintf('   Mean of background:   %8.4f\n', mns);

%-- Average outlier-removed data
apodIn  = mean(apodIn, 3);
apodOut = mean(apodOut, 3); 
apodBck = mean(apodBck, 3);

%-- Subtract Background 
apodIn_bk  = apodIn - apodBck;
apodOut_bk = apodOut - apodBck;

if isPlotFull
    %-- Display data [FULL FRAME]
    % In data
    figure;
    imagesc(apodIn_bk);
    axis image
    colorbar;
    title('Apodizer IN')
    % Out data
    figure;
    imagesc(apodOut_bk);
    axis image
    colorbar;
    title('Apodizer Out')
    % Background
    figure;
    imagesc(log10(apodBck));
    axis image
    colorbar;
    title('Apodizer LOG10(Background)')
end

%-- Display data [CROPPED FRAME]
% In data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
[~, cntI] = max(apodIn_bk(:));
[apod_cnt(1), apod_cnt(2)] = ind2sub(size(apodIn_bk), cntI);
crp(1) = max(apod_cnt(1)-apod_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(apod_cnt(1)+apod_crp, size(apodIn_bk, 1));
crp(3) = max(apod_cnt(2)-apod_crp, 0);
crp(4) = min(apod_cnt(2)+apod_crp, size(apodIn_bk, 2));
apodIn_crp = apodIn_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(apodIn_crp);
axis image
colormap(cmap)
colorbar;
title('Apodizer IN - Cropped')
if isSaveFigs
    export_fig([svfld 'ApodIn_Crp_BckSub.png'],'-r300', '-painters');
end
% Out data
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
[~, cntI] = max(apodOut_bk(:));
[cnt(1), cnt(2)] = ind2sub(size(apodOut_bk), cntI);
crp(1) = max(cnt(1)-apod_crp, 0);  % Make sure crp vals are within frame bounds
crp(2) = min(cnt(1)+apod_crp, size(apodOut_bk, 1));
crp(3) = max(cnt(2)-apod_crp, 0);
crp(4) = min(cnt(2)+apod_crp, size(apodOut_bk, 2));
apodOut_crp = apodOut_bk(crp(1):crp(2),crp(3):crp(4),:);
imagesc(apodOut_crp);
axis image
colormap(cmap)
colorbar;
title('Apodizer Out - Cropped')
if isSaveFigs
    export_fig([svfld 'ApodOut_Crp_BckSub.png'],'-r300', '-painters');
end
% Background
figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
apodBck_crp = apodBck(crp(1):crp(2),crp(3):crp(4),:); % Crop to "Out" dims
imagesc(log10(apodBck_crp));
axis image
colormap(cmap)
colorbar;
title('Apodizer LOG10(Background) - Cropped')
if isSaveFigs
    export_fig([svfld 'ApodBck_Crp.png'],'-r300', '-painters');
end

if isPlotLogs
    % In data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    apodIn_crpLog = apodIn_crp;
    apodIn_crpLog(apodIn_crpLog<=0) = min(apodIn_crp(apodIn_crp>0));
    imagesc(log10(apodIn_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Apodizer LOG10(IN) - Cropped')
    if isSaveFigs
        export_fig([svfld 'ApodIn_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
    % Out data
    figure('color','white','units', 'inches', 'Position', [0 0 6 6]);
    apodOut_crpLog = apodOut_crp;
    apodOut_crpLog(apodOut_crpLog<=0) = min(apodOut_crp(apodOut_crp>0));
    imagesc(log10(apodOut_crpLog));
    axis image
    colormap(cmap)
    colorbar;
    title('Apodizer LOG10(Out) - Cropped')
    if isSaveFigs
        export_fig([svfld 'ApodOut_Crp_BckSub_Log.png'],'-r300', '-painters');
    end
end