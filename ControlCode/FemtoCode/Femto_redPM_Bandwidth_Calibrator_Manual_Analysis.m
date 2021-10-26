
close all; clear all;

W2uW = 1e6;         % Conversion from Watts to microWatts
%% User inputs
% Folder with Data
fld = '/media/Data_Drive/VFN/TestbedData/210428_COV8/';

% Filename for data
sInd = 5;   % Index of sample to consider
flnm = sprintf('Femto2redPMCalibration_onSuperKEXT_BWs_RedoForNewData_Sample%d.mat',sInd);

% Femto gain factor
gainFact = 10.01;

%-- Manually enter biases -- must be a ROW vector if columns in data
% matrices correspond to different gain settings
biases = [0.007437, 0.008055, 0.008114, 0.008464, 0.008707, 0.008599, 0.010616];

%% Load the file
data_struct = load([fld, flnm]);

%% Extract the data we need
BWs   = data_struct.LASBWs;
lam0  = data_struct.LAScentL;
gainVals = data_struct.gainVals;
meas_FT  = squeeze(data_struct.meas_femto);
% RESCALE the redPM power from watts to microwatts
meas_rP_col  = data_struct.meas_redpm * W2uW;

%% Sort the data in ascending BW order
% This isn't conceptually or algorithmically important but it helps in 
% quickly interpretting and comparing results so that all datasets have the
% same, logical ordering
[BWs, sortind] = sort(BWs);
meas_FT = meas_FT(sortind,:);
meas_rP_col = meas_rP_col(sortind);

%% Clean up the femto data (to get a single, properly scaled voltage vector)
%-- Subtract the bias
meas_FT  = meas_FT - biases;

%-- Rescale all to gain of 6 
meas_FT = meas_FT.*(gainFact.^-(gainVals-6));

%-- Extract one value (gain) to use per row 
    % We only need one value per bandwidth so let's just take the rightmost
    % value in each row so that we have the highest non-nan voltage and are
    % thus in the best response region for the femto
% Preallocate the resulting vector
meas_FT_col = nan(size(meas_FT,1),1);
% Yes, for loop isn't efficient but it's fast enough for small matrices
for rowind = 1:size(meas_FT,1)
    % Find index of the last column that isn't nan in this row
    colind = find(~isnan(meas_FT(rowind,:)),1,'last');
    % In case this BW has no valid values (ie. the whole row is nans since
        % all gains were saturated), just take the first value (a nan)
    if isempty(colind)
        colind = 1;
    end
    % Extract the value and save it in results vector
    meas_FT_col(rowind) = meas_FT(rowind,colind);
end

%% Plot the Femto and redPM values against the sampled BW
figure(); 
yyaxis left
plot(BWs, meas_FT_col, '-o'); 
ylabel('Femto Power (relative to gain=6) [V]')
yyaxis right
plot(BWs, meas_rP_col, '-o');
legend({'Femto', 'redPM'}, 'Location', 'northwest')
xlabel('Bandwidth [nm]')
ylabel('redPM Power [uW]')
title({'Femto and redPM Bandwidth Response',sprintf('lam0 = %0.1f',lam0)})

%% Compute the ratio between redPM and Femto
uW2V = meas_FT_col./meas_rP_col;

%% Plot the resulting BW vs. ratio
figure();
plot(BWs, uW2V, '-o');
ylabel('Response Ratio [V/uW]')
xlabel('Bandwidth [nm]')
title({'Femto-to-redPM Ratio vs. Bandiwdth',sprintf('lam0 = %0.1f',lam0)})

%% Display results in command window
fprintf('\n--- Femto vs. redPM BW Response (Sample %d) ---\n',sInd)
for ind = 1:length(uW2V)
    fprintf('BW: %5.1f  | Femto: %10.6f [V @gain=6] | redPM: %10.6f [uW] | %0.2f [V/uW]\n', ...
            BWs(ind), meas_FT_col(ind), meas_rP_col(ind), uW2V(ind))
end