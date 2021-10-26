
close all; clear all;

W2uW = 1e6;         % Conversion from Watts to microWatts
%% User inputs
% Folder with Data
fld = '/media/Data_Drive/VFN/TestbedData/210428_COV8/';

% Filename for datasets
sInd = 1;   % Index of sample to consider
nmBWs_rP = sprintf('Femto2redPMCalibration_onSuperKEXT_BWs_redPM_Sample%d',sInd);
nmBWs_FT = sprintf('Femto2redPMCalibration_onSuperKEXT_BWs_Femto_Sample%d',sInd);

% Femto gain factor
gainFact = 10.01;

%-- Manually enter biases -- must be a ROW vector if columns in data
% matrices correspond to different gain settings
biases = [0.007437, 0.008055, 0.008114, 0.008464, 0.008707, 0.008599, 0.010616];

%% Load the data
redPM = load([fld, nmBWs_rP]);
Femto = load([fld, nmBWs_FT]);

%% Extract the data we need

% NOTE:: we'll assume the datasets are samples of identical conditions
  % If you don't want to make this assumption, you'll need to do something
  % like the following to check for every variable.
% % First check that they're all the same
% if ~isequal(redPM.LASBWs, Femto.LASBWs)
%      error('LASBWs are not equal')
% end

%-- Assuming the datasets are correct, pull out the pertinent values from
%   just one of the two as needed
LASBWs   = redPM.LASBWs;
LAScentL = redPM.LAScentL;
LASpower = redPM.LASpower;
gainVals = Femto.gainVals;
meas_FT  = squeeze(Femto.meas);
% RESCALE the redPM power from watts to microwatts
meas_rP_col  = redPM.meas * W2uW;

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
    % Extract the value and save it in results vector
    meas_FT_col(rowind) = meas_FT(rowind,colind);
end

%% Plot the Femto and redPM values against the sampled BW
figure(); 
yyaxis left
plot(LASBWs, meas_FT_col, '-o'); 
ylabel('Femto Power (relative to gain=6) [V]')
yyaxis right
plot(LASBWs, meas_rP_col, '-o');
legend({'Femto', 'redPM'}, 'Location', 'northwest')
xlabel('Bandwidth [nm]')
ylabel('redPM Power [uW]')
title({'Femto and redPM Bandwidth Response',sprintf('lam0 = %0.1f, superK pwr = %0.1f',LAScentL,LASpower)})

%% Compute the ratio between redPM and Femto
uW2V = meas_FT_col./meas_rP_col;

%% Plot the resulting BW vs. ratio
figure();
plot(LASBWs, uW2V, '-o');
ylabel('Response Ratio [V/uW]')
xlabel('Bandwidth [nm]')
title({'Femto-to-redPM Ratio vs. Bandiwdth',sprintf('lam0 = %0.1f, superK pwr = %0.1f',LAScentL,LASpower)})

%% Display results in command window
fprintf('\n--- Femto vs. redPM BW Response (Sample %d) ---\n',sInd)
for ind = 1:length(uW2V)
    fprintf('BW: %5.1f  | Femto: %9.6f [V @gain=6] | redPM: %9.6f [uW] | %0.2f [V/uW]\n', ...
            LASBWs(ind), meas_FT_col(ind), meas_rP_col(ind), uW2V(ind))
end