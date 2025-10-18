% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================
% This script serves as the main orchestrator for the analysis of toroidal
% momentum transport in the TCABR tokamak.
%
% Author: Douglas Oliveira Novaes
% Date: 18-Oct-2025
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

% Add function directories to the MATLAB path
addpath('src', 'plotting', 'utils');

% --- STAGE 1: Setup and Data Loading ---
constants = setup_constants();
exp_data = load_experimental_data(constants);

% --- Verification Step ---
% Display the loaded structures in the Command Window to verify success.
disp('Constants structure loaded:');
disp(constants);
disp('Experimental data structure loaded:');
disp(exp_data);

fprintf('Stage 1 complete. Data is ready for analysis.\n');