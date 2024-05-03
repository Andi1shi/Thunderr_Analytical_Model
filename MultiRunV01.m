%% Thunderstorm Simulation and Optimization Script
% Code Description

% This script simulates horizontal mean wind velocity caused by traveling downbursts, combining three velocity fields: 
% radial impinging jet, downburst translational, and ambient ABL wind velocities. It employs a parametric-analytical model 
% for spatiotemporal reconstruction of downburst events, optimizing downburst dimensions and kinematics against real measurements 
% using the Teaching-Learning-Based Optimization (TLBO) algorithm. 

% References:
% [1] - Xhelaj, A., Burlando, M., Solari, G. (2020). A general-purpose analytical model for 
% reconstructing the thunderstorm outflows of travelling downbursts immersed in ABL flows.
% Journal of Wind Engineering and Industrial Aerodynamics, 207, 104373.

%[2] -  Xhelaj, A., Burlando, M.: Application of metaheuristic optimization algorithms to evaluate the geometric and
% kinematic parameters of downbursts. Advances in Engineering Software. Volume 173, November 2022, 103203.,
% 2022.

% Author: [Andi Xhelaj]
% Date: [29/03/2024]
% Code Version: MultiRunV01


%% Initialization
clc;          % Clear the command window
clear;       % Clear the workspace
close all;  % Close all figures
FontSize = 14; % For plotting purposes

%% Load Data
% Load the thunderstorm data from a .mat file
load('Thunderstorm_Romania_MA.mat');
fs = 10;       % Sampling frequency of recorded data [Hz]
dtr = 1 / fs;  % Sampling period of the recorded data [s]


%% Simulation Parameters
% Section for defining simulation time and station location parameters.
% Users can modify these parameters as per their requirements.

% Simulation Time Settings
ti = 0;    % Initial simulation time [sec]
tf = 3600; % Final simulation time [sec]
dt = 1;    % Simulation time step [sec]

t = (ti:dt:tf-dt)';    % Time vector
nt = dt * fs;          % For time conversion with the recorded data

% Station Location
xp = 0;  % x-component of the station location [meters]
yp = 0;  % y-component of the station location [meters]


%% Parameter Bounds
% Define lower and upper bounds for various parameters involved in the simulation.
% These bounds are crucial for the optimization algorithm to find the best parameters.

%% Downburst Touchdown at t = 0
% Parameters for the initial downburst touchdown location
xc0_lb = -10000;  % Lower bound for x-component [meters]
xc0_ub = 10000;   % Upper bound for x-component [meters]
yc0_lb = -10000;  % Lower bound for y-component [meters]
yc0_ub = 10000;   % Upper bound for y-component [meters]

%% Stationary Downburst Parameters
R_lb = 400;       % Downburst downdraft radius lower bound [meters]
R_ub = 1000;      % Downburst downdraft radius upper bound [meters]
psi_lb = 1.6;     % Radius of maximum radial velocity lower bound (Rmax_lb/R)
psi_ub = 2.6;     % Radius of maximum radial velocity upper bound (Rmax_ub/R)
v1_lb = 0;        % Downburst maximum radial velocity lower bound [m/s]
vl_ub = 40;       % Downburst maximum radial velocity upper bound [m/s]

% Intensity Decay Function Parameters
T1_lb = 2 * 60;     % Period of linear intensification lower bound [seconds]
T1_ub = 15 * 60;  % Period of linear intensification upper bound [seconds]
Tf_lb = 15 * 60;    % Total duration of the downburst lower bound [seconds]
Tf_ub = 60 * 60;   % Total duration of the downburst upper bound [seconds]

%% Translating Downburst Parameters
% Parameters for the movement of the downburst
v2_lb = 0;            % Downburst translation velocity lower bound [m/s]
v2_ub = 40;         % Downburst translation velocity upper bound [m/s]
beta2_lb = 0;        % Downburst translation direction lower bound [degrees from East]
beta2_ub = 359;   % Downburst translation direction upper bound [degrees from East]

%% ABL Background Wind Parameters
% Ambient Boundary Layer (ABL) wind conditions outside the downburst
v3_lb = 0;              % ABL background wind speed lower bound [m/s]
v3_ub = 40;           % ABL background wind speed upper bound [m/s]
beta3_lb = 0;          % ABL wind direction lower bound [degrees from East]
beta3_ub = 359;      % ABL wind direction upper bound [degrees from East]

%% Aggregate Lower and Upper Bounds
% Compiling all lower and upper bounds into vectors for optimization algorithm input
lb = [xc0_lb, yc0_lb, R_lb, psi_lb, v1_lb, T1_lb, Tf_lb, v2_lb, beta2_lb, v3_lb, beta3_lb];  % Lower bounds vector
ub = [xc0_ub, yc0_ub, R_ub, psi_ub, vl_ub, T1_ub, Tf_ub, v2_ub, beta2_ub, v3_ub, beta3_ub];  % Upper bounds vector
D = length(lb); % Dimension of the problem based on the length of lower bounds vector

%% Optimization Algorithm Configuration
% This section configures parameters for the optimization algorithm, including weights for the objective function and settings for algorithm iterations.

%% Population and Iteration Settings
% Settings for the size of the population and the number of iterations, which are critical parameters for the optimization process.
Npop = 25;  % Size of population (class size), default is 50
T = 100;        % Total number of iterations in TLBO (Teaching-Learning-Based Optimization), default is 100
%% Multiple Runs of Optimization Algorithm
% This section executes the optimization algorithm multiple times to ensure robustness and to identify the best solution across runs.

% Configuration for multiple runs
Nruns = 2^5; % Total number of runs, here 2^4 = 16

%% Objective Function Weights
% Definition of weights for the components of the Normalized Mean Square Error (NMSE) in the objective function.
kv = 1;      % Weight coefficient for NMSE_v, range: 0 <= kv <= 1
kalpha = 1;  % Weight coefficient for NMSE_alpha, range: 0 <= kalpha <= 1

% Objective Function
% NMSE = kv * NMSE_v + kalpha * NMSE_alpha;
% The objective function combines NMSE_v and NMSE_alpha weighted by kv and kalpha, respectively.
FITNESSFCN = @THUNDERR_NMSE;  % Function handle of the fitness function

% Total Number of Objective Function Evaluations
% The total number of times the objective function will be evaluated.
% It is calculated based on the size of the population and the number of iterations.
FE = Npop + 2 * Npop * T;  % Total number of evaluations of the objective function

% Initialization of matrices to store results
BestSolMat = NaN(Nruns, D); % Matrix to store the best solution for each run
BestFitness = NaN(Nruns, 1); % Vector to store the best fitness value for each run
BestFitIterMat = NaN(Nruns, T+1); % Matrix to store the fitness over iterations for each run

% Parallel loop for multiple runs
% Start timing the runs
tic

parfor i = 1:Nruns
    rng(i, 'twister') % Ensure reproducibility of the results by setting a unique seed for each iteration
    [BestSolMat(i, :), BestFitness(i), BestFitIterMat(i, :), Pop, obj] = TLBO_OPTIMIZATION(FITNESSFCN, lb, ub, Npop, T, ti, dt, tf, xp, yp, vrt, alphart, fs, kv, kalpha);
end

% Stop timing the runs
toc

%% Post-processing and Analysis of Results
% Analyze and display the results from multiple runs for insights into the optimization performance.

% Sorting solutions from worst to best based on fitness
[BestFitness, idx] = sort(BestFitness, 'descend');
BestSolMat = BestSolMat(idx, :);

% Statistical analysis of the Best Fitness across runs
Stat(1) = min(BestFitness); % Minimum (best) fitness function value
Stat(2) = max(BestFitness); % Maximum (worst) fitness function value
Stat(3) = mean(BestFitness); % Mean fitness function value
Stat(4) = std(BestFitness); % Standard deviation of fitness function values
sol_min_NMSE = BestSolMat(end, :); % Solution with minimum NMSE

% Calculating the mean best fitness over iterations
BestFitIterMean = mean(BestFitIterMat);

% Displaying summary of results
disp(['Total number of fitness function evaluations = ', num2str(Nruns * FE)])
disp(['Best fitness function value = ', num2str(Stat(1))])
disp(['Worst fitness function value = ', num2str(Stat(2))])
disp(['Mean fitness function value = ', num2str(Stat(3))])
disp(['Standard deviation of fitness function values = ', num2str(Stat(4))])
disp('Best solution parameters:')
disp(['xc0 = ', num2str(sol_min_NMSE(1)), ' yc0 = ', num2str(sol_min_NMSE(2)), ' R = ', num2str(sol_min_NMSE(3)), ' Rmax = ', num2str(sol_min_NMSE(3) * sol_min_NMSE(4)), ' v1 = ', num2str(sol_min_NMSE(5)), ' T1 = ', num2str(sol_min_NMSE(6)), ' Tf = ', num2str(sol_min_NMSE(7)), ' v2 = ', num2str(sol_min_NMSE(8)), ' beta2 = ', num2str(sol_min_NMSE(9)), ' v3 = ', num2str(sol_min_NMSE(10)), ' beta3 = ', num2str(sol_min_NMSE(11))]);

%% Convergence with semilog plot - Algorithm Performance
CreateEnvelope(T,Nruns,BestFitness,BestFitIterMat,FontSize)

%% Schematic plotting and storing the best solution
M = 3; % Number of last best solutions to plot

% Loop over the last 3 best solutions
for k = Nruns-M+1:Nruns
[vxt_m,vyt_m,vt_m,alphat_m] = THUNDERR_PLOT(BestSolMat(k,:),BestFitness(k),ti,dtr,tf,xp,yp,vrt,alphart,fs,lb(1),ub(1),lb(2),ub(2),FontSize);
end
