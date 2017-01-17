clear; close all; clc;
tic

%% Define MCMC parameters for JAGS
nchains  = 4;       % How Many Chains?
nburnin  = 1000;    % How Many Burn-in Samples?
nsamples = 50000;   % How Many Recorded Samples?
nthin    = 5;       % How Many Samples to Thin by?

%% Load text file with model specification
model_name = 'multiBayesianMCMC.txt';

%% Load observed data
self = 0;
if self   
    group_data = csvread('DATA_LOCATION/DATA_FILENAME_CHOOSE_SELF.csv', 1, 0);
else
    group_data = csvread('DATA_LOCATION/DATA_FILENAME_CHOOSE_OTHER.csv', 1, 0);
end

%% Index experimental variables of interest
col_rating_self = 16;
col_p_self = 1;
col_rating_other = 17;
col_p_other = 2;

%% Set initial values of each parameter in each chain
for i = 1 : nchains
    S.alphagmean = 0;
    S.beta1gmean = 0;
    S.beta2gmean = 0;
    S.beta3gmean = 0;
    S.beta4gmean = 0;
    S.beta5gmean = 0;
    S.alphagprec = 0.001;
    S.beta1gprec = 0.001;
    S.beta2gprec = 0.001;
    S.beta3gprec = 0.001;
    S.beta4gprec = 0.001;
    S.beta5gprec = 0.001;
    S.pmod = [.25 .25 .25 .25];
    init0(i) = S; 
end

%% Parallel computing
doparallel = 1;

if doparallel
   if isempty(gcp('nocreate'))
      pool = parpool(4);
   end
end

%% Prepare input data

% Remove missed trials (RT = 0, column 9)
group_data(group_data(:, 9) == 0, :) = [];
% Number of data points
N = max(size(group_data));
% List of subject indices
subj = group_data(:, 22);
% Dirichlet priors for hierarchical model probabilites
dd = [1 1 1 1];

% Prepare input data with or without centering and scaling
    
do_zscore = 1;

if do_zscore
    if self
       datastruct = struct('y', group_data(:, 5), 'N', N, 'rating', zscore(group_data(:, col_rating_self)), 'prob', zscore(group_data(:, col_p_self)), 'ratingirr', zscore(group_data(:, col_rating_other)), 'probirr', zscore(group_data(:, col_p_other)), 'dd', dd, 'subj', subj);
    else
       datastruct = struct('y', group_data(:, 5), 'N', N, 'rating', zscore(group_data(:, col_rating_other)), 'prob', zscore(group_data(:, col_p_other)), 'ratingirr', zscore(group_data(:, col_rating_self)), 'probirr', zscore(group_data(:, col_p_self)), 'dd', dd, 'subj', subj);
    end
else
    if self
       datastruct = struct('y', group_data(:, 5), 'N', N, 'rating', (group_data(:, col_rating_self)), 'prob', (group_data(:, col_p_self)), 'ratingirr', (group_data(:, col_rating_other)), 'probirr', (group_data(:, col_p_other)), 'dd', dd, 'subj', subj);
    else
       datastruct = struct('y', group_data(:, 5), 'N', N, 'rating', (group_data(:, col_rating_other)), 'prob', (group_data(:, col_p_other)), 'ratingirr', (group_data(:, col_rating_self)), 'probirr', (group_data(:, col_p_self)), 'dd', dd, 'subj', subj);
    end
end

    
%% Run MCMC
[samples, stats, structArray] = matjags( ...
datastruct, ...                               % Observed data
fullfile(pwd, model_name), ...                % File that contains model definition
init0, ...                                    % Initial values for latent variables
'doparallel' , doparallel, ...                % Parallelization flag
'nchains', nchains,...                        % Number of MCMC chains
'nburnin', nburnin,...                        % Number of burnin steps
'nsamples', nsamples, ...                     % Number of samples to extract
'thin', nthin, ...                            % Thinning parameter
'monitorparams', {'mod', 'pmod'}, ...         % List of latent variables to monitor
'savejagsoutput' , 0 , ...                    % Save command line output produced by JAGS?
'verbosity' , 2 , ...                         % 0=do not produce any output; 1=minimal text output; 2=maximum text output
'cleanup' , 1);

fprintf('\n\n');

% %% Plot posterior distributions mod for each subject
% figure;
% eps=1;
% %range=1;
% %bins=-range:eps:range;
% bins=1:5;
% 
% for subj = 1 : 45
%     
%     count=hist(samples.mod(:,:,subj),bins);
%     count=count/sum(count)/eps;
%     
%     subplot(ceil(sqrt(45)), ceil(sqrt(45)), subj),
%     %plot(bins,count,'k-');
%     bar(count);
%     xlabel('Mod','fontsize', 12);
%     ylabel('Mod. Prob. (Density)','fontsize', 12);
%     title(['P' num2str(subj)]);
%     
% end
% 
% %% Plot MCMC chains for mod for each subject
% figure;
% for subj = 1 : 45
%     subplot(ceil(sqrt(45)), ceil(sqrt(45)), subj);
%     colors = {'r', 'k', 'b', 'y', 'g'};
%     for c = 1 : nchains
%         plot(samples.mod(c,:,subj), ['o-' colors{c}]);
%         ylim([0 6]);
%     end
% end

toc
