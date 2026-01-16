%% Efficient LOOCV posterior computation with parfor (fit once per trial)

clear trial_epochs;
spacing = 5; start_sliding = 1; length_trial_epoch = 10;
last_sliding = 6001 - length_trial_epoch;
num_epochs1 = (last_sliding - start_sliding)/spacing + 1;

trial_epochs = zeros(num_epochs1,2);
trial_epochs(1,:) = [start_sliding (start_sliding + length_trial_epoch - 1)];
for i = 2:num_epochs1
    trial_epochs(i,:) = trial_epochs(i-1,:) + spacing;
end

num_shuffles = 1;

for n = 1:6

    best_epoch = find(all_subjects_mean_accuracies(n,:) == max(all_subjects_mean_accuracies(n,:)));
    best_epoch = max(best_epoch);

    filter = find(subject(n).decision < 2 & ...
                  subject(n).electrode(1).trigger(3).good_trials == 1);

    response_variable = subject(n).decision(filter);

    electrodes = predictors_global{n,best_epoch};
    electrode_indices = electrodes(all_subjects_predictors_stability{n,best_epoch} > 0.9);

    num_trials = numel(filter);
    num_feats  = numel(electrode_indices);

    % ----------------------------------------------------------
    % Build TRAINING predictors once (fixed across epochs)
    % ----------------------------------------------------------
    predictors = zeros(num_trials, num_feats);

    tTrain1 = sliding_epochs(best_epoch,1);
    tTrain2 = sliding_epochs(best_epoch,2);

    for f = 1:num_feats
        tmp = subject(n).electrode(electrode_indices(f)).trigger(3).high_gamma_mat(filter, tTrain1:tTrain2);
        if isa(tmp,'gpuArray'), tmp = gather(tmp); end
        predictors(:,f) = mean(tmp, 2);
    end

    % Ensure CPU for GLM
    if isa(predictors,'gpuArray'), predictors = gather(predictors); end
    if isa(response_variable,'gpuArray'), response_variable = gather(response_variable); end

    % Full model (kept for saving / reference)
    glm_model = fitglm(predictors, response_variable, 'Distribution','binomial');

    % Outputs
    posterior_probabilities = NaN(num_trials, num_epochs1);
    posterior_probabilities_shuffled = cell(1, num_shuffles);
    for j = 1:num_shuffles
        posterior_probabilities_shuffled{j} = NaN(num_trials, num_epochs1);
    end

    % Pre-store RT (CPU) to avoid repeated struct access in parfor
    rt_all = subject(n).rt(filter);
    if isa(rt_all,'gpuArray'), rt_all = gather(rt_all); end

    % Also cache trial indices vector
    idx_all = (1:num_trials)';

    % ----------------------------------------------------------
    % PARFOR over trials: fit LOO once, then predict all epochs
    % ----------------------------------------------------------
    parfor x = 1:num_trials

        % local outputs for this trial
        post_x = NaN(1, num_epochs1);
        post_x_shuf = cell(1, num_shuffles);
        for j = 1:num_shuffles
            post_x_shuf{j} = NaN(1, num_epochs1);
        end

        % LOO training indices
        train_idx = idx_all(idx_all ~= x);

        % Fit real LOO model once
        glm_LOO = fitglm(predictors(train_idx,:), response_variable(train_idx), ...
                         'Distribution','binomial');

        % Fit shuffled LOO model(s) once
        glm_LOO_shuf = cell(1, num_shuffles);
        for j = 1:num_shuffles
            ysh = response_variable(train_idx);
            ysh = ysh(randperm(numel(ysh)));
            glm_LOO_shuf{j} = fitglm(predictors(train_idx,:), ysh, ...
                                     'Distribution','binomial');
        end

        % Compute epoch-wise posteriors for trial x
        rt_x = rt_all(x);

        for i = 1:num_epochs1

            if (6000 - trial_epochs(i,1)) > rt_x
                continue;
            end

            t1 = trial_epochs(i,1);
            t2 = trial_epochs(i,2);

            % Build predictors for this trial + epoch
            new_predictors = zeros(1, num_feats);

            for b = 1:num_feats
                sig = subject(n).electrode(electrode_indices(b)).trigger(3).high_gamma_mat(filter(x), t1:t2);
                if isa(sig,'gpuArray'), sig = gather(sig); end
                new_predictors(b) = mean(sig, 2);
            end

            % Predict (real)
            post_x(i) = predict(glm_LOO, new_predictors);

            % Predict (shuffled)
            for j = 1:num_shuffles
                post_x_shuf{j}(i) = predict(glm_LOO_shuf{j}, new_predictors);
            end
        end

        % write back to sliced outputs
        posterior_probabilities(x,:) = post_x;
        for j = 1:num_shuffles
            posterior_probabilities_shuffled{j}(x,:) = post_x_shuf{j};
        end
    end

    savename = strcat(['posteriors_' num2str(n)]);
    save(savename,'posterior_probabilities','posterior_probabilities_shuffled','glm_model');

    savename = strcat(['chosen_indices_' num2str(n)]);
    save(savename,'electrode_indices');
end
