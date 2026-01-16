% first compute sliding windows to try decoding from
spacing = 50;start_sliding = 5000;length_epoch = 300;
last_sliding = 6000 - length_epoch;num_epochs = (last_sliding - start_sliding)/spacing + 1;
sliding_epochs= zeros(num_epochs,2);
sliding_epochs(1,:) = [start_sliding (start_sliding +length_epoch-1)];

for i=2:num_epochs
    sliding_epochs(i,:) = [(sliding_epochs(i-1,1)+spacing) (sliding_epochs(i-1,2)+ spacing)];
end

num_subjects  = 6;
subject_array = 1:6;

% Initialize storage for results
all_subjects_accuracies               = cell(num_subjects, num_epochs);  % each cell: [1 x 1] mean LOOCV accuracy
all_subjects_mean_accuracies          = zeros(num_subjects, num_epochs); % same as above
all_subjects_shuffled_accuracies      = cell(num_subjects, num_epochs);  % mean shuffled LOOCV accuracy
all_subjects_mean_shuffled_accuracies = zeros(num_subjects, num_epochs);
all_subjects_selected_predictors      = cell(num_subjects, num_epochs);  % stable electrodes (global indices)
all_subjects_selection_counts         = cell(num_subjects, num_epochs);  % selection frequency per electrode

for j = 1:num_subjects
    n = subject_array(j);
    fprintf('Processing Subject %d...\n', n);

    % Trial-level filter (same across epochs)
    filter = subject(n).decision < 2 & ...
             subject(n).electrode(1).trigger(1).good_trials == 1 & ...
             subject(n).electrode(1).trigger(3).good_trials == 1;

    % --- Use ONLY OFC electrodes (olf_coor > 0) ---
    num_electrodes = length(subject(n).electrode);
    olf_mask       = arrayfun(@(el) el.olf_coor > 0, subject(n).electrode);
    olf_idx        = find(olf_mask);            % indices into subject(n).electrode
    num_olf        = numel(olf_idx);

    for i = 1:num_epochs

        % ----- Build predictors & response for this subject and epoch -----
        valid_trials = filter;
        n_trials     = sum(valid_trials);

        predictors_epoch = zeros(n_trials, num_olf);

        for k = 1:num_olf
            e = olf_idx(k);    % actual electrode index

            % Signal in current epoch
            signal = sum(subject(n).electrode(e).trigger(3).high_gamma_mat( ...
                valid_trials, sliding_epochs(i,1):sliding_epochs(i,2)), 2) / length_epoch;

            predictors_epoch(:, k) = signal;                             
        end

        response_epoch = subject(n).decision(valid_trials);   % 0/1

        % ----- Leave-One-Out CV -----
        % Each trial is held out once
        cv = cvpartition(n_trials, 'LeaveOut');
        K  = cv.NumTestSets;

        fold_acc      = zeros(K, 1);
        fold_acc_null = zeros(K, 1);

        % how often each OFC predictor is selected across all LOOCV folds
        selection_counts_olf = zeros(1, num_olf);

        for kfold = 1:K
            train_idx = training(cv, kfold);
            test_idx  = test(cv, kfold);

            X_train = predictors_epoch(train_idx, :);  % [nTrain x num_olf]
            y_train = response_epoch(train_idx);

            X_test  = predictors_epoch(test_idx, :);   % [1 x num_olf]
            y_test  = response_epoch(test_idx);        % scalar

            % ===== REAL MODEL (LASSO logistic) =====
            [B, FitInfo] = lassoglm(X_train, y_train, 'binomial', 'Alpha',0.4,'CV', 5);

            % choose lambda (either IndexMinDeviance or Index1SE)
            idxLambda = FitInfo.IndexMinDeviance;
            % idxLambda = FitInfo.Index1SE;   % try this if you want more regularization

            beta0 = FitInfo.Intercept(idxLambda);   % scalar
            beta  = B(:, idxLambda);                % [num_olf x 1]

            % Nonzero coefficients = selected predictors
            inModel = beta ~= 0;                    % logical [num_olf x 1]
            selection_counts_olf = selection_counts_olf + inModel';  % accumulate across folds

            % Predict on held-out trial
            eta   = beta0 + X_test * beta;
            p_hat = 1 ./ (1 + exp(-eta));
            y_pred = p_hat > 0.5;

            fold_acc(kfold) = (y_pred == y_test);

            % ===== SHUFFLED NULL MODEL (same LASSO, shuffled labels) =====
            y_train_shuff = y_train(randperm(length(y_train)));

            [B_null, FitInfo_null] = lassoglm(X_train, y_train_shuff, 'binomial', 'Alpha',0.3,'CV', 5);

            idxLambda_null = FitInfo_null.IndexMinDeviance;
            % idxLambda_null = FitInfo_null.Index1SE;   % keep parallel if you swap above

            beta0_null = FitInfo_null.Intercept(idxLambda_null);
            beta_null  = B_null(:, idxLambda_null);

            eta_null   = beta0_null + X_test * beta_null;
            p_hat_null = 1 ./ (1 + exp(-eta_null));
            y_pred_null = p_hat_null > 0.5;

            fold_acc_null(kfold) = (y_pred_null == y_test);
        end

        % LOOCV accuracies (mean over folds)
        repeat_accuracies          = mean(fold_acc);      % scalar
        shuffled_repeat_accuracies = mean(fold_acc_null); % scalar

        % Store per-epoch accuracy and mean accuracy (same here)
        all_subjects_accuracies{j, i}              = repeat_accuracies;
        all_subjects_mean_accuracies(j, i)         = repeat_accuracies;

        all_subjects_shuffled_accuracies{j, i}     = shuffled_repeat_accuracies;
        all_subjects_mean_shuffled_accuracies(j,i) = shuffled_repeat_accuracies;

        % ----- Selection frequencies and "stable" predictors -----
        total_models = K;  % total LOOCV folds
        selection_freq_olf = selection_counts_olf / total_models;   % 0–1 frequency

        % Map frequencies back to FULL electrode indexing
        selection_counts_full = zeros(1, num_electrodes);
        selection_counts_full(olf_idx) = selection_freq_olf;
        all_subjects_selection_counts{j,i} = selection_counts_full;

        % "Stable" predictors = selected in ≥50% of LOOCV folds
        stable_olf_idx_local     = find(selection_freq_olf >= 0.5);   % in OFC-only index
        stable_predictors_global = olf_idx(stable_olf_idx_local);     % in original electrode index
        all_subjects_selected_predictors{j, i} = stable_predictors_global;
    end
end

% save accuracies

save('accuracies.mat','all_subjects_selected_predictors','all_subjects_accuracies',...
    'all_subjects_mean_accuracies','all_subjects_shuffled_accuracies','all_subjects_mean_shuffled_accuracies','sliding_epochs');
