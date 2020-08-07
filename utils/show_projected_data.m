function show_projected_data(train_data, train_ts, Vs, class_inds)
%SHOW_PROJECTED_DATA Summary of this function goes here
%   Detailed explanation goes here
if isempty(Vs)
    K = max(class_inds);
    for k = 1:K
        Vs{k} = eye(size(train_data{1}, 1));
    end
end

K = length(Vs);
M = size(Vs{1},2);

figure;
for k = 1:K
    train_data_k = train_data(class_inds == k);
    train_ts_k = train_ts(class_inds == k);

    projs = {};
    for n = 1:length(train_data_k)
        for j = 1:M
            v1 = Vs{k}(:,j);
            if v1(1) < 0, v1 = -v1; end

            proj1 = v1' * train_data_k{n};
            t = train_ts_k{n}/365;

            subplot(M, K, k + (j-1)*K);
            plot(t, proj1); hold on
            xlabel({'Year',['Cluster ',num2str(k)]});
            ylabel(['Projected to v_',num2str(j)]);

            projs{n}(j,:) = proj1;
        end
    end

    proj_data(class_inds == k) = projs;
    
end
% subject_ids = train_ids;

% save('eigenvector_projected_data.mat','proj_data', ...
%     'proj_ts','subject_ids','first_dates');
% saveR('eigenvector_projected_data.Rdata','proj_data_1','proj_data_4', ...
%     'proj_data_ts','subject_ids','first_dates');
end

