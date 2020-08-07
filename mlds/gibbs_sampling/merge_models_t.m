function [model, r_sk, extra] = merge_models_t(models, r_sks, extras)
%MERGE_MODELS_T Summary of this function goes here
%   Detailed explanation goes here
num_runs = length(models);

for run = 1:num_runs
    nll(:,run) = extras{run}.nll;
end

if 1
    figure;
    colors = distinguishable_colors(num_runs);
    
    for i = 1:num_runs
        semilogx(nll(:,i), 'color', colors(i,:)); 
        hold on;
    end
    xlabel('Iteration number'); ylabel('NLL');
end

chosen = find_optimal_converged_chains(nll);

if isempty(chosen)
    model = [];
    r_sk = [];
    extra = [];
    return;
end

if 1
    % choose one of them
    model = models{chosen(1)};
    r_sk = r_sks{chosen(1)};
    extra = extras{chosen(1)};
else
    % merge different results
    [model, r_sk, extra] = merge_chosen(chosen, models, r_sks, extras);
end

end

function [model, r_sk, extra] = merge_chosen(chosen, models, r_sks, extras)
template_ind = chosen(1);
[D,D,K] = size(models{template_ind}.As);
template_As = reshape(models{template_ind}.As, [D*D K])';

As_all_merged = extras{template_ind}.As_all;
sigma2s_all_merged = extras{template_ind}.sigma2s_all;
nus_all_merged = extras{template_ind}.nus_all;
pis_all_merged = extras{template_ind}.pis_all;
cluster_sizes_all_merged = extras{template_ind}.cluster_sizes_all;
zs_all_merged = extras{template_ind}.zs_all;

for model_ind = 2:length(chosen)
    extra = extras{chosen(model_ind)};
    model = models{chosen(model_ind)};
    As1 = reshape(model.As, [D*D K])';
    
    % how to handle nan in As is not yet solved in the permute function
    [P, As2] = permute_endmembers(template_As, As1);
    
    nus_all = extra.nus_all(keep_ind:end,:);
    sigma2s_all = extra.sigma2s_all(keep_ind:end,:);
    pis_all = extra.pis_all(keep_ind:end,:);
    cluster_sizes_all = extra.cluster_sizes_all(keep_ind:end,:);
    zs_all = extra.zs_all(keep_ind:end,:);
    As_all = extra.As_all(:,:,:,keep_ind:end);
    
    nus_all_merged = cat(1, nus_all_merged, nus_all);
    sigma2s_all_merged = cat(1, sigma2s_all_merged, sigma2s_all * P');
    pis_all_merged = cat(1, pis_all_merged, pis_all * P');
    cluster_sizes_all_merged = cat(1, cluster_sizes_all_merged, ...
        cluster_sizes_all * P');
    
    zs_all1 = one_hot_encoding(zs_all, 'to_one_hot');

    for i = 1:size(pis_all,1)
        As1 = reshape((P * reshape(extra.As_all(:,:,:,i), [D*D K])')', [D D K]);
        As_all_merged(:,:,:,end+1) = As1;
        zs1 = squeeze(zs_all1(i,:,:)) * P';
        zs1 = one_hot_encoding(reshape(zs1, [1 size(zs1,1) size(zs1,2)]), 'from_one_hot');
        zs_all_merged(end+1,:) = zs1;
    end
end

extra.keep_ind = 1;
extra.As_all = As_all_merged;
extra.sigma2s_all = sigma2s_all_merged;
extra.nus_all = nus_all_merged;
extra.pis_all = pis_all_merged;
extra.cluster_sizes_all = cluster_sizes_all_merged;
extra.zs_all = zs_all_merged;

mus_first = models{template_ind}.mus_first;
Sigmas_first = models{template_ind}.Sigmas_first;
[model, r_sk] = calc_model(extra, mus_first, Sigmas_first);

end

function [model, r_sk] = calc_model(extra, mus_first, Sigmas_first)
keep_ind = extra.keep_ind;
r_ik = calc_prob_z(extra.zs_all(keep_ind:end,:));

model = [];
model.As = mean(extra.As_all(:,:,:,keep_ind:end), 4);
model.sigmas = mean(sqrt(extra.sigma2s_all(keep_ind:end,:)), 1);
model.nus = mean(extra.nus_all(keep_ind:end,:), 1);
model.pis = mean(extra.pis_all(keep_ind:end,:), 1);
model.mus_first = mus_first;
model.Sigmas_first = Sigmas_first;

end

