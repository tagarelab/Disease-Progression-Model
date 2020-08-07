function [features_full, labels_fea, raw_full, labels_raw, ts] = load_updrs(dataset_name)
%LOAD_UPDRS Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    dataset_name = 'processed_images';
end

% load data
dataset = load_dataset(dataset_name, [], [], []);

train_ids = dataset.train_ids;       

% load updrs
if 0 % load data from Johnathan (2019)
    [raw_full1, labels_raw] = load_data_by_combining_two_sources(train_ids);
    load('../updrs/full_updrs_data.mat');
    [raw_full, ts] = load_data_directly(train_ids, data);
    assert(nansum(nansum(nansum(abs(raw_full(:,2:end,:) - raw_full1)))) == 0);
else % load data (csv file) from Johnathan (2020.2)
    if 0
        [data, labels] = load_data_from_csv('../updrs/updrs2020_sorted.csv');
        save('../updrs/updrs2020_sorted.mat','data','labels');
    end
    load('../updrs/updrs2020_sorted.mat');
    [raw_full, ts] = load_data_directly(train_ids, data);
    labels_raw = labels;
end

[features_full, labels_fea] = create_features(raw_full);

% convert date to datenum using the 15th of the month
for i = 1:numel(ts)
    if ~isnan(ts(i))
        month = floor(ts(i) / 10000);
        year = ts(i) - month * 10000;
        dv = [year,month,15];
        ts(i) = datenum(dv);
    end
end

% check time is monotonically increasing
all_increasing = zeros(size(ts,1), 1);
for i = 1:size(ts,1)
    ts_i = ts(i,:);
    ts_i(isnan(ts_i)) = [];
    all_increasing(i) = all(ts_i(:,2:end) - ts_i(:,1:end-1) >= 0, 2);
end
assert(all(all_increasing));


% show_data(features_full_1, labels, raw_full_1);

end

function [data, labels] = load_data_from_csv(csvfile)
M = read_mixed_csv(csvfile, ',');
labels = M(1,4:end);

patient_ids_str = unique(M(2:end,2));
patient_ids_num = cellfun(@str2num, patient_ids_str);
n = length(patient_ids_str);
d = size(M,2) - 1;
t = 17; % ts max is 17
data = nan(n, t, d);

ts = [];
for i = 1:n
    inds = find(strcmp(M(1:end,2), patient_ids_str{i}));
    data_i = convert_str_to_num(M(inds, 2:end));
    ts(i) = size(data_i, 1);
    data(i,1:ts(i),:) = data_i;
end

end

function mat = convert_str_to_num(str_cellarray)
mat = nan(size(str_cellarray));

for i = 1:numel(str_cellarray)
    if ~isempty(str_cellarray{i})
        mat(i) = str2num(str_cellarray{i});
    end
end
    
end

function [raw_full, ts] = load_data_directly(train_ids, data)
ts = [];
for i = 1:length(train_ids)
    ind = find(data(:,1,1) == str2num(train_ids{i}));
    data_i = data(ind,:,:);
    data_i(1, data_i(1,:,1) == 0, :) = nan;
    raw_full(i,:,:) = data_i(1,:,3:end);
    ts(i,:) = data_i(1,:,2); 
end
end

function [raw_full, labels_raw] = load_data_by_combining_two_sources(train_ids)
% load('../updrs/feature_patient_data_corrected.mat');
% features_first = features;
load('../updrs/raw_patient_data.mat');
raw_first = data;
% load('../updrs/feature_data_full_timeseries.mat');
% features_full_ori = features;
load('../updrs/raw_data_full_timeseries.mat');
raw_full_ori = features;


[raw_full, labels_raw] = find_corresponding_raw_data( ...
    train_ids, raw_first, raw_full_ori);
end

function [features_full, labels_fea] = create_features(raw_full)
if 0
    labels_fea = {'Speech','Rigidity','Bradykinesia','Posture','Gait','Tremor'};

    speech_inds = [1];
    rigidity_inds = [3,4,5,6,7];
    bradykinesia_inds = [2,8,9,10,11,12,13,14,15,16,17,23];
    posture_inds = [18,21,22];
    gait_inds = [19,20];
    tremor_inds = [24,25,26,27,28,29,30,31,32,33];

    inds = {speech_inds, rigidity_inds, bradykinesia_inds, posture_inds, ...
        gait_inds, tremor_inds};
else
    labels_fea = {'Rigidity L', 'Rigidity R', 'Bradykinesia L', ...
        'Bradykinesia R', 'Tremor L', 'Tremor R'};
    rigidity_R_inds = [4,6];
    rigidity_L_inds = [5,7];
    bradykinesia_R_inds = [8,10,12,14,16];
    bradykinesia_L_inds = [9,11,13,15,17];
    tremor_R_inds = [24,26,28,30];
    tremor_L_inds = [25,27,29,31];
    
    inds = {rigidity_R_inds, rigidity_L_inds, bradykinesia_R_inds, ...
        bradykinesia_L_inds, tremor_R_inds, tremor_L_inds};
end

for k = 1:length(inds)
    features_full(:,:,k) = mean(raw_full(:,:,inds{k}), 3);
end

end

function [raw_full, labels_raw] = find_corresponding_raw_data( ...
    train_ids, raw_first, raw_full_ori)
labels_raw = {'Speech', 'Facial expression', 'Rigidity neck', ... % 1 - 3
    'Rigidity RUE', 'Rigidity LUE', 'Rigidity RLE', 'Rigidity LLE', ... % 4 - 7 
    'Finger tapping R', 'Finger tapping L', 'Hand movements R', ... % 8 - 10
    'Hand movements L', 'Pronation-supination of hands R', ... % 11 - 12
    'Pronation-supination of hands L', 'Toe tapping R', 'Toe tapping L', ... % 13 - 15
    'Leg agility R', 'Leg agility L', 'Arising from chair', 'Gait', ... % 16 - 19
    'Freezing of gait', 'Postural stability', 'Posture', ... % 20 - 22
    'Body bradykinesia', 'Postural tremor of hands R', ... % 23 - 24
    'Postural tremor of hands L', 'Kinetic tremor of hands R', ... % 25 - 26
    'Kinetic tremor of hands L', 'Rest tremor amplitude RUE', ... % 27 - 28
    'Rest tremor amplitude LUE', 'Rest tremor amplitude RLE', ... % 29 - 30
    'Rest tremor amplitude LLE', 'Rest tremor amplitude LipJaw', ... % 31 - 32
    'Constancy of rest tremor', 'Hoehn and yahr stage'}; % 33 - 34

% pick the corresponding patient ids 
for i = 1:length(train_ids)
    id = str2num(train_ids{i});
    raw = find_full_series_data(raw_first, raw_full_ori, id);
    raw_full(i,:,:) = raw;
end

% remove the last score (Hoehn and Yahr stage) and the 5 following scores
raw_full(:,:,34:end) = [];
labels_raw(end) = [];

% set missing data as nan and normal as 0
raw_full(raw_full == 0) = nan;
raw_full = raw_full - 1;

end

function [raw] = find_full_series_data(raw_first, raw_full, id)
x = raw_first(raw_first(:,1) == id, 2:end);
inds = all(squeeze(raw_full(:,1,:)) == repmat(x, [size(raw_full,1),1]), 2);
ks = find(inds);
if length(ks) > 1 % more than 1 matched raw data
    error('Multiple patients have same baseline raw data');
else
    raw = raw_full(ks,:,:);
end
end

function show_data(features_full_1, labels, raw_full_1)
features_full_1(features_full_1 == 0) = nan;
figure,
for k = 1:size(features_full_1,3)
    subplot(2,3,k);
    plot(features_full_1(:,:,k)');
    xlabel(labels{k});
end
end
