function [occMask,LSmask,RSmask,LCmask1,RCmask1,LPmask1,RPmask1, ...
    WSmask, LGPmask1, RGPmask1, LTHmask1, RTHmask1] = ...
    load_dilated_masks(dilate_size)
%LOAD_DILATED_MASKS dilate_size negative means erosion
if nargin < 1
    dilate_size = 1;
end

if 1
    [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
        LGPmask, RGPmask, LTHmask, RTHmask] = load_masks_original();
else
    load('../masks/flipped_registered_masks.mat');
end

masks_to_dilate = {LCmask, RCmask, LPmask, RPmask, LGPmask, RGPmask, ...
    LTHmask, RTHmask};

% calculate the mask of all the regions
mask_all_regions = zeros(size(LCmask));
for i = 1:length(masks_to_dilate)
    mask_all_regions = mask_all_regions + double(masks_to_dilate{i} > 0.5);
end

% calculate the positions of pixels outside the masks
k = find(mask_all_regions == 0);
[r,c,s] = ind2sub(size(mask_all_regions), k);
non_mask_pixel_pos = [c,r,s];

% calculate the positions of pixels inside the masks
for i = 1:length(masks_to_dilate)
    k = find(masks_to_dilate{i} > 0.5);
    [r,c,s] = ind2sub(size(masks_to_dilate{i}), k);
    mask_pixel_poses{i} = [c,r,s];
end

% calculate the distance of pixels outside the masks to the pixels inside
for i = 1:length(mask_pixel_poses)
    d = pdist2(mask_pixel_poses{i}, non_mask_pixel_pos, 'Euclidean', 'smallest', 1);
    D(:,i) = d';
end

[M,ind] = min(D,[],2);

% separate the whole brain pixels into different regions
for i = 1:length(mask_pixel_poses)
    pos = non_mask_pixel_pos(ind == i,:);
    k = sub2ind(size(mask_all_regions), pos(:,2), pos(:,1), pos(:,3));
    mask_tmp = zeros(size(mask_all_regions));
    mask_tmp(k) = 1;
    regions_for_masks{i} = double(mask_tmp > 0.5 | masks_to_dilate{i} > 0.5);
end

if 0
    fh = figure;
    hold on;
    plot_mask(regions_for_masks{1}, [1 1 0], 0);
    plot_mask(regions_for_masks{2}, [0.8 0.8 0], 0);
    plot_mask(regions_for_masks{3}, [0 1 0], 0);
    plot_mask(regions_for_masks{4}, [0 0.8 0], 1);
end

%% dilate the masks inside the regions for each mask
for i = 1:length(masks_to_dilate)
    mask_i = dilate_mask(masks_to_dilate{i}, dilate_size);
    masks_dilated{i} = mask_i .* regions_for_masks{i};
end

masks_to_dilate = {LCmask, RCmask, LPmask, RPmask, LGPmask, RGPmask, ...
    LTHmask, RTHmask};

LCmask1 = masks_dilated{1};
RCmask1 = masks_dilated{2};
LPmask1 = masks_dilated{3};
RPmask1 = masks_dilated{4};
LGPmask1 = masks_dilated{5};
RGPmask1 = masks_dilated{6};
LTHmask1 = masks_dilated{7};
RTHmask1 = masks_dilated{8};


if 0
    fh = figure;
    hold on;
    plot_mask(LCmask1, [1 1 0], 0);
    plot_mask(RCmask1, [0.8 0.8 0], 0);
    plot_mask(LPmask1, [0 1 0], 0);
    plot_mask(RPmask1, [0 0.8 0], 1);
end

% figure, imshow(LCmask(:,:,20));
% figure, imshow(LCmask1(:,:,20));

end
