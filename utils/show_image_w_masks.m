function show_image_w_masks(img_ori,list_oper_slices,list_masks,...
    list_colors,fh_image_w_mask)
% Example:
%   list_oper_slices{1} = [zLimL:zLimU];
%   list_masks{1} = cat(4,LCmask+RCmask,LPmask+RPmask,occMask);
%   list_colors{1} = [0.2,0,0;0,0.2,0;0,0,0.2];
%   list_oper_slices{2} = [1:91];
%   list_masks{2} = 0;
%   list_colors{2} = [0.1,0.1,0];


[rows,cols,slices] = size(img_ori);

img_color = reshape(img_ori,[rows,cols,1,slices]);
img_color = repmat(img_color, [1,1,3,1]);

for i = 1:length(list_oper_slices)
    oper_slices = list_oper_slices{i};
    masks = list_masks{i};
    colors = list_colors{i};
    
    img_per_oper = img_color(:,:,:,oper_slices);

    for j = 1:size(masks,4)
        img_per_oper(:,:,1,:) = squeeze(img_per_oper(:,:,1,:)) + ...
            masks(:,:,:,j)*colors(j,1);
        img_per_oper(:,:,2,:) = squeeze(img_per_oper(:,:,2,:)) + ...
            masks(:,:,:,j)*colors(j,2);
        img_per_oper(:,:,3,:) = squeeze(img_per_oper(:,:,3,:)) + ...
            masks(:,:,:,j)*colors(j,3);
    end

    img_color(:,:,:,oper_slices) = img_per_oper;
end

ax_h = findobj(fh_image_w_mask,'type','axes');
if isempty(ax_h)
    figure(fh_image_w_mask);
    montage(img_color);
else
    montage(img_color,'parent',ax_h);
end

end

