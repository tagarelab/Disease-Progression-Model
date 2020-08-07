function [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
    LGPmask, RGPmask, LTHmask, RTHmask] = load_masks()
%RETRIEVE_META_DATA This function will be used throughout the data creation
% step so changing the implementation here will affect all the steps.

if 0 % load LR flipped and registered masks
    load('../masks/flipped_registered_masks.mat');
elseif 1 % load dilated masks from the original masks
    [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask, ...
        WSmask, LGPmask, RGPmask, LTHmask, RTHmask] = ...
        load_dilated_masks(1);
elseif 0 % load original masks
    [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
        LGPmask, RGPmask, LTHmask, RTHmask] = load_masks_original();
else % load Hemant's masks
    load('..\masks\occMask');
    
    load('..\masks\leftStriatum');
    LSmask=mask;
    load('..\masks\rightStriatum');
    RSmask=mask;
    load('..\masks\leftCaudate');
    LCmask=mask;
    load('..\masks\rightCaudate');
    RCmask=mask;
    load('..\masks\leftPutamen');
    LPmask=mask;
    load('..\masks\rightPutamen');
    RPmask=mask;
    clear mask
end

end

function show_mask(mask)
figure;
ind = 1; 
for i = 1:5
    for j = 1:6
        if ind <= size(mask,3)
            subplot(5,6,ind); 
            imshow(mask(:,:,ind),[]); 
            ind = ind+1; 
        end
    end 
end

end