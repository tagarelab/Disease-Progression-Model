function zs1 = one_hot_encoding(zs, mode)
%ONE_SHOT_ENCODING Summary of this function goes here
%   Detailed explanation goes here
if strcmp(mode, 'to_one_hot')
    [num_iter, N] = size(zs);
    K = max(zs(:));
    
    zs1 = zeros(num_iter,N,K);

    for k = 1:K
        zs1(:,:,k) = double(zs == k);
    end

elseif strcmp(mode, 'from_one_hot')
    [num_iter,N,K] = size(zs);
    zs1 = zeros(num_iter, N);
    
    for k = 1:K
        zs1(zs(:,:,k) == 1) = k;
    end
end

end

