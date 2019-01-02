function circen=find_circen(accum,th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *summary: find circle center by outlier analysis
% *input:
% accum - accumulator arry
% radrange - range of radius
% *output:
% circen: possible circle center
% *special data needed: no
% *function needed:no
% *author: Yao Zhenjie
% *email: yaozhenjie@gmail.com
% *2010.6.23@Chinese Acadamy of Sciences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fltr4LM_R = 8;
prm_fltrLM_s = 1.35;
prm_fltrLM_r = ceil(fltr4LM_R * 0.6 );
prm_fltrLM_npix = max([ 6, ceil((fltr4LM_R/2)^1.8)]);

circen = [];
[m,n] = size(accum);

% Thresholding of 'accum' by a lower bound
prm_LM_LoBnd = max(accum(:)) * 0.2;

% check out whether there is a circle
if max(accum(:))<mean(accum(:))+th*std(accum(:));
    return;
end

% Thresholding of 'accum' by a lower bound
accumaoi_LBMask = ...
    ( accum > max(prm_LM_LoBnd,mean(accum(:))+th*std(accum(:))) );

candLM_mask = accumaoi_LBMask;%( candLM > mean(aoi_s)+3*std(aoi_s));

% Clear the margins of 'candLM_mask'
candLM_mask([1:fltr4LM_R, (end-fltr4LM_R+1):end], :) = 0;
candLM_mask(:, [1:fltr4LM_R, (end-fltr4LM_R+1):end]) = 0;


% Group the local maxima candidates by adjacency, compute the
% centroid position for each group and take that as the center
% of one circle detected
[candLM_label, candLM_nRgn] = bwlabel( candLM_mask, 8 );

for ilabel = 1 : candLM_nRgn,
    % Indices (to current AOI) of the pixels in the group
    candgrp_masklin = find( candLM_label == ilabel );
    [candgrp_IdxI, candgrp_IdxJ] = ...
        ind2sub( size(candLM_label) , candgrp_masklin );

    % Indices (to 'accum') of the pixels in the group
    candgrp_idx2acm = ...
        sub2ind( size(accum) , candgrp_IdxI , candgrp_IdxJ );

    % Minimum number of qualified pixels in the group
    if sum(accumaoi_LBMask(candgrp_masklin)) < prm_fltrLM_npix,
        continue;
    end

    % Compute the centroid position
    candgrp_acmsum = sum( accum(candgrp_idx2acm) );
    cc_x = sum( candgrp_IdxJ .* accum(candgrp_idx2acm) ) / ...
        candgrp_acmsum;
    cc_y = sum( candgrp_IdxI .* accum(candgrp_idx2acm) ) / ...
        candgrp_acmsum;
    circen = [circen; cc_x, cc_y];
end

