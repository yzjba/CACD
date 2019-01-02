function cirparam = rad_est(cirparam,circen,radrange,grdx,grdy,multirad,BW,arclen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *summary: Estimation of the Radii of Circles by local support
% *input:
% cirparam - old circle parameters
% circen - possible new circle center
% radrange - range of radius
% grdx - horizonal component of gradient
% grdy - vertical component of gradient
% multirad - parameter to check wether need to find mutiple radius
% BW - dilated edge image
% *output:
% cirparam - new circle parameters
% *special data needed: no
% *function needed:no
% *author: Yao Zhenjie
% *email: yaozhenjie@gmail.com
% *2010.6.23@Chinese Acadamy of Sciences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the estimation of the radii of circles
fltr4SgnCv = [2 1 1];
fltr4SgnCv = fltr4SgnCv / sum(fltr4SgnCv);

[m,n] = size(BW);
m = min(m,n);


for k = 1 : size(circen,1),
    % Neighborhood region of the circle for building the sgn.curve
    circen_round = round(circen(k,:));

    SCvR_I0 = circen_round(2) - radrange(2) - 1;
    if SCvR_I0 < 1,
        SCvR_I0 = 1;
    end
    SCvR_I1 = circen_round(2) + radrange(2) + 1;
    if SCvR_I1 > size(grdx,1),
        SCvR_I1 = size(grdx,1);
    end
    SCvR_J0 = circen_round(1) - radrange(2) - 1;
    if SCvR_J0 < 1,
        SCvR_J0 = 1;
    end
    SCvR_J1 = circen_round(1) + radrange(2) + 1;
    if SCvR_J1 > size(grdx,2),
        SCvR_J1 = size(grdx,2);
    end

    % Build the sgn. curve
    SgnCvMat_dx = repmat( (SCvR_J0:SCvR_J1) - circen(k,1) , ...
        [SCvR_I1 - SCvR_I0 + 1 , 1] );
    SgnCvMat_dy = repmat( (SCvR_I0:SCvR_I1)' - circen(k,2) , ...
        [1 , SCvR_J1 - SCvR_J0 + 1] );
    SgnCvMat_r = sqrt( SgnCvMat_dx .^2 + SgnCvMat_dy .^2 );
    SgnCvMat_rp1 = round(SgnCvMat_r) + 1;

    f4SgnCv = abs( ...
        double(grdx(SCvR_I0:SCvR_I1, SCvR_J0:SCvR_J1)) .* SgnCvMat_dx + ...
        double(grdy(SCvR_I0:SCvR_I1, SCvR_J0:SCvR_J1)) .* SgnCvMat_dy ...
        ) ./ SgnCvMat_r;
    SgnCv = accumarray( SgnCvMat_rp1(:) , f4SgnCv(:) );

    
    ee = double(BW(SCvR_I0:SCvR_I1,SCvR_J0:SCvR_J1));
    SgnCv_Cnt = accumarray( SgnCvMat_rp1(:) , ee(:) );
    SgnCv_Cnt = SgnCv_Cnt + (SgnCv_Cnt == 0);

    % Suppress the undesired entries in the sgn. curve
    % -- Radii that correspond to short arcs
    rran = [0:(numel(SgnCv_Cnt)-1)]';
    
    % minimal circumference to be detected as a circle
%     th=0.7+log(5)./log(rran);
%     th=min(th,1.8);
%     th(find(th<0.9)) = 0.6;

%     bb = round(numel(SgnCv_Cnt)/10);
%     ee = round(numel(SgnCv_Cnt)*9/10);
%     th = 1.2+0.7*cos(pi*(rran-bb)/(ee-bb));
%     th(1:bb) = 1.9;
%     th(ee:end) = 0.5;
    
    rs = numel(SgnCv_Cnt);
%     rs = max(m,n);
    th = -1.4/(rs)^2*(rran.^2)+arclen;
    th(1:6) = 2;
    th(6:13) = 1.8;
    th = max(1,th);
    
    SgnCv = SgnCv .* ( SgnCv_Cnt >= (pi * th.*[0:(numel(SgnCv_Cnt)-1)]') );
    % -- Radii that are out of the given range
    SgnCv( 1 : (round(radrange(1))+1) ) = 0;
    SgnCv( (round(radrange(2))+1) : end ) = 0;

    % Get rid of the zero radius entry in the array
    SgnCv = SgnCv(2:end);
    % Smooth the sgn. curve
    SgnCv = filtfilt( fltr4SgnCv , [1] , SgnCv );

    % Get the maximum value in the sgn. curve
    SgnCv_max = max(SgnCv);
    if SgnCv_max <= 0,
        continue;
    end

    % Find the local maxima in sgn. curve by 1st order derivatives
    % -- Mark the ascending edges in the sgn. curve as 1s and
    % -- descending edges as 0s
    SgnCv_AscEdg = ( SgnCv(2:end) - SgnCv(1:(end-1)) ) > 0;
    % -- Mark the transition (ascending to descending) regions
    SgnCv_LMmask = [ 0; 0; SgnCv_AscEdg(1:(end-2)) ] & (~SgnCv_AscEdg);
    SgnCv_LMmask = SgnCv_LMmask & [ SgnCv_LMmask(2:end) ; 0 ];

    % Incorporate the minimum value requirement
    SgnCv_LMmask = SgnCv_LMmask & ...
        ( SgnCv(1:(end-1)) >= (multirad * SgnCv_max) );
    % Get the positions of the peaks
    SgnCv_LMPos = sort( find(SgnCv_LMmask) );

    % Save the detected radii
    if isempty(SgnCv_LMPos),
        continue;
    else
        cirparam = [cirparam;circen(k,:),SgnCv_LMPos(end)];
        for i_radii = (length(SgnCv_LMPos) - 1) : -1 : 1,
            cirparam = [cirparam;circen(k,:),SgnCv_LMPos(i_radii)];
        end
    end
end