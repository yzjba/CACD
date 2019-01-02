function accum = curvat_accum(curvat_R,radrange,grdmag,grdx,grdy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *summary: calculate the accumulator with gradient direction and curature heuristic
% *input:
% BW - edge image
% radrange - range of radius
% grdmag - magitude of gradient
% grdx - horizonal component of gradient
% grdy - vertical component of gradient
% *output:
% accum: accumulator of the image
% *special data needed: no
% *function needed:no
% *author: Yao Zhenjie
% *email: yaozhenjie@gmail.com
% *2010.6.23@Chinese Acadamy of Sciences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(curvat_R);
% only qualified curature points stays
curvat_R_l = curvat_R.*(curvat_R>radrange(1) & curvat_R<radrange(2));

numofedge_lier = length(find(curvat_R_l));
if numofedge_lier<0.5*pi*radrange(1)
    % if no enough points to form a circle
    accum = [];
    return;
else
    % fill the point with curarure radius and direction heuristic
    accum = zeros(m,n);
%     figure; imshow(curvat_R_l);
%     close;
    for i = 1:m
        for j = 1:n
            if curvat_R_l(i,j)~=0
                % the move vector
                dx = curvat_R_l(i,j)*grdy(i,j)/grdmag(i,j);
                dy = curvat_R_l(i,j)*grdx(i,j)/grdmag(i,j);
                
                % direction 1
                sup_cen = round([i+dx,j+dy]);
                if sup_cen(1)>0 && sup_cen(1)<=m && sup_cen(2)>0 && sup_cen(2)<=n
                    accum(sup_cen(1),sup_cen(2)) = accum(sup_cen(1),sup_cen(2))+1;%grdmag(i,j);
                end
                % direction 2
                sup_cen = round([i-dx,j-dy]);
                if sup_cen(1)>0 && sup_cen(1)<=m && sup_cen(2)>0 && sup_cen(2)<=n
                    accum(sup_cen(1),sup_cen(2)) = accum(sup_cen(1),sup_cen(2))+1;%grdmag(i,j);
                end
            end
        end
    end
    
    % smooth the accumulator, make each point support a region not a lonely point
    supR = 5;%max(min(50,0.4*radrange(1)),5);
    sup = (fspecial('disk',supR)>0);
    accum = filter2(double(sup),accum,'same');
%     figure; surf(accum, 'EdgeColor', 'none'); axis ij;
end