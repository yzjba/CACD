function [cirparam,sd_accum]=curature_HT(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *summary: curature aided circle detection
% *input:
% I - image
% Other - choosable
% *output:
% cirparam: parameters of circles, circle center and radius in each row
% *special data needed: no
% *function needed: 
% parse_inputs,extract_curve,curvature_calcul,curvat_accum\oht_accum,
% find_circen,rad_est
% *author: Yao Zhenjie 
% *email: yaozhenjie@gmail.com
% *2010.6.23@Chinese Acadamy of Sciences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some parameters
multirad = 0.3;
sizinc = 7;

[I,th,arclen,snr,C,T_angle,sig,H,L,Endpoint,Gap_size] = parse_inputs(varargin{:});

[m,n,k] = size(I);


% edge detection
BW=edge(I,'canny');
figure; imshow(BW)

BW = double(BW);
BW = imnoise(BW,'salt & pepper',snr);
figure; imshow(BW)
BW([1:3,m-2:m],:) = 0;
BW(:,[1:3,n-2:n]) = 0;

% gradient calculation
[grdx, grdy] = gradient(single(I));
grdmag = sqrt(grdx.^2 + grdy.^2);

% extract curves
[curve,curve_start,curve_end,curve_mode,curve_num,BW_edge]=extract_curve(BW,Gap_size);  % Extract curves
figure; imshow(BW_edge)
BW_edge_imd = imdilate(BW_edge,strel('disk',1));
figure; imshow(BW_edge_imd)

% calculate curvature
curvat_im=curvature_calcul(curve,curve_start,curve_end,curve_mode,curve_num,BW,sig,Endpoint,C,T_angle); % Detect corners
curvat_R = 1./(curvat_im);
curvat_R(find(curvat_im==0)) = 0;
curvat_R = curvat_R.*(curvat_R>1 & curvat_R<1.5*max(m,n));


liern = floor(max(m,n)/sizinc);
Rrange = zeros(liern,2);
Rrange(1:3,:) = [3,15;8,20;13,33];

for i = 4:liern
    if Rrange(i-1,2)<1.5*max(m,n)
        Rrange(i,1) = Rrange(i-1,1)+sizinc;
        Rrange(i,2) = Rrange(i,1)+max(floor(0.5*Rrange(i,1)),15);
    else
        Rrange = Rrange(1:i-1,:);
        break;
    end
end
cirparam = [];
liern = size(Rrange,1);

sd_accum = zeros(m,n);

for lier=1:liern
    radrange = Rrange(lier,:);
    
    % accumulator calculation
    % curvature aided circle detection
    accum = curvat_accum(curvat_R,radrange,grdmag,grdx,grdy);
    
    if isempty(accum)
        continue;
    else
        % calculate circle center from accumulator
        s_accum = (accum-mean(accum(:)))/std(accum(:));
        sd_accum = max(sd_accum,s_accum);
        circen = find_circen(accum,th);
        if ~isempty(circen)
            % estimate radius
            cirparam = rad_est(cirparam,circen,radrange,grdx,grdy,multirad,BW_edge_imd,arclen);
        end
    end
end

return






function [curve,curve_start,curve_end,curve_mode,cur_num,BW_edge]=extract_curve(BW,Gap_size)

%   Function to extract curves from binary edge map, if the endpoint of a
%   contour is nearly connected to another endpoint, fill the gap and continue
%   the extraction. The default gap size is 1 pixles.

[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1);
cur_num=0;

while size(r,1)>0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[cur;point];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
        
    if size(cur,1)>10
        cur_num=cur_num+1;
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end

for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
        (curve_start(i,2)-curve_end(i,2))^2<=32
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1;
end
% figure(1)
% imshow(~BW_edge)
% title('Edge map')
% imwrite(~BW_edge,'edge.jpg');


function curvat_im=curvature_calcul(curve,curve_start,curve_end,curve_mode,curve_num,BW,sig,Endpoint,C,T_angle)
% caltulate the curvature

curvat_im = zeros(size(BW));

corner_num=0;
cout=[];

GaussianDieOff = .0001; 
pw = 1:30;
ssq = sig*sig;
width = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
if isempty(width)
    width = 1;  
end
t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
gau=gau/sum(gau);

for i=1:curve_num;
    x=curve{i}(:,1);
    y=curve{i}(:,2);
    W=width;
    L=size(x,1);
    if L>W
        
        % Calculate curvature
        if curve_mode(i,:)=='loop'
            x1=[x(L-W+1:L);x;x(1:W)];
            y1=[y(L-W+1:L);y;y(1:W)];
        else
            x1=[ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
            y1=[ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end
       
        xx=conv(x1,gau);
        xx=xx(W+1:L+3*W);
        yy=conv(y1,gau);
        yy=yy(W+1:L+3*W);
        Xu=[xx(2)-xx(1) ; (xx(3:L+2*W)-xx(1:L+2*W-2))/2 ; xx(L+2*W)-xx(L+2*W-1)];
        Yu=[yy(2)-yy(1) ; (yy(3:L+2*W)-yy(1:L+2*W-2))/2 ; yy(L+2*W)-yy(L+2*W-1)];
        Xuu=[Xu(2)-Xu(1) ; (Xu(3:L+2*W)-Xu(1:L+2*W-2))/2 ; Xu(L+2*W)-Xu(L+2*W-1)];
        Yuu=[Yu(2)-Yu(1) ; (Yu(3:L+2*W)-Yu(1:L+2*W-2))/2 ; Yu(L+2*W)-Yu(L+2*W-1)];
        K=abs((Xu.*Yuu-Xuu.*Yu)./((Xu.*Xu+Yu.*Yu).^1.5));
        K=conv(K,gau);

        for j = 1:L
            curvat_im(x(j),y(j)) = K(W+j);
        end
    end
end




function ang=curve_tangent(cur,center)

for i=1:2
    if i==1
        curve=cur(center:-1:1,:);
    else
        curve=cur(center:size(cur,1),:);
    end
    L=size(curve,1);
    
    if L>3
        if sum(curve(1,:)~=curve(L,:))~=0
            M=ceil(L/2);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M,1);
            y2=curve(M,2);
            x3=curve(L,1);
            y3=curve(L,2);
        else
            M1=ceil(L/3);
            M2=ceil(2*L/3);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M1,1);
            y2=curve(M1,2);
            x3=curve(M2,1);
            y3=curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2*(-y1*x2^2+y3*x2^2-y3*y1^2-y3*x1^2-y2*y3^2+x3^2*y1+y2*y1^2-y2*x3^2-y2^2*y1+y2*x1^2+y3^2*y1+y2^2*y3)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            y0 = -1/2*(x1^2*x2-x1^2*x3+y1^2*x2-y1^2*x3+x1*x3^2-x1*x2^2-x3^2*x2-y3^2*x2+x3*y2^2+x1*y3^2-x1*y2^2+x3*x2^2)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_direction=angle(complex(x0-x1,y0-y1));
            adjacent_direction=angle(complex(x2-x1,y2-y1));
            tangent_direction=sign(sin(adjacent_direction-radius_direction))*pi/2+radius_direction;
        end
    
    else % very short line
        tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    direction(i)=tangent_direction*180/pi;
end
ang=abs(direction(1)-direction(2));



function [I,th,arclen,snr,C,T_angle,sig,H,L,Endpoint,Gap_size] = parse_inputs(varargin);

error(nargchk(0,11,nargin));

Para=[6,2,30,1.5,162,3,0.2,0.01,1,1]; %Default experience value;

if nargin>=2
    I=varargin{1};
    for i=2:nargin
        if size(varargin{i},1)>0
            Para(i-1)=varargin{i};
        end
    end
end

if nargin==1
    I=varargin{1};
end

if nargin==0 | size(I,1)==0
    [fname,dire]=uigetfile('*.bmp;*.jpg;*.gif','Open the image to be detected');
    I=imread([dire,fname]);
end

th = Para(1);
arclen = Para(2);
snr = Para(3);
C = Para(4);
T_angle=Para(5);
sig=Para(6);
H=Para(7);
L=Para(8);
Endpoint=Para(9);
Gap_size=Para(10);
