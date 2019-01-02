function cirparam = circheck(cirparam)    

if isempty(cirparam)
    return;
end

cirparam(find(cirparam(:,3)==0),:) = [];  %°ë¾¶Îª0µÄÔ²È¥µô

k = 1;
while k<= size(cirparam,1)
    n = size(cirparam,1)
    ref = repmat(cirparam(k,:),n,1);
    diff = abs(cirparam-ref);
    ind = [];
    for j=k+1:n
        if diff(j,1)<4 & diff(j,2)<4 & diff(j,3)<4
            ind = [ind,j];
        end
    end
    cirparam(ind,:) = [];
    k = k+1;
end

cirparam