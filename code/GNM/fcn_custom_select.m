function xnew = fcn_custom_select(x,fval,ndraw,pow)

weights = fval.^-pow;

[~,points2sample_ind] = datasample(fval,ndraw,'Replace',true,'Weights',weights);

nParam = size(x,2);

xnew = zeros(ndraw,nParam);

for i = 1:ndraw
    point = x(points2sample_ind(i),:);
    nearest_per_param = point - x;
    nearest_per_param(nearest_per_param==0)= nan;
    dist2nearest_param = nanmin(nearest_per_param);
    dist2nearest_param(isnan(dist2nearest_param)) = 0;
    R = rand(nParam,1);
    xnew(i,:) = find_point_on_line(point',point' - dist2nearest_param',R);
end