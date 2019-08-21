function out=distortion_3(in)
%% 一些系数与图像初始化
    A =[ 895.8176,0,606.7259;0,898.4794,447.4418;0,0,1];
    D = [-0.0550,0.0574,0.0011,0.0108,0];
    fx = A(1,1);
    fy = A(2,2);
    cx = A(1,3);
    cy = A(2,3);
    k1 = D(1);
    k2 = D(2);
    p1 = D(3);
    p2 = D(4);
    tmp1 = im2double(in);
    tmp2 = rgb2gray(tmp1);
    tmp3 = zeros(size(tmp2));
%% 计算2D点对应的3D点
    [v u] = find(~isnan(tmp3));
    DDDpoints= inv(A)*[u v ones(length(u),1)]';
    r2 = DDDpoints(1,:).^2+DDDpoints(2,:).^2;
    x = DDDpoints(1,:);
    y = DDDpoints(2,:);

%% 按照公式去畸变
    x = x.*(1+k1*r2 + k2*r2.^2) + 2*p1.*x.*y + p2*(r2 + 2*x.^2);
    y = y.*(1+k1*r2 + k2*r2.^2) + 2*p2.*x.*y + p1*(r2 + 2*y.^2);
    ud = reshape(fx*x + cx,size(tmp3));
    vd = reshape(fy*y + cy,size(tmp3));
%% 插值重建原始图像
    tmp3(:,:,1) = interp2(tmp1(:,:,1), ud, vd);
    tmp3(:,:,2) = interp2(tmp1(:,:,2), ud, vd);
    tmp3(:,:,3) = interp2(tmp1(:,:,3), ud, vd);
    out=tmp3;
end