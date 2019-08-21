function [x,y,z] = triangulation_7(p1,p2,R,t)
% get idea from https://blog.csdn.net/yg838457845/article/details/81293286
A =[ 895.8176,0,606.7259;0,898.4794,447.4418;0,0,1];
n1=A^(-1)*[p1,ones(size(p1,1),1)]'; % 由世界坐标系恢复摄像机坐标系
n2=A^(-1)*[p1,ones(size(p2,1),1)]'; % 由世界坐标系恢复摄像机坐标系
M1=[1,0,0,0;0,1,0,0;0,0,1,0];       % 由摄像机坐标系1到摄像机坐标系1的矩阵
H=[R,t;[zeros(1,3),1]];
H_i=inv(H);
M2=H_i(1:3,1:4);                    % 由摄像机坐标系2到摄像机坐标系1的矩阵
for item=1:size(p1,1)
    tmp1=[0,-n1(3,item),n1(2,item); % 计算每个点对应的反对称矩阵
        n1(3,item),0,-n1(1,item);
        -n1(2,item),n1(1,item),0];
    tmp2=[0,-n2(3,item),n2(2,item);
        n2(3,item),0,-n2(1,item);
        -n2(2,item),n2(1,item),0];
    Q=[tmp1*M1;tmp2*M2];            % tmp2*M2为摄像机坐标系2对应到摄像机坐标系1中
    [U,D,V]=svd(Q);                 % SVD分解
    P=V(:,4);
    P=P/P(4);                       % 归一，z=1
    x(item)=P(1);                   % 计算xyz
    y(item)=P(2);
    z(item)=P(3);
end

end

