function F=eightpoints(x10,x20,y10,y20)
% get idea from https://blog.csdn.net/kokerf/article/details/72630863?locationNum=2&fps=1
A=[];
%% 归一化
x1=x10-min(x10);
x1=x1/max(x1);
x2=x20-min(x20);
x2=x2/max(x2);
y1=y10-min(y10);
y1=y1/max(y1);
y2=y20-min(y20);
y2=y2/max(y2);
%% T为记录的缩小倍数，计算F时使用
T1=[x1;y1;ones(1,8)]/[x10;y10;ones(1,8)];
T2=[x2;y2;ones(1,8)]/[x20;y20;ones(1,8)];
%% 计算A矩阵
for  i=1:8      
    A=[A;[x1(i)*x2(i),x2(i)*y1(i),x2(i),x1(i)*y2(i),y1(i)*y2(i),y2(i),x1(i),y1(i),1]];
end
%% SVD分解，然后重建F
[u,d,v]=svd(A);
F=reshape(v(:,end),3,3)';
[U,D,V]=svd(F);
D(3,3)=0;
F=T2'*U*D*V'*T1;
end




