function [R,t] = P3P_9(X,Y,Z,u,v)

% get idea from https://www.cnblogs.com/mafuqiang/p/8302663.html
A =[ 895.8176,0,606.7259;0,898.4794,447.4418;0,0,1];
syms f
fx = A(1,1);
fy = A(2,2);
cx = A(1,3);
cy = A(2,3);
u=((u)-cx)/fx;% 归一化
v=((v)-cy)/fy;
n=sqrt(u.^2+v.^2+1);% 归一化
x=u./n;
y=v./n;
z=1./n;
xp=[X(1:3);Y(1:3);Z(1:3)];
xp=xp*diag(1./sqrt(sum(xp.^2,1)));  % P3P函数之前先归一，减少计算量
poses=p3p([x(1:3);y(1:3);z(1:3)],xp);
t=[];
R=[];
for k=1:4
    k=poses(:,1+(k-1)*4:4+(k-1)*4);
    if isreal(k)   % 保留实数的k，然后分解出对应的R与t
        t=[t,k(:,1)];
        R(:,:,size(t,2))=k(:,2:4);
    end
end
end

function poses = p3p( worldp, imagep )
poses = zeros(3,16);
P1 = worldp(:,1);
P2 = worldp(:,2);
P3 = worldp(:,3);
% 三点共线情况下提前停止迭代
if norm(cross(P2 - P1,P3 - P1)) == 0
    return
end
f1 = imagep(:,1);
f2 = imagep(:,2);
f3 = imagep(:,3);
e1 = f1;
% 归一化
e3 = cross(f1,f2);
norm_e3 = sqrt(sum(e3.*e3));
e3 = e3 ./ repmat(norm_e3,3,1);
e2 = cross(e3,e1);
T = [ e1'; e2'; e3' ];
f3 = T*f3;
if ( f3(3,1) > 0 )
    f1 = imagep(:,2);
    f2 = imagep(:,1);
    f3 = imagep(:,3);
    e1 = f1;
    e3 = cross(f1,f2);
    norm_e3 = sqrt(sum(e3.*e3));
    e3 = e3 ./ repmat(norm_e3,3,1);
    e2 = cross(e3,e1);
    T = [ e1'; e2'; e3' ];
    f3 = T*f3;
    P1 = worldp(:,2);
    P2 = worldp(:,1);
    P3 = worldp(:,3);
end
% 计算对应的系数
n1 = P2-P1;
norm_n1 = sqrt(sum(n1.*n1));
n1 = n1 ./ repmat(norm_n1,3,1);
n3 = cross(n1,(P3-P1));
norm_n3 = sqrt(sum(n3.*n3));
n3 = n3 ./ repmat(norm_n3,3,1);
n2 = cross(n3,n1);
N = [ n1'; n2'; n3' ];
P3 = N*(P3-P1);
ab = sqrt(sum((P2 - P1).^2));% 距离
v = f3(1,1)/f3(3,1);
w = f3(2,1)/f3(3,1);
c1 = P3(1,1);   % 三个cos
c2 = P3(2,1);
c3 = f1'*f2;
b = 1/(1-c3^2) - 1;
if c3 < 0
    b = -sqrt(b);
else
    b = sqrt(b);
end
% 参考opencv代码，减少计算量
v2 = v^2;
w2 = w^2;
c12 = c1^2;
c13 = c12 * c1;
c14 = c13 * c1;
c22 = c2^2;
c23 = c22 * c2;
c24 = c23 * c2;
d2 = ab^2;
b2 = b^2;
%% 计算对应的系数
factor_4 = -w2*c24 ...
    -c24*v2 ...
    -c24;
factor_3 = 2*c23*ab*b ...
    +2*w2*c23*ab*b ...
    -2*w*c23*v*ab;
factor_2 = -w2*c22*c12 ...
    -w2*c22*d2*b2 ...
    -w2*c22*d2 ...
    +w2*c24 ...
    +c24*v2 ...
    +2*c1*c22*ab ...
    +2*v*w*c1*c22*ab*b ...
    -c22*c12*v2 ...
    +2*c1*c22*w2*ab ...
    -c22*d2*b2 ...
    -2*c12*c22;
factor_1 = 2*c12*c2*ab*b ...
    +2*w*c23*v*ab ...
    -2*w2*c23*ab*b ...
    -2*c1*c2*d2*b;
factor_0 = -2*w*c22*v*c1*ab*b ...
    +w2*c22*d2 ...
    +2*c13*ab ...
    -c12*d2 ...
    +w2*c22*c12 ...
    -c14 ...
    -2*w2*c22*c1*ab ...
    +c22*v2*c12 ...
    +w2*c22*d2*b2;
x=solveQuartic([factor_4 factor_3 factor_2 factor_1 factor_0]);
%% 计算返回的poses矩阵（poses=[t R])
for i=1:4
    cot_a = (-v*c1/w-real(x(i))*c2+ab*b)/(-v*real(x(i))*c2/w+c1-ab);
    cos_t = real(x(i));
    sin_t = sqrt(1-real(x(i))^2);
    sin_a = sqrt(1/(cot_a^2+1));
    cos_a = sqrt(1-sin_a^2);
    if cot_a < 0
        cos_a = -cos_a;
    end
    C = [ ab*cos_a*(sin_a*b+cos_a);
        cos_t*ab*sin_a*(sin_a*b+cos_a);
        sin_t*ab*sin_a*(sin_a*b+cos_a) ];
    C = P1 + N'*C;
    R = [ -cos_a -sin_a*cos_t -sin_a*sin_t;
        sin_a -cos_a*cos_t -cos_a*sin_t;
        0 -sin_t cos_t ];
    R = N'*R'*T;
    poses(1:3,(i-1)*4+1) = C;
    poses(1:3,(i-1)*4+2:(i-1)*4+4) = R;
end
end
% 求解4次方程对应的根，该子函数来自网络
function roots = solveQuartic( factors )

A = factors(1,1);
B = factors(1,2);
C = factors(1,3);
D = factors(1,4);
E = factors(1,5);

A_pw2 = A*A;
B_pw2 = B*B;
A_pw3 = A_pw2*A;
B_pw3 = B_pw2*B;
A_pw4 = A_pw3*A;
B_pw4 = B_pw3*B;

alpha = -3*B_pw2/(8*A_pw2)+C/A;
beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

alpha_pw2 = alpha*alpha;
alpha_pw3 = alpha_pw2*alpha;

P = -alpha_pw2/12-gamma;
Q = -alpha_pw3/108+alpha*gamma/3-beta^2/8;
R = -Q/2+sqrt(Q^2/4+P^3/27);
U = R^(1/3);

if U == 0
    y = -5*alpha/6-Q^(1/3);
else
    y = -5*alpha/6-P/(3*U)+U;
end

w = sqrt(alpha+2*y);

roots(1,1) = -B/(4*A) + 0.5*(w+sqrt(-(3*alpha+2*y+2*beta/w)));
roots(2,1) = -B/(4*A) + 0.5*(w-sqrt(-(3*alpha+2*y+2*beta/w)));
roots(3,1) = -B/(4*A) + 0.5*(-w+sqrt(-(3*alpha+2*y-2*beta/w)));
roots(4,1) = -B/(4*A) + 0.5*(-w-sqrt(-(3*alpha+2*y-2*beta/w)));

end
