function [E,R1,R2,t]=getE_6(F)
% get ideas from https://blog.csdn.net/kokerf/article/details/72911561
Z=[0,1,0;-1,0,0;0,0,0];
W=[0,-1,0;1,0,0;0,0,1];
A =[ 895.8176,0,606.7259;0,898.4794,447.4418;0,0,1];
E=A'*F*A;
[U,D,V]=svd(E);
% R有两种情况，而t只有一种
R1=U*W*V';
R2=U*W'*V';
t=U*[0,0,1]';

end

