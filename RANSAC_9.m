function [T] = RANSAC_9(x,y,z,u,v,rate,err)
A =[ 895.8176,0,606.7259;0,898.4794,447.4418;0,0,1];
fx = A(1,1);
fy = A(2,2);
cx = A(1,3);
cy = A(2,3);
    
l=length(x);
n=0;
targetn=floor(rate*l);  % rate为符合的点数满足一定的比例时停止迭代，targetn为此时符合的点数
while n<targetn
    t0=randperm(l,3);
    xx=[];              % 计算每组点在3D下的xyz和第三张图中的uv
    yy=[];
    zz=[];
    uu=[];
    vv=[];
    
    for i=1:3           % 随机产生3组点
        xx=[xx,x(t0(i))];
        yy=[yy,y(t0(i))];
        zz=[zz,z(t0(i))];
        vv=[vv,v(t0(i))];
        uu=[uu,u(t0(i))];
    end
    [R,t]=P3P_9(xx,yy,zz,uu,vv);% P3P计算R和t，可能有多重情况
    for j=1:size(t,2)
        T=[R(:,:,j),t(:,j)];
    for i=1:l
        p=[(u(i)-cx)/fx;(v(i)-cy)/fy;1];    % 先归一化，再计算误差
        p=p./sqrt(sum(p.^2));
        merr=norm(p-T*[x(i);y(i);z(i);1]);
        if merr<err
            n=n+1;
        end

    end
    if (n>=targetn) 
            break;
    else 
            if n>1
                n
            end
            n=0;
    end
    end
end
result=T*[x;y;z;ones(1,l)];
for i=1:l
    result(:,i)=result(:,i)/result(3,i);    % z=1
end
%plot(result(1,:)*fx+cx,result(2,:)*fy+cy,'*');
%figure,plot(u,v,'*');
end
