function [F] = RANSAC_5(p1,p2,targetn,err)
% GET IDEAS FROM https://blog.csdn.net/kokerf/article/details/72630863
l=size(p1,1);                   % 点的数量
n=0;
% 符合的点数大于targetn时停止迭代
while n<targetn
    t=randperm(l,8);            % 产生1-l的8个随机整数
    x1=[];
    x2=[];
    y1=[];
    y2=[];
    n=0;
    for item=1:8                % x1 x2 y1 y2为随机提取的点对的坐标
        x1=[x1,p1(t(item),1)];
        y1=[y1,p1(t(item),2)];
        x2=[x2,p2(t(item),1)];
        y2=[y2,p2(t(item),2)];
    end
    F=eightpoints_5(x1,x2,y1,y2);  % 8点法计算F
    x1=[];
    x2=[];
    y1=[];
    y2=[];
    for item=1:l
        tmp(item)=abs([p1(item,1),p1(item,2),1]*F*[p2(item,1);p2(item,2);1]);% 每个点计算误差
        if tmp(item)<err        % 误差小于err时将该点记录下来
            n=n+1;
            x1=[x1,p1(item,1)]; % 此处记录x1 x2 y1 y2是为了将符合要求的点在循环技术后重新计算一次
            y1=[y1,p1(item,2)];
            x2=[x2,p2(item,1)];
            y2=[y2,p2(item,2)];
        end
        
    end
end
F1=eightpoints_5(x1(1:8),x2(1:8),y1(1:8),y2(1:8));
F2=eightpoints_5(x1(9:16),x2(9:16),y1(9:16),y2(9:16));
F3=eightpoints_5(x1(17:24),x2(17:24),y1(17:24),y2(17:24));
F=F1+F2+F3/3;           % 使用前三组点
end



