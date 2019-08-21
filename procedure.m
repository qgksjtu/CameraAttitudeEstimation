clc
clear
%% load images
im1=imread('IMG_0963.JPG');
im2=imread('IMG_0962.JPG');
im3=imread('IMG_0961.JPG');
%% step 3
im1_3=distortion_3(im1);
im2_3=distortion_3(im2);
im3_3=distortion_3(im3);
%% step 4
% �õ��Ҷ�ͼ
I1=rgb2gray(im1_3);
I2=rgb2gray(im2_3);
% ���������
points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2); 
[f1, vpts1] = extractFeatures(I1, points1);
[f2, vpts2] = extractFeatures(I2, points2);
% ������ƥ��
[indexPairs,matchmetric] = matchFeatures(f1, f2, 'Prenormalized', true) ;
indexPairs=indexPairs(matchmetric>0.6*max(matchmetric),:);
matched_pts1 = vpts1(indexPairs(:, 1));
matched_pts2 = vpts2(indexPairs(:, 2));
%% step 5
while true
% ����������ж�ȡ����
p1=matched_pts1.Location;
p2=matched_pts2.Location;
F=RANSAC_5(p1,p2,24,2e-3);
%% step 6
[E,R1,R2,t]=getE_6(F);
%% step 7
[x1,y1,z1]=triangulation_7(p1,p2,R1,t); 
[x2,y2,z2]=triangulation_7(p1,p2,R2,t);
if min(z2)<0                            % ����z>0��һ��㣬��ʱ���ܾ�<0�������ٴε���
    if min(z1)>0
        x2=x1;
        y2=y1;
        z2=z1;
        break;
    end
else
    break;
end
end
%% step 8
I3=rgb2gray(im3_3);
points3 = detectSURFFeatures(I3); 
[f3, vpts3] = extractFeatures(I3, points3);
% MaxRation=1��Ŀ���Ǿ����ܱ��������������
[indexPairs_8,matchmetric_8] = matchFeatures(f1, f3, 'Prenormalized',true,'MaxRatio',1) ;
tmp=indexPairs(:,1);
indexPairs_new=[];
l=[];
% Ѱ��ͬʱ��indexPairs_new��indexPairs�е�������
for i=1:size(indexPairs_8,1)
    if ~isempty(find(tmp==indexPairs_8(i,1), 1))
        indexPairs_new=[indexPairs_new;indexPairs_8(i,:)];
        l=[l;find(tmp==indexPairs_8(i,1), 1)];
    end
end
matched_pts81 = vpts1(indexPairs_new(:, 1));
matched_pts82 = vpts3(indexPairs_new(:, 2));
%% step 9
% ������ȷ��T=[R t] �ù��̺�ʱ�ܳ�
T = RANSAC_9(x2(l),y2(l),z2(l),matched_pts82.Location(:,1),matched_pts82.Location(:,2),0.5,1.5e-1);
%% plot
color=rand(length(l),3);    % ���еĵ������ѡ��ɫ
figure,subplot(121),
for i=1:length(l)
place=-T(:,1:3)'*T(:,4);    % 3Dͼ
plot3(x2(l(i)),y2(l(i)),z2(l(i)),'*','Color',color(i,:)),hold on;  
end
hold on,plot3(place(1),place(2),place(3),'.','markersize',20)
for i=1:length(l)           % ԭʼͼ��
subplot(122),im1_3(floor(matched_pts81.Location(i,2))+(-5:5),floor(matched_pts81.Location(i,1))+(-5:5),1)=color(i,1);
im1_3(floor(matched_pts81.Location(i,2))+(-10:10),floor(matched_pts81.Location(i,1))+(-10:10),2)=color(i,2);
im1_3(floor(matched_pts81.Location(i,2))+(-10:10),floor(matched_pts81.Location(i,1))+(-10:10),3)=color(i,3);
end
imshow(im1_3);

