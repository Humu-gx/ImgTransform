clear all;
clc;
fn = 'car2.txt';                                                            % 换成自己的文本文件名
P = textread(fn,'%s');                                                      % 以字符串的形式读取  fn 
rs = size(P);                                                               % 获取 字符串大小 以此来获取 矩阵大小
R = zeros(rs(1),length(P{1})/2);                                            % 创建了一个 [rs(1),length(P{1})/2] 的全零矩阵
for r = 1:rs(1)                  
    ss = P{r};
    for k = 1:length(ss)/2
        pp = fix(k/2)+1;
        R(r,k) = hex2dec(ss(pp:pp+1));
    end
end
size(R)                                                                     %显示大小
IMG  = reshape(R,160,120);                                                  %将 R 转化为 160*120 的矩阵
IMGB = imrotate(IMG,-90);                                                   %图像旋转逆时针90度
img = fliplr(IMGB);                                                         %将矩阵IMGB的列绕垂直轴进行左右翻转
[M, N] = size(img);                                                         %获取图片的长宽
imshow(img);                                                                %显示图片
dot=ginput();                                                               %取四个点，依次是左上，右上，左下，右下,


y=[dot(1,1) dot(2,1) dot(3,1) dot(4,1)];                        %四个原顶点
x=[dot(1,2) dot(2,2) dot(3,2) dot(4,2)];

w=round(sqrt((dot(1,1)-dot(2,1))^2+(dot(1,2)-dot(2,2))^2));     %从原四边形获得新矩形宽
h=round(sqrt((dot(1,1)-dot(3,1))^2+(dot(1,2)-dot(3,2))^2));     %从原四边形获得新矩形高

%这里是新的顶点，我取的矩形,也可以做成其他的形状
%大可以原图像是矩形，新图像是从dot中取得的点组成的任意四边形.:)
Y=[dot(1,1) dot(1,1) dot(1,1)+h*2/3 dot(1,1)+h*2/3];
X=[dot(1,2) dot(1,2)+w*2/3  dot(1,2) dot(1,2)+w*2/3];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %变换后的四个顶点，方程右边的值
%联立解方程组，方程的系数
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);
   0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
   0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
   0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
   0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=A\B;                                         %用四点求得的方程的解，也是全局变换系数
                                                %等价于 inv(A)*B  其中 inv(A)矩阵求逆
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
g=fa(7);h=fa(8);

rot=[
     d e f;
     a b c;
     g h 1];        %公式中第一个数是x,Matlab第一个表示y，所以我矩阵1,2行互换了

pix1=rot*[1 1 1]'/(g*1+h*1+1);  %变换后图像左上点
pix2=rot*[1 N 1]'/(g*1+h*N+1);  %变换后图像右上点
pix3=rot*[M 1 1]'/(g*M+h*1+1);  %变换后图像左下点
pix4=rot*[M N 1]'/(g*M+h*N+1);  %变换后图像右下点
% round()四舍五入取整
height=round(max([pix1(1) pix2(1) pix3(1) pix4(1)])-min([pix1(1) pix2(1) pix3(1) pix4(1)]));     %变换后图像的高度
width=round(max([pix1(2) pix2(2) pix3(2) pix4(2)])-min([pix1(2) pix2(2) pix3(2) pix4(2)]));      %变换后图像的宽度

if min([pix1(1) pix2(1) pix3(1) pix4(1)]) >= 0
    delta_y = -round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])))+1;
else
    delta_y = round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])))-1;             %取得y方向的负轴超出的偏移量
end
if min([pix1(2) pix2(2) pix3(2) pix4(2)]) >= 0
    delta_x = -round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));
else
    delta_x = round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));               %取得x方向的负轴超出的偏移量
end
inv_rot=inv(rot);

AIR=zeros(height,width,'uint8');
gx1=zeros(height,width,'uint8');                        % 预先留出空间 可加快编译速度
gx2=zeros(height,width,'uint8');                        

for z = 1-delta_y:height-delta_y                        %从变换图像中反向寻找原图像的点，以免出现空洞，和旋转放大原理一样
    for j = 1-delta_x:width-delta_x-1
        pix=inv_rot*[z j 1]';       %求原图像中坐标，因为[YW XW W]=fa*[y x 1],所以这里求的是[YW XW W],W=gy+hx+1;
        pix=inv([g*pix(1)-1 h*pix(1);g*pix(2) h*pix(2)-1])*[-pix(1) -pix(2)]'; %相当于解[pix(1)*(gy+hx+1) pix(2)*(gy+hx+1)]=[y x],这样一个方程，求y和x，最后pix=[y x];
        
         if pix(1)>=0.5 && pix(2)>=0.5 && pix(1)<=M && pix(2)<=N
            AIR(z+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));
            gx1(z+delta_y,j+delta_x) = round(pix(1));                       % imgchang（i,j）= img(gx1[i,j],gx2[i,j]) 
            gx2(z+delta_y,j+delta_x) = round(pix(2));                       % 通过 gx1 gx2 可查表对应出逆透视图形
            %AIR(z+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));     %最邻近插值,也可以用双线性或双立方插值
            %if pixell(round(pix(1)),round(pix(2))).x~=0
            %pixelll(z+delta_y,j+delta_x)=round(pix(1))+ round(pix(2))*i ; 
            %AIR(z+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));
            %pixel(z+delta_y,j+delta_x)=(pixell(round(pix(1)),round(pix(2))).x)+ (pixell(round(pix(1)),round(pix(2))).y)*i ;
            %else
            %   AIR(z+delta_y,j+delta_x)=128;
                %pixel(z+delta_y,j+delta_x)=0+ 0*i ;
           %floor(round(pix(2)))+floor(round(pix(1)))*i;   
           %pixel(z+delta_y,j+delta_x).x=floor(round(pix(2)))
           %pixel(z+delta_y,j+delta_x).y=floor(round(pix(1)))
           %end
        else
            AIR(z+delta_y,j+delta_x)=128;
        end  
    end
end
imshow(AIR);


