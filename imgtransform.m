clear all;
clc;
fn = 'car2.txt';                                                            % �����Լ����ı��ļ���
P = textread(fn,'%s');                                                      % ���ַ�������ʽ��ȡ  fn 
rs = size(P);                                                               % ��ȡ �ַ�����С �Դ�����ȡ �����С
R = zeros(rs(1),length(P{1})/2);                                            % ������һ�� [rs(1),length(P{1})/2] ��ȫ�����
for r = 1:rs(1)                  
    ss = P{r};
    for k = 1:length(ss)/2
        pp = fix(k/2)+1;
        R(r,k) = hex2dec(ss(pp:pp+1));
    end
end
size(R)                                                                     %��ʾ��С
IMG  = reshape(R,160,120);                                                  %�� R ת��Ϊ 160*120 �ľ���
IMGB = imrotate(IMG,-90);                                                   %ͼ����ת��ʱ��90��
img = fliplr(IMGB);                                                         %������IMGB�����ƴ�ֱ��������ҷ�ת
[M, N] = size(img);                                                         %��ȡͼƬ�ĳ���
imshow(img);                                                                %��ʾͼƬ
dot=ginput();                                                               %ȡ�ĸ��㣬���������ϣ����ϣ����£�����,


y=[dot(1,1) dot(2,1) dot(3,1) dot(4,1)];                        %�ĸ�ԭ����
x=[dot(1,2) dot(2,2) dot(3,2) dot(4,2)];

w=round(sqrt((dot(1,1)-dot(2,1))^2+(dot(1,2)-dot(2,2))^2));     %��ԭ�ı��λ���¾��ο�
h=round(sqrt((dot(1,1)-dot(3,1))^2+(dot(1,2)-dot(3,2))^2));     %��ԭ�ı��λ���¾��θ�

%�������µĶ��㣬��ȡ�ľ���,Ҳ����������������״
%�����ԭͼ���Ǿ��Σ���ͼ���Ǵ�dot��ȡ�õĵ���ɵ������ı���.:)
Y=[dot(1,1) dot(1,1) dot(1,1)+h*2/3 dot(1,1)+h*2/3];
X=[dot(1,2) dot(1,2)+w*2/3  dot(1,2) dot(1,2)+w*2/3];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %�任����ĸ����㣬�����ұߵ�ֵ
%�����ⷽ���飬���̵�ϵ��
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);
   0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
   0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
   0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
   0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=A\B;                                         %���ĵ���õķ��̵Ľ⣬Ҳ��ȫ�ֱ任ϵ��
                                                %�ȼ��� inv(A)*B  ���� inv(A)��������
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
g=fa(7);h=fa(8);

rot=[
     d e f;
     a b c;
     g h 1];        %��ʽ�е�һ������x,Matlab��һ����ʾy�������Ҿ���1,2�л�����

pix1=rot*[1 1 1]'/(g*1+h*1+1);  %�任��ͼ�����ϵ�
pix2=rot*[1 N 1]'/(g*1+h*N+1);  %�任��ͼ�����ϵ�
pix3=rot*[M 1 1]'/(g*M+h*1+1);  %�任��ͼ�����µ�
pix4=rot*[M N 1]'/(g*M+h*N+1);  %�任��ͼ�����µ�
% round()��������ȡ��
height=round(max([pix1(1) pix2(1) pix3(1) pix4(1)])-min([pix1(1) pix2(1) pix3(1) pix4(1)]));     %�任��ͼ��ĸ߶�
width=round(max([pix1(2) pix2(2) pix3(2) pix4(2)])-min([pix1(2) pix2(2) pix3(2) pix4(2)]));      %�任��ͼ��Ŀ��

if min([pix1(1) pix2(1) pix3(1) pix4(1)]) >= 0
    delta_y = -round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])))+1;
else
    delta_y = round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])))-1;             %ȡ��y����ĸ��ᳬ����ƫ����
end
if min([pix1(2) pix2(2) pix3(2) pix4(2)]) >= 0
    delta_x = -round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));
else
    delta_x = round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));               %ȡ��x����ĸ��ᳬ����ƫ����
end
inv_rot=inv(rot);

AIR=zeros(height,width,'uint8');
gx1=zeros(height,width,'uint8');                        % Ԥ�������ռ� �ɼӿ�����ٶ�
gx2=zeros(height,width,'uint8');                        

for z = 1-delta_y:height-delta_y                        %�ӱ任ͼ���з���Ѱ��ԭͼ��ĵ㣬������ֿն�������ת�Ŵ�ԭ��һ��
    for j = 1-delta_x:width-delta_x-1
        pix=inv_rot*[z j 1]';       %��ԭͼ�������꣬��Ϊ[YW XW W]=fa*[y x 1],�������������[YW XW W],W=gy+hx+1;
        pix=inv([g*pix(1)-1 h*pix(1);g*pix(2) h*pix(2)-1])*[-pix(1) -pix(2)]'; %�൱�ڽ�[pix(1)*(gy+hx+1) pix(2)*(gy+hx+1)]=[y x],����һ�����̣���y��x�����pix=[y x];
        
         if pix(1)>=0.5 && pix(2)>=0.5 && pix(1)<=M && pix(2)<=N
            AIR(z+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));
            gx1(z+delta_y,j+delta_x) = round(pix(1));                       % imgchang��i,j��= img(gx1[i,j],gx2[i,j]) 
            gx2(z+delta_y,j+delta_x) = round(pix(2));                       % ͨ�� gx1 gx2 �ɲ���Ӧ����͸��ͼ��
            %AIR(z+delta_y,j+delta_x)=img(round(pix(1)),round(pix(2)));     %���ڽ���ֵ,Ҳ������˫���Ի�˫������ֵ
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


