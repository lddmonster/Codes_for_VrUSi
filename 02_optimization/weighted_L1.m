function y = weighted_L1(x)
y=0;
y1=0;
y2=0;
y3=0;

% mip figure
Nmpz=285;
xstack = reshape(x,Nmpz,128,[]);
x_mip2 = std(xstack,0,3);


y2 = norm(x_mip2,1);% �Ե���ͼϡ����Լ������Լ��Ѫ�ܽṹϡ��

y3 = norm(x(:),1) ; % 3άʱ�վ����L1��������Լ��΢�ݹ켣ϡ��
y3_scale = y3 / size(x,2);


y = y2 +  y3_scale;% ������ʹ�õ�����������




