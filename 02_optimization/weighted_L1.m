function y = weighted_L1(x)
y=0;
y1=0;
y2=0;
y3=0;

% mip figure
Nmpz=285;
xstack = reshape(x,Nmpz,128,[]);
x_mip2 = std(xstack,0,3);


y2 = norm(x_mip2(:),1);% 对叠加图稀疏性约束，即约束血管结构稀疏

y3 = norm(x(:),1) ; % 3维时空矩阵的L1范数，即约束微泡轨迹稀疏
y3_scale = y3 / size(x,2);


y = y2 +  y3_scale;% 即本文使用的两种正则项




