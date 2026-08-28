function y = soft2(x,T)

y = max(abs(x) - T, 0);
y = y./(y+T) .* x;


% xstack = reshape(x,325,128,1000);
% x_mip2 = std(xstack(:,:,100),0,3);

% figure(2)
% subplot(121)
% I = reshape(abs(x_mip2).^.5,285,128);
% imagesc(I);
% colormap(gray);


