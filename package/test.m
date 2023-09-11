load lena

ir = ones(3,3) / 9;

lena_blurred = conv2(lena, ir, 'same') + randn(size(lena));
lena_deconvolued = udeconv(lena_blurred, ir);

figure(1)
subplot(221)
imagesc(lena, [min(lena(:)) max(lena(:))])
colormap gray
axis image
title('True')

subplot(222)
imagesc(lena_blurred, [min(lena(:)) max(lena(:))])
colormap gray
axis image
title('Data')

subplot(223)
imagesc(lena_deconvolued, [min(lena(:)) max(lena(:))])
colormap gray
axis image
title('Deconvolued')


figure(2)
plot(lena(68,:),'LineWidth',2)
hold on
plot(lena_blurred(68,:),'LineWidth',2)
plot(lena_deconvolued(68,:),'LineWidth',2)
legend('True','Data','Deconv')

delta_2 = norm(lena_deconvolued - lena,2)/norm(lena,2);
s1 = 'Error in norm 2';
delta_1 = norm(lena_deconvolued - lena,1)/norm(lena,1);
s2 = 'Error in norm 1';
delta_inf = norm(lena_deconvolued - lena,inf)/norm(lena,inf);
s3 = 'Error in norm infinity';


figure(3)
semilogy(delta_1,'-*', 'DisplayName',s1);
hold on
semilogy(delta_2,'-*', 'DisplayName',s2)
semilogy(delta_inf,'-*', 'DisplayName',s3)
