%plot_rho_imtermediate_state

load('para4_data1e4.mat')

%compare different \lambda
xxx = -1+hx1:hx1:1;
yyy =-1+ hx2:hx2:1;
[X,Y] = meshgrid(xxx,yyy);
x = [-1 1];
y = [-1 1];
figure;


rho_lambda_a = m;

p = 1;
xxx2 = -1+hx1/p:hx1/p:1;
yyy2 = -1+hx2/p:hx2/p:1;
[X2,Y2] = meshgrid(yyy2,xxx2);
Q_landscape = zeros(p*M1,p*M2);

for i = 1:M1*p
    x1 = (i-1)*hx1/p-domainx/2;
    for j = 1:M2*p
        x2 = (j-1)*hx2/p-domainy/2;
        if abs(x1-c_x1)<=kuan && abs(x2- 0*c_x2)<=chang
            Q_landscape(i,j) = 100;
        end
        if abs(x1-(-1*c_x1))<=kuan && abs(x2- (0*c_x2))<=chang
            Q_landscape(i,j) = 100;
        end
        if abs(x1-c_x1-1/4)<=kuan && abs(x2- 0*c_x2)<=chang
            Q_landscape(i,j) = 100;
        end
        if abs(x1-(-1*c_x1)+1/4)<=kuan && abs(x2- (0*c_x2))<=chang
            Q_landscape(i,j) = 100;
        end
    end
end
u = 0*Q_landscape;
% largest_inensity = [max(max(rho_lambda_a(:,:,ind))),max(max(rho_lambda_b(:,:,ind))),max(max(rho_lambda_c(:,:,ind))),max(max(rho_lambda_d(:,:,ind)))];
ind = [3,16,29]; %N32 0.1,0.5,0.9
largest_inensity = [max(max(rho_lambda_a(:,:,ind(1)))),max(max(rho_lambda_a(:,:,ind(2)))),max(max(rho_lambda_a(:,:,ind(3))))];

lll = max(largest_inensity);
lim = [0,lll];
h(1) =  subplot(1,3,1);
z=rho_lambda_a(:,:,ind(1)) + obs(:,:,ind(1));
outData = interp2(z,2);
imagesc(x,y,outData);

xlabel('x_2');
ylabel('x_1');
caxis(lim);
axis square
axis xy
ax =gca;
ax.FontSize = 12;
ax.XTick = [-1 0 1];
ax.YTick = [-1 0 1];
hold on
h(2) = subplot(1,3,2);
z=rho_lambda_a(:,:,ind(2))+ obs(:,:,ind(2));

outData = interp2(z,2);
imagesc(x,y,outData);

caxis(lim);

xlabel('x_2');
ylabel('x_1');
axis square
axis xy
ax =gca;
ax.FontSize = 12;
ax.XTick = [-1 0 1];
ax.YTick = [-1 0 1];

h(3) = subplot(1,3,3);
z=rho_lambda_a(:,:,ind(3))+ obs(:,:,ind(3));
outData = interp2(z,2);
imagesc(x,y,outData);
caxis(lim);

xlabel('x_2');
ylabel('x_1');
axis square
axis xy
ax =gca;
ax.FontSize = 12;
ax.XTick = [-1 0 1];
ax.YTick = [-1 0 1];
chang1 = 0.225%changkuan
kuan1 = 0.95
get(h(1),'Position');%[left bottom width height]
set(h(1),'Position',[.07 .05 chang1 kuan1]);
set(h(2),'Position',[.37 .05 chang1 kuan1]);
set(h(3),'Position',[.67 .05 chang1 kuan1]);
% set(h(4),'Position',[.74 .05 chang kuan]);
x0=100;
y0=100;
width=1000;
height=350;
set(gcf,'position',[x0,y0,width,height])
colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);
ax.FontSize =12;
ax.XTick = [-1 0 1];
ax.YTick = [-1 0 1];
hold off