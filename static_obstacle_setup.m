%% obstacle setup
chang = 0.5; %half of it
kuan = 1/32; %half of
c_x1 = 1/8;
c_x2 = 0.5;
%%obstacle set up
%cound be time dependent
obs = zeros(M1,M2);
Q_landscape = zeros(M1,M2);
sQ=0.002;
for i = 1:M1
    x1 = (i-1)*hx1-domainx/2;
    for j = 1:M2
        x2 = (j-1)*hx2-domainy/2;
         if abs(x1-c_x1)<=kuan && abs(x2- 0*c_x2)<=chang 
             Q_landscape(i,j) =2000;
         end
         if abs(x1-(-1*c_x1))<=kuan && abs(x2- (0*c_x2))<=chang 
             Q_landscape(i,j) = 2000;
         end
         if abs(x1-c_x1-1/4)<=kuan && abs(x2- 0*c_x2)<=chang 
             Q_landscape(i,j) = 2000;
         end
         if abs(x1-(-1*c_x1)+1/4)<=kuan && abs(x2- (0*c_x2))<=chang 
             Q_landscape(i,j) = 2000;
         end
    end
end