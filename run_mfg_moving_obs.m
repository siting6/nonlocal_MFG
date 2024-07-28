function [a,m,u,w1,w2,w3,w4,rho0] = run_mfg_moving_obs(para,totalitr)
%solve mfg with non-local interaction
%kernel coefficient fixed!
%para: kernel parameters

%splitting method
%non-local mean-field game
%fourier technique

%with G-prox
% clear
%parameters for pdhg algorithm

stepsize1 =5e-1; %for primal variable m,w
stepsize2 =0.5; %for dual variable u
tau3 = 0.2; %for variable a
stepsize_rho0 = 0.2;
stepsize_m0 = 0.2;
max_iteration =totalitr;

%domain
x1domain =2.0;
x2domain =2.0;
domainx =2;
domainy =2;
time =1.0;
%mesh
M1 =64;
M2 = M1;
N = 32;
hx1 = domainx/M1;
hx2 = domainy/M2;
hx = x1domain/(M1);
ht = time/N;

kernel_sigma1=para(1);
kernel_sigma2=para(2);
kernel_mu = para(3);
% basis functions
f0 = zeros(M1,M2);

for i = 1:M1
    x1 = (i-1)*hx - 0.5*x1domain;
    for j = 1:M2
        x2 = (j-1)*hx- 0.5*x1domain;
        f0(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)));
    end
end
f1=f0;
f2 = f0;
for i = 1:M1
    x1 = (i-1)*hx- 0.5*x1domain;
    for j = 1:M2
        x2 = (j-1)*hx- 0.5*x1domain;
        f1(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x1/kernel_sigma1;
        f2(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2/kernel_sigma2;
    end
end
f11=zeros(M1,M2);
f22 = f11;
f12=f11;

for i = 1:M1
    x1 = (i-1)*hx- 0.5*x1domain;
    for j = 1:M2
        x2 = (j-1)*hx- 0.5*x1domain;
        f11(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x1*x1/(sqrt(2)*kernel_sigma1^2);
        f22(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2*x2/(sqrt(2)*kernel_sigma2^2);
        f12(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2*x1/(kernel_sigma2*kernel_sigma1);
    end
end

f3_1=zeros(M1,M2);
f3_2 = f3_1;
f3_12 = f3_1;
f3_21 = f3_1;
for i = 1:M1
    x1 = (i-1)*hx- 0.5*x1domain;
    for j = 1:M2
        x2 = (j-1)*hx- 0.5*x1domain;
        f3_1(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x1*x1*x1/(sqrt(6)*kernel_sigma1^3);
        f3_2(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2*x2*x2/(sqrt(6)*kernel_sigma2^3);
        f3_12(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2*x2*x1/(sqrt(2)*kernel_sigma2^2*kernel_sigma1);
        f3_21(i,j) =sqrt(kernel_mu)*exp(-0.5*((x1)^2/(kernel_sigma1^2) + (x2)^2/(kernel_sigma2^2)))*x2*x1*x1/(sqrt(2)*kernel_sigma2^1*kernel_sigma1^2);
    end
end


f= cat(3,f0,f1,f2,f11,f22,f12,f3_1,f3_2,f3_12,f3_21);

r = 3;
rr = (r+1)*(r+2)/2; %number of 2d basis

v = ones(rr,1); %diagonal enrites of K^{-1}
K = diag(v);

a = zeros(rr,N); %a_i(t) i =1..rr
a_bar=a;

%% running cost matrix
runcost= ones(M1,M2);

%% initial boundary conditions
g = zeros(M1,M2);
for i = 1:M1
    x1 = i*hx -  hx- 0.5*x1domain;
    for j = 1:M2
        x2 = j*hx - hx- 0.5*x1domain;
        g(i,j) =2*exp(-5*((x1-0)^2+0.05*(x2-0.9)^2))*((x2-0.9)^2 -1)+ 1*x1^2;
    end
end

m0 = zeros(M1,M2);
sigma=0.2^2;
for i = 1:M1
    x1 = i*hx -  hx- 0.5*x1domain;
    for j = 1:M2
        x2 = j*hx - hx- 0.5*x1domain;
        m0(i,j) =  1/sigma^2/(2*pi)*exp(-0.5*(((x2+0.9)^2 + (x1-0)^2)/(sigma^2))) + 1/sigma^2/(2*pi)*exp(-0.5*(((x2+0.9)^2 + (x1-0.4)^2)/(sigma^2))) +1/sigma^2/(2*pi)*exp(-0.5*(((x2+0.9)^2 + (x1-0.8)^2)/(sigma^2)))+1/sigma^2/(2*pi)*exp(-0.5*(((x2+0.9)^2 + (x1+0.4)^2)/(sigma^2)))+1/sigma^2/(2*pi)*exp(-0.5*(((x2+0.9)^2 + (x1+0.8)^2)/(sigma^2)));
    end
end
m0 = m0/sum(sum(m0))/hx/hx;
lambda_rho0 = zeros(M1,M2);
rho0 = m0;
m =repmat(1*m0,[1,1,N]);
u = repmat(g,[1,1,N]);

w1 =  zeros(M1,M2,N);
w2 = w1;
w3 = w1;
w4 = w1;

%% obstacle setup
chang = 0.5; %half of it
kuan = 1/32; %half of
c_x1 = 1/8;
c_x2 = 0.5;

%%obstacle set up
%cound be time dependent
obs = zeros(M1,M2,N);
d_x1 = 0.3;
d_x2 = 0.05;
for l = 1:N
    c_x1 = 0.5-1*l*ht;
    c1_x1 = 0.5-1.5*l*ht;
    for i = 1:M1
        x1 = (i-1)*hx1-domainx/2;
        for j = 1:M2
            x2 = (j-1)*hx2-domainy/2;
            
            if abs(x1-c_x1)<=d_x1 && abs(x2-0.4)<=d_x2
                obs(i,j,l) =500;
            end
            if abs(x1-(c1_x1+0.1))<=d_x1 && abs(x2-0.1)<=d_x2
                obs(i,j,l) =500;
            end
            if abs(x1+c1_x1+0.2)<=d_x1 && abs(x2+0.2)<=d_x2
                obs(i,j,l) =500;
            end
            if abs(x1+c1_x1+0.2)<=d_x1 && abs(x2+0.5)<=d_x2
                obs(i,j,l) =500;
            end
        end
    end
end

%% pdhg
%fourier coefficient for Spatial Laplacian
fLapalacian = zeros(M1,M2);
for i = 1:(M1)
    for j = 1:M2
        fLapalacian(i,j) =-1/hx1/hx1*(2*sin(pi*(i-1)/2/(M1-0)))^2 -1/hx2/hx2*(2*sin(pi*(j-1)/2/(M2-0)))^2;
    end
end




% record_diffm = zeros(max_iteration,1);
% record_diffw1 = zeros(max_iteration,1);
% record_diffu = zeros(max_iteration,1);

% record_residual = zeros(max_iteration,1);
% record_residual2 = zeros(max_iteration,1);
% record_residual3 = zeros(max_iteration,1);

tic
for itr = 1:max_iteration
    %     if mod(itr,500) ==0
    %                 residual_cache = calculate_residual_neumann(m0,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
    %                 record_residual(itr) = residual_cache;
    %                 residual_cache2 = calculate_residual_neumann_uw_runcost(runcost,u,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
    %                 record_residual2(itr) = residual_cache2;
    %                residual_cache3 =   calculate_residual_neumann_hjb_runcost_obstacles(runcost,u,g,a,f,obs,ht,hx,M1,M2,N)
    %                record_residual3(itr) = residual_cache3;
    %     end
    
    
    
    % use -bar to save previous iteration values
    u_bar = u;
    a_bar = a;
    lambda_rho0_bar = lambda_rho0;
    
    %update a
    %calculate the integration
    for r_k = 1:rr
        for l = 1:(N)
            a(r_k,l) = (a(r_k,l) + tau3*hx*hx*(sum(sum(m(:,:,l).*f(:,:,r_k)))))/ (1.0 + tau3*v(r_k));
        end
    end
    lambda_rho0 = lambda_rho0 - stepsize_rho0*( m0 - rho0 );
    
    %dual update
    
    w1_pada = padarray(w1,[1],0,'both');
    diff_w1 = 1/hx*( w1_pada- circshift(w1_pada,1,1) );
    
    w2_pada = padarray(w2,[1],0,'both');
    diff_w2 = 1/hx*( circshift(w2_pada,-1,1) - w2_pada);
    
    w3_pada = padarray(w3,[0,1],0,'both');
    diff_w3 = 1/hx*( w3_pada- circshift(w3_pada,1,2) );
    
    w4_pada = padarray(w4,[0,1],0,'both');
    diff_w4 = 1/hx*( circshift(w4_pada,-1,2) - w4_pada);
    
    matrix_m_w =diff_w1(2:M1+1,:,:) + diff_w2(2:M1+1,:,:)+diff_w3(:,2:M2+1,:) +diff_w4(:,2:M2+1,:);
    
    for l = 2:(N)
        matrix_m_w(:,:,l) = matrix_m_w(:,:,l)+(m(:,:,l)- m(:,:,l-1))/ht;
    end
    l=1;
    matrix_m_w(:,:,l) = matrix_m_w(:,:,l)+(m(:,:,l)- m0)/ht;
    
    
    %     record_diffu(itr) = sum(sum(sum((u-u_bar).^2)));
    %using Gprox
    u_update =solvePoisson_time_space_2d_neumann_simple_for_mex(M1,M2,N,matrix_m_w,ht,fLapalacian);
    u = u + stepsize2 * u_update;
    
    %L2
    %     u = u -  stepsize2*matrix_m_w;
    
    %extra interpolation step
    u_bar = 2*u - u_bar;
    a_bar = 2*a - a_bar;
    lambda_rho0_bar = 2 * lambda_rho0 - lambda_rho0_bar;
    %primal update
    m0 = m0 - stepsize_m0 * (1.0/ht * u_bar(:,:,1) -lambda_rho0_bar);
    B_star1 = -(circshift(u_bar,-1,1) - u_bar)/hx;
    B_star2 = - (u_bar - circshift(u_bar,1,1) )/hx;
    B_star3 = - (circshift(u_bar,-1,2) - u_bar)/hx;
    B_star4 = - (u_bar - circshift(u_bar,1,2) )/hx;
    
    B_star1(M1,:,:) = 0*B_star1(M1,:,:);
    B_star2(1,:,:) = 0*B_star2(1,:,:);
    B_star3(:,M2,:) = 0*B_star3(:,M2,:);
    B_star4(:,1,:) = 0*B_star4(:,1,:);
       
    for l = 1:N
        for i = 1:M1
            for j = 1:M2
                
                if l < N
                    a_star_u = -(u_bar(i,j,l+1)-u_bar(i,j,l))/ht;
                else
                    a_star_u = 1/ht*u_bar(i,j,l);
                end
                p = [B_star1(i,j,l);B_star2(i,j,l);B_star3(i,j,l);B_star4(i,j,l)];
                % note F should include final cost at time T
                %now g = 1, G = m, linear in m.
                
                temp_m = m(i,j,l) + stepsize1*(a_star_u - obs(i,j,l));
                temp_p =  [w1(i,j,l);w2(i,j,l);w3(i,j,l);w4(i,j,l)] + stepsize1*p;
                
                
                if l <N
                    %solve for cubic root
                    pp2 = zeros(4,1);
                    pp2(1) = temp_p(1)*(temp_p(1)>=0);
                    pp2(2) = temp_p(2)*(temp_p(2)<=0);
                    pp2(3) = temp_p(3)*(temp_p(3)>=0);
                    pp2(4) = temp_p(4)*(temp_p(4)<=0);
                    norm_projection = sum(pp2.^2);
                    coeff0 = -stepsize1/2* norm_projection;
                    
                    zz = cubic_poly_solve((stepsize1* dot(squeeze(a_bar(:,l)),squeeze(f(i,j,:)))-temp_m-stepsize1),0,coeff0);
                    zz = zz - stepsize1;
                    if zz<=0
                        m(i,j,l) = 0;
                        w1(i,j,l) = 0;
                        w2(i,j,l) = 0;
                        w3(i,j,l) = 0;
                        w4(i,j,l) = 0;
                    else
                        coeff_p = zz/(zz+stepsize1);
                        m(i,j,l) = zz;
                        w1(i,j,l) = coeff_p*pp2(1);
                        w2(i,j,l) = coeff_p*pp2(2);
                        w3(i,j,l) = coeff_p*pp2(3);
                        w4(i,j,l) = coeff_p*pp2(4);
                    end
                else
                    %solve for cubic root
                    pp2 = zeros(4,1);
                    pp2(1) = temp_p(1)*(temp_p(1)>=0);
                    pp2(2) = temp_p(2)*(temp_p(2)<=0);
                    pp2(3) = temp_p(3)*(temp_p(3)>=0);
                    pp2(4) = temp_p(4)*(temp_p(4)<=0);
                    norm_projection = sum(pp2.^2);
                    coeff0 = -stepsize1/2* norm_projection;
                    zz = cubic_poly_solve((stepsize1* dot(squeeze(a_bar(:,l)),squeeze(f(i,j,:))) + stepsize1*(g(i,j)/ht)-temp_m-stepsize1),0,coeff0);
                    zz = zz - stepsize1;
                    if zz<=0
                        m(i,j,l) = 0;
                        w1(i,j,l) = 0;
                        w2(i,j,l) = 0;
                        w3(i,j,l) = 0;
                        w4(i,j,l) = 0;
                    else
                        coeff_p = zz/(zz+stepsize1);
                        m(i,j,l) = zz;
                        w1(i,j,l) = coeff_p*pp2(1);
                        w2(i,j,l) = coeff_p*pp2(2);
                        w3(i,j,l) = coeff_p*pp2(3);
                        w4(i,j,l) = coeff_p*pp2(4);
                    end
                end
                
            end
        end
    end
    
end
toc

residual_cache1 = calculate_residual_neumann(m0,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
residual_cache2 = calculate_residual_neumann_uw_runcost(runcost,u,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
residual_cache3 =   calculate_residual_neumann_hjb_runcost_obstacles(runcost,u,g,a,f,obs,ht,hx,M1,M2,N)

end

