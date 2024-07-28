function residual = calculate_residual_neumann_hjb_runcost_obstacles(runcost,u,g,a,f,obs,ht,hx,M1,M2,N)
% matrix_m_w = 1/hx*(w1 - circshift(w1,1,1)) + 1/hx*( circshift(w2,-1,1) - w2) + 1/hx*(w3 - circshift(w3,1,2))+ 1/hx*( circshift(w4,-1,2) - w4) ;
matrix_ut = zeros(M1,M2,N);
for l = 1:(N-1)
    matrix_ut(:,:,l) = matrix_ut(:,:,l)+(u(:,:,l+1)- u(:,:,l))/ht;
end
l=N;
matrix_ut(:,:,l) = matrix_ut(:,:,l)+(g - u(:,:,l))/ht;

% matrix_du = zeros(M1,M2,N);
B_star1 = (circshift(u,-1,1) - u)/hx;
B_star2 = (u - circshift(u,1,1) )/hx;
B_star3 = (circshift(u,-1,2) - u)/hx;
B_star4 = (u - circshift(u,1,2) )/hx;

B_star1(M1,:,:) = 0*B_star1(M1,:,:);
B_star2(1,:,:) = 0*B_star2(1,:,:);
B_star3(:,M2,:) = 0*B_star3(:,M2,:);
B_star4(:,1,:) = 0*B_star4(:,1,:);

component1 =min( B_star1,0);
component2 =max( B_star2,0);
component3 =min( B_star3,0);
component4 =max( B_star4,0);

matrix_du = 0.5*((component1.^2)+(component2.^2)+(component3.^2) + (component4.^2));
for l  = 1:N
    matrix_du(:,:,l) = matrix_du(:,:,l).*runcost;
end
matrix_af =zeros(M1,M2,N);
for l = 1:N
    for i = 1:M1
        for j = 1:M2
            matrix_af(i,j,l) =  dot(squeeze(a(:,l)),squeeze(f(i,j,:)));
        end
    end
end
%
out =  max(0,-matrix_ut+ matrix_du-matrix_af-obs);
%since for the hjb, we only guarantee the inequality

residual = sqrt(sum(sum(sum(out.* out)))*hx^2*ht);

end

