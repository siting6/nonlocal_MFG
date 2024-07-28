function residual = calculate_residual_neumann_uw_runcost(runcost,u,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
% matrix_m_w = 1/hx*(w1 - circshift(w1,1,1)) + 1/hx*( circshift(w2,-1,1) - w2) + 1/hx*(w3 - circshift(w3,1,2))+ 1/hx*( circshift(w4,-1,2) - w4) ;

matrix_m_u = zeros(M1,M2,N);
B_star1 = (circshift(u,-1,1) - u)/hx;
B_star2 = (u - circshift(u,1,1) )/hx;
B_star3 = (circshift(u,-1,2) - u)/hx;
B_star4 = (u - circshift(u,1,2) )/hx;

B_star1(M1,:,:) = 0*B_star1(M1,:,:);
B_star2(1,:,:) = 0*B_star2(1,:,:);
B_star3(:,M2,:) = 0*B_star3(:,M2,:);
B_star4(:,1,:) = 0*B_star4(:,1,:);

component1 =m.* B_star1;
matrix_m_u =matrix_m_u + min(component1,0);
component1 = circshift(m,1,1).*B_star2;
matrix_m_u =matrix_m_u - min(component1,0);

component1 =circshift(m,-1,1).* B_star1;
matrix_m_u =matrix_m_u + max(component1,0);
component1 = m.*B_star2;
matrix_m_u =matrix_m_u - max(component1,0);

component1 =m.* B_star3;
matrix_m_u =matrix_m_u + min(component1,0);
component1 = circshift(m,1,2).*B_star4;
matrix_m_u =matrix_m_u - min(component1,0);

component1 =circshift(m,-1,2).* B_star3 ;
matrix_m_u =matrix_m_u + max(component1,0);
component1 = m.*B_star4;
matrix_m_u =matrix_m_u - max(component1,0);

w1_pada = padarray(w1,[1],0,'both');
diff_w1 = 1/hx*( w1_pada- circshift(w1_pada,1,1) );

w2_pada = padarray(w2,[1],0,'both');
diff_w2 = 1/hx*( circshift(w2_pada,-1,1) - w2_pada);

w3_pada = padarray(w3,[0,1],0,'both');
diff_w3 = 1/hx*( w3_pada- circshift(w3_pada,1,2) );

w4_pada = padarray(w4,[0,1],0,'both');
diff_w4 = 1/hx*( circshift(w4_pada,-1,2) - w4_pada);

matrix_m_w =diff_w1(2:M1+1,:,:) + diff_w2(2:M1+1,:,:)+diff_w3(:,2:M2+1,:) +diff_w4(:,2:M2+1,:);
out =  matrix_m_w+ 1/hx*matrix_m_u;

residual = sqrt(sum(sum(sum(out.* out)))*hx^2*ht);

end

