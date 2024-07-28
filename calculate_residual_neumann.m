function residual = calculate_residual_neumann(m0,m,w1,w2,w3,w4,ht,hx,M1,M2,N)
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
residual = sqrt(sum(sum(sum(matrix_m_w.* matrix_m_w)))*hx^2*ht);

end

