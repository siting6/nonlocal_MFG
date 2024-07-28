function [F_phi_updates] = solvePoisson_time_space_2d_neumann_simple_for_mex(M1,M2,N,rho_m_vec,ht,fv)
F_phi_updates = (zeros(M1,M2,N));
F_rho_m_vec = complex(zeros(M1,M2,N));
% do fourier transform
for l =1:(N)
     F_rho_m_vec(:,:,l) =  dct2_simple(rho_m_vec(:,:,l));
end

phi_fouir_part = complex(zeros(M1,M2,N));
negative_onesa =(1*1/ht/ht)*ones(N-1,1);
negative_onesc =negative_onesa ;
for i = 1: (M1) %j is the corresponding fourier mode
    for j = 1:M2
        %prepare the RHS, including t=0 added to first term
            f =  squeeze(F_rho_m_vec(i,j,:));
            cc =  (fv(i,j))- 2/ht/ht;
            thomas_b = cc*ones(N,1);  
            thomas_b(1) =  thomas_b(1)+1/ht/ht -1/ht/ht;
            thomas_b(N) =  thomas_b(N)+0/ht/ht;
            thomas_n = N;
            s =ThomasAlgorithm_complex(negative_onesa,thomas_b,negative_onesc,f,thomas_n);
            phi_fouir_part(i,j,:) = s;
    end
end

for l =1:(N)
     F_phi_updates(:,:,l) = real(idct2_simple(phi_fouir_part(:,:,l)));
end

end
