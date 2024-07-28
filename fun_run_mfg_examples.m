para1 = [8e-1,8e-1,1e-1];
para2 = [2e-1,2e-1, 5];
para3 = [5e-1,2e-1, 5];
para4 = [2e-1,5e-1, 5];
filename = "para2_data1e4";
% [a,m,u,w1,w2,w3,w4,rho0] = run_mfg_moving_obs(para2,5e2);
[a,m,u,w1,w2,w3,w4,rho0] = run_mfg_moving_obs_mex(para2,1e4);

save(filename);
