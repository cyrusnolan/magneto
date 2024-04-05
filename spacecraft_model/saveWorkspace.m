% Create a MatFile object
matObj = matfile('48h_omega_neg.mat','Writable',true);
matObj.("sc") = sc;
