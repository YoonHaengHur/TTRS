load('gl_continuous5.mat')

err5_a = err_a;
err5_e_avg = mean(err_e);
err5_e_std = std(err_e);
err5_t = sqrt(err5_a.^2 + err5_e_avg.^2);


load('gl_continuous10.mat')

err10_a = err_a;
err10_e_avg = mean(err_e);
err10_e_std = std(err_e);
err10_t = sqrt(err10_a.^2 + err10_e_avg.^2);

load('gl_continuous15.mat')

err15_a = err_a;
err15_e_avg = mean(err_e);
err15_e_std = std(err_e);
err15_t = sqrt(err15_a.^2 + err15_e_avg.^2);
