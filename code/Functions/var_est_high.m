function [var_m, M_s_ax, M_s_rad] = var_est_high(H_app,M_ax,M_rad,N_end)
M_ax_ud = flipud(M_ax);
M_rad_ud = flipud(M_rad);

dof = 2;
var_m_ud = zeros(length(H_app),1);
M_s_ax = max(M_ax_ud);
M_s_rad = max(M_rad_ud);


for i = 1 : length(H_app)
    M_vec_diff = [M_ax_ud(1:i)-M_s_ax; M_rad_ud(1:i)-M_s_rad];
    var_m_ud(i) = (1/length(M_vec_diff))*sum(M_vec_diff.^2);
end

var_m = flipud(var_m_ud);
end

