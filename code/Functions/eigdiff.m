function [eig_diff, eig_measured] = eigdiff(x,ax_ax,ax_rad,rad_ax,rad_rad,chi_a_sim)

rad1_rad2 = x(1);
rad2_rad1 = x(2);
chi_a = [rad_rad rad2_rad1 ax_rad;...
    rad1_rad2 rad_rad ax_rad;...
    rad_ax rad_ax ax_ax];


eig_m = sort(real(eig(chi_a)));
eig_sim = sort(real(eig(chi_a_sim)));

eig_diff = norm((eig(chi_a) - eig(chi_a_sim))./eig(chi_a_sim));
%eig_diff = norm(eig_m(1:2) - eig_sim(1:2));
%eig_diff = norm(eig_m - eig_sim);

eig_measured = eig_m;

end

