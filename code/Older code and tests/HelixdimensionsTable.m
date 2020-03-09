
clc
for i = 1 : 10
[h, D] = Helixinfo(i);

disp([num2str(i) ' & ' num2str(round(h*1e6,2)) ' & ' num2str(round(D*1e6,2)) ' \\ \hline'])
end

