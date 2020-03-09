close all
clear all
clc
%%
%mphstart()

%% REMANENT M
model = mphload('Helix1-10_(Remanent M).mph')
%%
%model.result.create('pg2',3)
%%
%model.result('pg2').feature.create('arwv2','ArrowVolume')

%%
load('DemagFactors.mat')
%%

writerObj = VideoWriter('Plots\H1-10.avi');
writerObj.FrameRate = 1;
open(writerObj);


for i = 1 : 9
N_M = DemagFactors.(['H' num2str(i)]).N_simulationM
[V, D] = eig(N_M);

%fixing manually the directions
if(i == 5)
   V(1,2) = -V(1,2);
   V(2,2) = -V(2,2);
   V(3,2) = -V(3,2);
end

if(i == 7)
   V(1,1) = -V(1,1);
   V(2,1) = -V(2,1);
   V(3,1) = -V(3,1);
end
%

if(i == 8)
V = -V;
end



%N_B = DemagFactors.(['H' num2str(i)]).N_simulation
model.result('pg2').feature('arwv1').set('expr',{num2str(V(1,1)/D(1,1)), num2str(V(2,1)/D(1,1)), num2str(V(3,1)/D(1,1))})
model.result('pg2').feature('arwv2').set('expr',{num2str(V(1,2)/D(2,2)), num2str(V(2,2)/D(2,2)), num2str(V(3,2)/D(2,2))})
model.result('pg2').feature('arwv3').set('expr',{num2str(V(1,3)/D(3,3)), num2str(V(2,3)/D(3,3)), num2str(V(3,3)/D(3,3))})
model.result('pg2').feature('arwv1').set('scale','3e-6')
model.result('pg2').feature('arwv2').set('scale','3e-6')
model.result('pg2').feature('arwv3').set('scale','3e-6')
model.result('pg2').feature('arwv1').set('arrowlength','logarithmic')
model.result('pg2').feature('arwv2').set('arrowlength','logarithmic')
model.result('pg2').feature('arwv3').set('arrowlength','logarithmic')

model.result('pg2').set('outersolnum',num2str(i))


f1 = figure(1);
set(f1,'Position', [0 0 1400 400])
s1 = subplot(1,4,1)
s1 = mphplot(model, 'pg2')
view([-1,0,0])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(i)],''})



s2 = subplot(1,4,2)
s2 = mphplot(model, 'pg2')
view([0,-1,0])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(i)],''})




s3 = subplot(1,4,3)
s3 = mphplot(model, 'pg2')
view([0,0,-1])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(i)],''})



s4 = subplot(1,4,4)
s4 = mphplot(model, 'pg2')
view([-1,-1,-1])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(i)],''})

print(['Plots/H' num2str(i) '.png'],'-dpng')

F(i) = getframe(f1);
writeVideo(writerObj,F(i));

end

close(writerObj)

%%
close all

movie(F_x,1,1)
