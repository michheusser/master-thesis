close all
clear all
clc
%%
%mphstart()

%% REMANENT M
model = mphload('Helix1-10_(Remanent M)_Mesh.mph')
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

Helix_vec = [1.5 1.65 1.75 1.9 2 2.2 2.3 2.35 2.45 2.55 2.6 2.7 2.8 3 3.2 3.3 4 5 6 7.2 8 9.1 10]';

for i = 1 : length(Helix_vec)
N_M = DemagFactors.(['S' num2str(i)]).N_simulationM;
[V, D] = eig(N_M);

%fixing manually the directions

n_helix = SimulationM.(['S' num2str(i)]).n_helix

for j = 1 : 3
    V(:,j) = StdDirection(V(:,j));
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
set(f1,'Position', [0 0 1500 800])

s1 = tight_subplot(1,4,[.001 .001],[.1 .01],[.01 .01])

axes(s1(1))
%s1 = subplot(1,4,1)
mphplot(model, 'pg2')
view([1,0,0])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(n_helix)],''},'FontSize',20)



%s2 = subplot(1,4,2)
%s2 = mphplot(model, 'pg2')
axes(s1(2))
mphplot(model, 'pg2')
view([0,-1,0])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(n_helix)],''},'FontSize',20)


N_M = round(N_M*100)/100;
N_M_str = cell(3,3);

for ii = 1 : 3
    
    for jj = 1 : 3
        
        N_M_str{ii,jj} = num2str(N_M(ii,jj));
        if(strcmp(N_M_str{ii,jj},'0'))
            N_M_str{ii,jj} = '0.00';
        end
        
    end
    
end


% s3 = subplot(1,4,3)
% s3 = mphplot(model, 'pg2')
axes(s1(3))
mphplot(model, 'pg2')
view([0,0,-1])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({'N = ',[N_M_str{1,1},' ' ,N_M_str{1,2}, ' ',N_M_str{1,3}]...
    ,[N_M_str{2,1},' ' ,N_M_str{2,2}, ' ',N_M_str{2,3}],...
    [N_M_str{3,1},' ' ,N_M_str{3,2}, ' ',N_M_str{3,3}],'',''},'HorizontalAlignment', 'center'...
    ,'FontSize',20)


%s4 = subplot(1,4,4)
%s4 = mphplot(model, 'pg2')
axes(s1(4))
mphplot(model, 'pg2')
view([-1,-1,-1])
grid on
axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title({['Helix: H' num2str(n_helix)],''},'FontSize',20)

print(['Plots/H' num2str(n_helix) '.png'],'-dpng')

F(i) = getframe(f1);
writeVideo(writerObj,F(i));
close(f1)

end

close(writerObj)

