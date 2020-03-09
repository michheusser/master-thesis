close all
clear all
clc
%%
%mphstart()

%% REMANENT M
model = mphload('Helix1-10_Detailed_(Remanent M).mph')


%%
load('DemagFactors.mat')

for i = 1 : 1%length(fieldnames(DemagFactors))
    
    
    n_helix = DemagFactors.(['S' num2str(i)]).n_helix;
    
    N_M = DemagFactors.(['S' num2str(i)]).N_simulationM;
    [V, D] = eig(N_M);
    
    for j = 1 : 3
        V(:,j) = StdDirection(V(:,j));
    end
    
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
    set(f1,'Position', [0 0 800 1400])
    f1 = mphplot(model, 'pg2')
    view([1,0,0])
    grid on
    axis([-2e-5 2e-5 -2e-5 2e-5 -5e-5 5e-5])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    title({['Helix: H' num2str(n_helix)],''})
    
    print(['Plots/H' num2str(n_helix) '_x.png'],'-dpng')
    
    
    
end



