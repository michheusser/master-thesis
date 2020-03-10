load('Dataset_H2_H9.mat');

figure;
plot([2 3 4 5 6 7 8 9],mean(MisAngle(:,:)),'--*','MarkerSize',10);hold all;
plot([2 3 4 5 6 7 8 9],(MisAngle(1,:)),'.');hold all;
plot([2 3 4 5 6 7 8 9],(MisAngle(2,:)),'.');hold all;
plot([2 3 4 5 6 7 8 9],(MisAngle(3,:)),'.');hold all;
plot([2 3 4 5 6 7 8 9],(MisAngle(4,:)),'.');
grid on;

ylabel('Misalignment Angle - degrees');xlabel('Helix Number');