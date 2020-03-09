close all
clear all
clc

format long

Data_old = dlmread('CobaltSteelVacoflux50_BH.txt');
Data_new = Data_old;

m_s_old = Data_old(end,2);
m_s_new = m_s_old.*0.7;

Data_new(:,2) = Data_new(:,2)*m_s_new/m_s_old;

plot(Data_old(:,1),Data_old(:,2), Data_new(:,1), Data_new(:,2))

dlmwrite('CobaltSteelVacoflux50_BH_ScaledTest.txt', Data_new, ' ')


Data_new_switched = [Data_new(:,2) Data_new(:,1)];
%plot(Data_new_switched(:,1), Data_new_switched(:,2))
dlmwrite('CobaltSteelVacoflux50_HB_ScaledTest.txt', Data_new_switched, ' ')