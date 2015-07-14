clear all; close all; clc;
T = load('str.dat');
figure;
hold all;
plot(T(:,1), T(:,2), 'k.');
plot(T(:,1), -T(:,2), 'k.');
plot(-T(:,1), T(:,2), 'k.');
plot(-T(:,1), -T(:,2), 'k.');
