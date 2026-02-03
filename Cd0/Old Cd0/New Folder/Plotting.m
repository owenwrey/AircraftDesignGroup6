clc; clear; close all

M = 0:0.01:2;

figure
hold on

fuselage = D_q150Gfuselage(M);
plot(M,fuselage,DisplayName="150Fuselage")

cluster = D_qClusterRack(M);
plot(M,cluster,DisplayName="ClusterRack")

KBwing = D_q2KBwing(M);
plot(M,KBwing,DisplayName="2KB wing")

G150wing = D_q150Gwing(M);
plot(M,G150wing,DisplayName="")



legend
