% 21.02.2017 1:54
% 21.02.2017 5:20 kinda works
% пробую сделать сведение побитовых масок 
clear all
close all
clc

t1=tt_tensor([0 1 1 1;0 1 1 1;0 1 1 1;0 0 0 0],1e-12,[2 2 2 2]);
t2=tt_tensor([0 0 0 0;1 1 1 0;1 1 1 0;1 1 1 0],1e-12,[2 2 2 2]);
t3=tt_tensor([0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 1],1e-12,[2 2 2 2]);

b1=classTTBinaryMask(2,t1) %idea: in this toy case we use TTField instead of Boxes
b2=classTTBinaryMask(2,t2)
b3=classTTBinaryMask(2,t3)



general_boxes=[b1 ;b2; b3]

Part1=classPart();
Part1.boxes=general_boxes; %the shortest way to pass links to boxes
Part1.composition='(boxmask(2))|boxmask(3)'; 
Part1.overall_TTBinaryMask.visualizeField
Part1.eps=3;


