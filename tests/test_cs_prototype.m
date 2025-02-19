% 21.02.2017 8:33
% прототип ск
clear all
close all
clc
% примитивы
% детали
% координатная сетка


%--local grid in this cs--
t1=tt_tensor([0 1 1 1;0 1 1 1;0 1 1 1;0 0 0 0],1e-12,[2 2 2 2]);
t2=tt_tensor([0 0 0 0;1 1 1 0;1 1 1 0;1 1 1 0],1e-12,[2 2 2 2]);
t3=tt_tensor([0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 1],1e-12,[2 2 2 2]);
b1=classTTBinaryMask(2,t1); %idea: in this toy case we use TTField instead of Boxes
b2=classTTBinaryMask(2,t2);
b3=classTTBinaryMask(2,t3);

%--boxes--
%в коробочках хранятся коробочные маски. при перезаписи грида трется.
general_boxes=[b1 ;b2; b3];

%--parts--
%в деталях хранятся маски уровня деталей. при перезаписи коробочки трется. При
%перезаписи грида трется
Part1=classPart();
Part1.boxes=general_boxes; %the shortest way to pass links to boxes
Part1.composition='(boxmask(2))|boxmask(3)'; 
Part1.eps=3+1j;
Part1.overall_TTBinaryMask.visualizeField;

