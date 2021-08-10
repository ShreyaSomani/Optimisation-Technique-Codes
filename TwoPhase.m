format short
clear all
clc
%% 2 PHASE METHOD %%%
% To maximize
% Z = -7.5x1 + 3x2
% 3x1 - x2 - x3>=3
% x1 - x2 + x3 >=2
% x1,x2,x3>=0

%% Input parameters 

Variables = {'x_1','x_2','x_3','s_1','s_2','a_1','a_2','Sol'};
OriginalVariables = {'x_1','x_2','x_3','s_1','s_2','Sol'};

OrigC = [-7.5 3 0 0 0 -1 -1 0] % Z= -7.5x1 + 3x2 + 0x3 + 0s1 + 0s2 - a1 - a2 +0Sol 1 extra 0 at end for cost

info = [3 -1 -1 -1 0 1 0; 1 -1 1 0 -1 0 1];
b = [3; 2];
A = [info b];

BV = [6 7] % index of a1,a2 in A

cost = [0 0 0 0 0 -1 -1 0]
startBV = find(cost<0) %a1,a2

%% Calculate Zj-Cj
ZjCj = cost(BV)*A - cost;

%% Printing first table 

ZCj = [ZjCj; A];
Table  = array2table(ZCj);
Table.Properties.VariableNames(1:size(ZCj,2)) = Variables 

%% PHASE - 1 BEGINS : FEASIBILTY CHECK
[BFS,A] =simp(A,BV,cost,Variables);

%% PHASE - 2 BEGINS : OPTIMALITY CHECK

A(:,startBV)=[] %remove artificial variables column
OrigC(:,startBV)=[] %remove artificial variables cost
[OptBFS,OptA] =simp(A,BFS,OrigC,OriginalVariables);

%% Printing solution
 
        FINAL_BFS = zeros(1,size(A,2));
        FINAL_BFS(OptBFS) = OptA(:,end);
        FINAL_BFS(end)= sum(FINAL_BFS.*OrigC);
        OptimalBFS  = array2table(FINAL_BFS);
        OptimalBFS.Properties.VariableNames(1:size(FINAL_BFS,2)) = OriginalVariables

        
