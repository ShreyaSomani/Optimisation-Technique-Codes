format short
clear all
clc

%% SHREYA SOMANI
%101917127
%CSE -5

%% DUAL SIMPLEX METHOD %%%
% To maximize
% Z = -2x1 -x2
% 3x1 + x2 >= 3
% 4x1 + 3x2 >= 6
% x1 + 2x2 >= 3
% x1,x2 >=0

%% Input parameters 

Variables = {'x_1','x_2','s_1','s_2','s_3','Sol'};
cost = [-2 -1 0 0 0 0 ]; % 1 extra 0 at end for cost
info = [-3 -1;-4 -3; -1 -2]; % negative of coeff
b = [-3;-6;-3]; % negative of RHS
s = eye(size(info,1));
A = [info s b];

%% Constraint BV
BV = [];

for j=1:size(s,2)
    for i=1:size(A,2)
        if A(:,i)==s(:,j)
            BV = [BV i];
        end
    end
end

%% Calculate Zj-Cj
ZjCj = cost(BV)*A - cost;

%% Printing first table 

ZCj = [ZjCj; A];
SimplexTable  = array2table(ZCj);
SimplexTable.Properties.VariableNames(1:size(ZCj,2)) = Variables

%% DUAL SIMPLEX ITERATIONS BEGIN

RUN = true;
while RUN 
%% Finding leaving variable   
        sol = A(:,end);
        if any (sol <0)  % Xb <0
            fprintf('BFS is not feasible \n');
            [leave_Val,pvt_row]=min(sol);
            fprintf('The leaving variable is %d \n',pvt_row);
     
                      
%% Finding entering variable

          ZJ = ZjCj(:,1:end-1);
          row = A(pvt_row,1:end-1);
          for i=1:size(row,2)
                if(row(i)<0)
                    ratio(i) = abs(ZJ(i)./row(i));
                else
                    ratio(i) = inf;
                end
          end
          
           [MinVal,pvt_col] = min(ratio);
           fprintf('Min val corresponding to pivot col is %d \n',MinVal);
           fprintf('entering variable is %d \n',pvt_col);

        
        BV(pvt_row) = pvt_col;
        disp('New basic variables =');
        disp(BV);
        
        %% Finding pivot element      
 
        pvt_key = A(pvt_row,pvt_col);
        
        %%%% row transformations
        A(pvt_row,:)= A(pvt_row,:)./pvt_key;
        
        for i=1:size(A,1)
            if i~=pvt_row
                A(i,:) = A(i,:) - A(i,pvt_col).*A(pvt_row,:);
            end
        end
        
        ZjCj= cost(BV)*A - cost;
        
        
        %% printing 
        ZCj = [ZjCj; A];
        SimplexTable  = array2table(ZCj);
        SimplexTable.Properties.VariableNames(1:size(ZCj,2)) = Variables
        
 
        
    else
        RUN = false;
        fprintf('Current BFS is optimal');
        end
end

%% Printing solution
 
        BFS = zeros(1,size(A,2));
        BFS(BV) = A(:,end);
        BFS(end)= sum(BFS.*cost);
        CurrentBFS  = array2table(BFS);
        CurrentBFS.Properties.VariableNames(1:size(BFS,2)) = Variables


            
            
