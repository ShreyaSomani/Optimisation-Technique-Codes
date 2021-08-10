format short
clear all
clc
%% 2 PHASE METHOD %%%
% To maximize
% Z = -2x1 -x2
% 3x1 + x2 =3
% 4x1 + 3x2 >= 6
% x1 + 2x2 <=3
% x1,x2>=0

%% Input parameters 

Variables = {'x_1','x_2','x_3','a_1','a_2','Sol'};
M = 1000;
cost = [-1 2 3 -M -M 0]; %1 extra 0 at end for cost
info = [-2 1 3 1 0; 2 3 4 0 1];
b = [2; 1];
A = [info b];
s = eye(size(A,1));

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

%% Simplex iterations start

RUN = true;
while RUN 
    ZC = ZjCj(1:end-1);
    if any(ZC < 0) %check for negative values
        fprintf('Current BFS is not optimal');
        fprintf('\n================Next Iteration Result====================\n');
        
        disp('Old BV= ');
        disp(BV);
        
%% Finding entering variable

        ZC = ZjCj(1:end-1);
        [Enter_col,pvt_col]=min(ZC);
        fprintf('The min element in Zj-Cj row is %d correspnding to column %d \n',Enter_col,pvt_col);
        fprintf('The entering variable is %d \nn',pvt_col);
        
%% Finding leaving variable   

        if all(A(:,pvt_col)<=0)   %no leaving variable all αij < 0
            error ('LPP is unbounded');
            
        %%%% finding min ratios
        else
            sol = A(:,end);
            column = A(:,pvt_col);
            for i=1:size(A,1)
                if(column(i)>0)
                    ratio(i) = sol(i)./column(i);
                else
                    ratio(i) = inf; 
                end
            end
            
            %%%% finding the minimum 
            [MinRatio,pvt_row] = min(ratio);
            fprintf('Min ratio corresponding to pivot row is %d \n',MinRatio);
            fprintf('leaving variable is %d \n',BV(pvt_row));
            
        end
        
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
        
        ZjCj=  ZjCj - ZjCj(pvt_col).*A(pvt_row,:);
        
        
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
