format short
clear all
clc
%% SIMPLEX METHOD %%%
% To maximize
% Z = -x1 + 3x2 - 2x3 
% 3x1 - x2 + 2x3 <=7
% -2x1 + 4x2 + 0x3 <=12
% -4x1 + 3x2 + 8x3 <=10
% x1,x2,x3>=0

%% Input parameters 

NoOfVariables = 3;
C = [-1 3 -2];
info = [3 -1 2;-2 4 0; -4 3 8];
b = [7; 12; 10];

s = eye(size(info,1)); %storing slag variables

A = [info s b];

cost = zeros(1,size(A,2));
cost(1:NoOfVariables) = C;

%% Constraint BV
BV = NoOfVariables + 1:1:size(A,2)-1;

%% Calculate Zj-Cj
ZjCj = cost(BV)*A - cost;

%% Printing first table 

ZCj = [ZjCj; A];
SimplexTable  = array2table(ZCj);
SimplexTable.Properties.VariableNames(1:size(ZCj,2)) = {'x_1','x_2','x_3','s_1','s_2','s_3','Sol'}


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
        
 %% row transformations
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
        SimplexTable.Properties.VariableNames(1:size(ZCj,2)) = {'x_1','x_2','x_3','s_1','s_2','s_3','Sol'}
        
 
        
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
        CurrentBFS.Properties.VariableNames(1:size(BFS,2)) = {'x_1','x_2','x_3','s_1','s_2','s_3','Sol'}
            
        
        
