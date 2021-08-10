function[BFS,A] =simp(A,BV,cost,Variables)

ZjCj = cost(BV)*A - cost;

%% Simplex iterations for AUX LPP start

RUN = true;
while RUN 
    ZC = ZjCj(1:end-1);
    if any(ZC < 0) %check for negative values
        fprintf('Current BFS is not optimal');
        fprintf('\n================Next Iteration Result====================\n');
        
        disp('Old BV= ');
        disp(BV);
        
%% Finding entering variable

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
        Table  = array2table(ZCj);
        Table.Properties.VariableNames(1:size(ZCj,2)) = Variables
        
        BFS(BV) = A(:,end);
        
        
        else
        RUN = false;
        fprintf('Current BFS is optimal');
        BFS = BV;
    end
        
end
end
