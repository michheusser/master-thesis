function MatrixLatex( Matrix )

[row, col] = size(Matrix);

disp('\left[\begin{matrix}')

for ii = 1 : row
    
    curr_line = [];
    for i = 1 : col
        
        if(i == col)
            
            curr_line = [curr_line num2str(Matrix(ii,i)) ' \\'];
            
        else
            
            curr_line = [curr_line num2str(Matrix(ii,i)) ' & '];
            
        end
    end

    disp(curr_line)
    
end

disp('\end{matrix}\right]')

end

