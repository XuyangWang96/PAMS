function Vec=W_vec(matrix)
% vectorize the upper triangle of a symmetrical matrix
[x,y]=size(matrix);
for row=1:x
    for col=1:y
        if row==col
            matrix(row,col)=0;
        end
    end
end
Vec=squareform(matrix);
end