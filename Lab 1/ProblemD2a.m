tic
load("ELE532_Lab1_Data.mat")
rows = size(B,1); %number of rows
cols = size(B,2); %number of columns

for i=1:1:rows
    for j =1:1:cols
        if (abs(B(i,j)) < 0.01)
            B(i,j)=0;
        end
    end
end
toc