load('ELE532_Lab1_Data.mat')
A(:) %D.1 (a) List all elements of the matrix array starting from the leftmost column
A([2 4 7]) %D.1 (b) Lists the 2nd, 4th and 7th element of the matrix array
[A >= 0.2] %D.1 (c) Creates a 5x4 logical array that has 0 in place of the 
% elements less than 0.2 and has 1 in place of the elements equal or greater than 0.2
A([A >= 0.2]) %D.1 (d) Lists all elements from the matrix array that are greater or 
% equal to 0.2
A([A >= 0.2]) = 0 %D.1 (e) Replaces all elements in the matrix array with 0 if they are 
% greater or equal to 0.2




