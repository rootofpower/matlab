clc;
clear;

X = [0 1 2 3];
Y = [0,5 2,1 8,2 25,6];
format short 

A = [X.^2; X.^1; X.^0];
n = length(A);
I = [];
B = log(Y);

for i = 1:n
    I(i,:) = x.^i;
    for j = 1:n
        
        Res = A(i) .* A(j);
    end

end

disp(A);
disp(I);
disp(B);