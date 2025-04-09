clear;
clc;
format("rational");
% x = [0 1.5 6.8]; 
% y = [1.45 3.14 4.11];

% x = [0 1 3 5];
% y = [-6 -3 0 1];

 x = [-4 -2 0 1 3]; 
 y = [208 6 -4 3 131];

% x = [0 1 5];
% y = [2 3 147];

% x = [0.38 0.40 0.81 1.25];
% y = [1.462 1.491 2.247 3.490];

% x = [1 2 4];
% y = [3 -5 4];

% x = [-2 1 2 4];
% y = [25 -8 -15 -23];
c = -4; % введи тут той ікс, який треба порахувати

n = length(x);
syms X;
L = 0;
for k = 1:n
    Lk = 1;
    for j = 1:n
        if j ~= k
            Lk = Lk * (X - x(j)) / (x(k) - x(j));
        end
    end
    L = L + y(k) * Lk;
end

%% TABLE
fprintf("Matrica pre Langracny polynom:\n\n")
fprintf("     ")
for i = 0:n - 1
    fprintf("    x%d    ", i);
end
fprintf("\n ----");
for i = 1:n
    fprintf("|---------");
end
fprintf("|\n  x  |");
for i = 1:n
    fprintf("    %g    ", x(i));
end
fprintf("\n ----");
for i = 1:n
    fprintf("|---------");
end
fprintf("|\n f(x)|");
for i = 1:n
    fprintf("    %g    ", y(i));
end
fprintf("\n     ");
fprintf("\n\n");
%% FUNKCII
fprintf("Rovnice pre Langracny polynom:\n")
fprintf("L%d(x)  = \n", n);
for i = 1:n
    num_expr = 1;
    denom_expr = 1;

    num_str = "";
    denom_str = "";

    for j = 1:n
        if i ~= j
            num_expr = num_expr * (X - x(j));
            num_str = num_str + sprintf("(X - %g)", x(j));

            denom_expr = denom_expr * (x(i) - x(j));
            denom_str = denom_str + sprintf("(%g - %g)", x(i), x(j));
        end
    end

    num_simplified = expand(num_expr);
    num_simplified_vpa = vpa(num_simplified, 6);
    
    fprintf("\n         %s                   %s\n", num_str, char(num_simplified_vpa));
    fprintf("  %g *    %s       = %g *    %s\n", y(i), repmat('-', 1, strlength(num_str)), y(i), repmat('-', 1, strlength(char(num_simplified_vpa))));
    fprintf("         %s                    %g\n", denom_str, denom_expr);

    fprintf("\n");
end
fprintf("\n\n");
%% RESULT 
format("rational")
fprintf("Pripravene rovnice:\n\n")
L = expand(L);
fprintf("L(x) = ");
L_vpa = vpa(L, 6);
disp(L_vpa);
L_func = matlabFunction(L);
%% Solve function
format("rational")
L_value = subs(L, X, c);
fprintf("L(%g) = ", c);
disp(vpa(L_value), 6);