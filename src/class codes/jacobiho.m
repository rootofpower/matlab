% % 10*x1 + x2 + x3 = -9
% % x1 + 8*x2 - x3 = -2
% % -x1 -x2 + 8*x3 = 9
% 
% X = [10 1 1;
%      1 8 -1;
%      -1 -1 8];
% 
% Y = [-9;
%      -2;
%       9];
% 
% mat_size = size(X);
% n = mat_size(1);
% is_dominant = true;
% 
% for i = 1:n
%     diagonal_element = abs(X(i, i));
%     row_sum = sum(abs(X(i, :))) - diagonal_element;
% 
%     if row_sum >= diagonal_element
%         is_dominant = false;
%         return;
%     end
% end
% 
% disp(is_dominant);
% 
% %x1 = -(x2 / 10) - (x3 / 10) - (9 / 10);
% %x2 = -(x2 / 8) + (x3 / 8) - (2 / 8);
% %x3 = (x1 / 8) + (x2 / 8) + (9 /  8);
% 
% A = [0 -0.1 -0.1;     % -(x2 / 10) - (x3 / 10)
%      -0.125 0 0.125;  % -(x2 / 8) + (x3 / 8)
%       0.125 0.125 0]; %  (x1 / 8) + (x2 / 8)
% 
% B = [-0.9;   % - (9 / 10)
%      -0.25;  % - (2 / 8)
%      1.125]; % + (9 / 8)
% 
%       % x1   x2    x3
% x0 = [-0.9 -0.25 1.125]; % <- B
% 
% x1 = x0(1,1);
% x2 = x0(1,2);
% x3 = x0(1,3);
% 
% for i=1:3
%     x1 = -(1 / 10) * x2 - (1/10) * x3 + B(1,1);
%     x2 = -(1/8) * x1 + (1/8) * x3 + B(2,1);
%     x3 = (1/8) * x1 + (1/8) * x2 + B(3,1);
%     disp("Iteration: " + i + ": x1 = " + x1 + ", x2 = " + x2 + ", x3 = " + x3);
% end
% 
% 
% 
% 
% 
% 
% 
% 








%КРУТОЙ МЕТОД ОТ ГПТ

format short;

X = [20 5 4 7 ;
     4 13 4 1;
     5 5 22 8
     1 1 2 5];

Y = [44;
     2;
     68
     19];

n = size(X, 1);

A = zeros(n, n);
B = zeros(n, 1);

for i = 1:n
    A(i, :) = -X(i, :) / X(i, i);
    A(i, i) = 0;
    B(i) = Y(i) / X(i, i);
end

disp('Matica A:');
disp(A);

disp('Vektor B:');
disp(B);

format long;

is_dominant = true;
for i = 1:n
    diagonal_element = abs(X(i, i));
    row_sum = sum(abs(X(i, :))) - diagonal_element;
    if row_sum >= diagonal_element
        is_dominant = false;
        break;
    end
end

if is_dominant
    disp('Matica je riadkovo dominantná.');
else
    disp('Matica NIE JE riadkovo dominantná. Metóda Jacobiho sa nemusí zísť!');
end

alpha_norm = norm(A, inf); 
x_old = B; 
max_iter = 100;
tolerance = 1e-6;

disp('Riešenie metódou Jacobiho:');

for iter = 1:max_iter
    x_new = A * x_old + B;

    error_estimate1 = (alpha_norm / (1 - alpha_norm)) * norm(x_new - x_old, inf);
    error_estimate2 = (alpha_norm^iter / (1 - alpha_norm)) * norm(x_new - B, inf);

    disp(['Iterácia ', num2str(iter), ': x = ', num2str(x_new'), ...
          ', Odhad chyby 1: ', num2str(error_estimate1), ...
          ', Odhad chyby 2: ', num2str(error_estimate2)]);

    if norm(x_new - x_old, inf) < tolerance
        disp(['Sústava sa zbieha po ', num2str(iter), ' iteráciách.']);
        break;
    end

    x_old = x_new;
end

disp('Výsledok:');
disp(x_new);