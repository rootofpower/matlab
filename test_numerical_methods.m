%% ТЕСТУВАННЯ ЧИСЕЛЬНИХ МЕТОДІВ ДЛЯ РОЗВ'ЯЗАННЯ РІВНЯНЬ ТА СИСТЕМ
% Цей скрипт демонструє роботу різних чисельних методів на конкретних прикладах
% та порівнює їх ефективність.

close all;
clear;
clc;

fprintf('====== ТЕСТУВАННЯ ЧИСЕЛЬНИХ МЕТОДІВ ======\n\n');

%% 1. ТЕСТУВАННЯ МЕТОДІВ ДЛЯ НЕЛІНІЙНИХ РІВНЯНЬ
fprintf('1. МЕТОДИ ДЛЯ НЕЛІНІЙНИХ РІВНЯНЬ\n');
fprintf('--------------------------------------\n');

% Визначення тестової функції для методів пошуку коренів
f = @(x) x^3 - x - 2;  % f(x) = x^3 - x - 2 (має корінь x ≈ 1.5214)
df = @(x) 3*x^2 - 1;   % Похідна f'(x) = 3x^2 - 1

% Функція для методу простої ітерації
g = @(x) (x + 2)^(1/3);  % x = (x + 2)^(1/3)

% Параметри для всіх методів
tol = 1e-6;     % Точність
max_iter = 100; % Максимальна кількість ітерацій

% 1.1 Метод бісекції
fprintf('\n1.1 Метод бісекції\n');
a = 1; b = 2;  % Відрізок, на якому шукаємо корінь
tic;
[root_bisection, iter_bisection, ~] = bisection_method(f, a, b, tol, max_iter);
time_bisection = toc;
fprintf('Корінь: %.10f\n', root_bisection);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_bisection));
fprintf('Кількість ітерацій: %d\n', iter_bisection);
fprintf('Час виконання: %.6f сек\n', time_bisection);

% 1.2 Метод регула-фалсі
fprintf('\n1.2 Метод регула-фалсі\n');
a = 1; b = 2;
tic;
[root_falsi, iter_falsi, ~] = regula_falsi(f, a, b, tol, max_iter);
time_falsi = toc;
fprintf('Корінь: %.10f\n', root_falsi);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_falsi));
fprintf('Кількість ітерацій: %d\n', iter_falsi);
fprintf('Час виконання: %.6f сек\n', time_falsi);

% 1.3 Метод січних
fprintf('\n1.3 Метод січних\n');
x0 = 1; x1 = 1.5;
tic;
[root_secant, iter_secant, ~] = secant_method(f, x0, x1, tol, max_iter);
time_secant = toc;
fprintf('Корінь: %.10f\n', root_secant);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_secant));
fprintf('Кількість ітерацій: %d\n', iter_secant);
fprintf('Час виконання: %.6f сек\n', time_secant);

% 1.4 Метод Ньютона
fprintf('\n1.4 Метод Ньютона\n');
x0 = 1.5;
tic;
[root_newton, iter_newton, ~] = newton_method(f, df, x0, tol, max_iter);
time_newton = toc;
fprintf('Корінь: %.10f\n', root_newton);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_newton));
fprintf('Кількість ітерацій: %d\n', iter_newton);
fprintf('Час виконання: %.6f сек\n', time_newton);

% 1.5 Метод простої ітерації
fprintf('\n1.5 Метод простої ітерації\n');
x0 = 1.5;
tic;
[root_fixed, iter_fixed, ~] = fixed_point_iteration(g, x0, tol, max_iter);
time_fixed = toc;
fprintf('Корінь: %.10f\n', root_fixed);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_fixed));
fprintf('Кількість ітерацій: %d\n', iter_fixed);
fprintf('Час виконання: %.6f сек\n', time_fixed);

% Порівняння методів для нелінійних рівнянь
fprintf('\nПорівняння методів для нелінійних рівнянь:\n');
fprintf('----------------------------------------------\n');
fprintf('Метод          | Корінь      | Ітерації | Час (сек)\n');
fprintf('----------------------------------------------\n');
fprintf('Бісекція      | %.10f | %8d | %.6f\n', root_bisection, iter_bisection, time_bisection);
fprintf('Регула-фалсі  | %.10f | %8d | %.6f\n', root_falsi, iter_falsi, time_falsi);
fprintf('Січних        | %.10f | %8d | %.6f\n', root_secant, iter_secant, time_secant);
fprintf('Ньютона       | %.10f | %8d | %.6f\n', root_newton, iter_newton, time_newton);
fprintf('Простої ітер. | %.10f | %8d | %.6f\n', root_fixed, iter_fixed, time_fixed);

%% 2. ТЕСТУВАННЯ МЕТОДІВ ДЛЯ СИСТЕМ ЛІНІЙНИХ РІВНЯНЬ
fprintf('\n\n2. МЕТОДИ ДЛЯ СИСТЕМ ЛІНІЙНИХ РІВНЯНЬ\n');
fprintf('--------------------------------------\n');

% Визначення тестової системи рівнянь
A = [10, -1, 2; -1, 11, -1; 2, -1, 10];
b = [6; 25; -11];
x0 = zeros(3, 1);

% Параметри для методів
tol = 1e-6;
max_iter = 100;

% 2.1 Метод Якобі
fprintf('\n2.1 Метод Якобі\n');
tic;
[x_jacobi, iter_jacobi, ~] = jacobi_method(A, b, x0, tol, max_iter);
time_jacobi = toc;
fprintf('Розв''язок:\n');
disp(x_jacobi);
fprintf('Нев''язка ||Ax-b||: %.10e\n', norm(A*x_jacobi - b));
fprintf('Кількість ітерацій: %d\n', iter_jacobi);
fprintf('Час виконання: %.6f сек\n', time_jacobi);

% 2.2 Метод Гаусса-Зейделя
fprintf('\n2.2 Метод Гаусса-Зейделя\n');
tic;
[x_gauss_seidel, iter_gauss_seidel, ~] = gauss_seidel(A, b, x0, tol, max_iter);
time_gauss_seidel = toc;
fprintf('Розв''язок:\n');
disp(x_gauss_seidel);
fprintf('Нев''язка ||Ax-b||: %.10e\n', norm(A*x_gauss_seidel - b));
fprintf('Кількість ітерацій: %d\n', iter_gauss_seidel);
fprintf('Час виконання: %.6f сек\n', time_gauss_seidel);

% 2.3 Вбудований метод MATLAB (для порівняння)
fprintf('\n2.3 Вбудований метод MATLAB\n');
tic;
x_matlab = A \ b;
time_matlab = toc;
fprintf('Розв''язок:\n');
disp(x_matlab);
fprintf('Нев''язка ||Ax-b||: %.10e\n', norm(A*x_matlab - b));
fprintf('Час виконання: %.6f сек\n', time_matlab);

% Порівняння методів для систем лінійних рівнянь
fprintf('\nПорівняння методів для СЛАР:\n');
fprintf('----------------------------------------------\n');
fprintf('Метод          | Нев''язка     | Ітерації | Час (сек)\n');
fprintf('----------------------------------------------\n');
fprintf('Якобі         | %.10e | %8d | %.6f\n', norm(A*x_jacobi - b), iter_jacobi, time_jacobi);
fprintf('Гаусса-Зейделя| %.10e | %8d | %.6f\n', norm(A*x_gauss_seidel - b), iter_gauss_seidel, time_gauss_seidel);
fprintf('MATLAB        | %.10e |        - | %.6f\n', norm(A*x_matlab - b), time_matlab);

%% 3. ТЕСТУВАННЯ МЕТОДУ НАЙМЕНШИХ КВАДРАТІВ
fprintf('\n\n3. МЕТОД НАЙМЕНШИХ КВАДРАТІВ\n');
fprintf('--------------------------------------\n');

% Створення тестових даних
n = 20;
x = linspace(0, 10, n)';
y_true = 2*x.^2 - 3*x + 1;
rng(42); % Для відтворюваності результатів
noise = randn(size(x)) * 5;
y = y_true + noise;

% Тестування методу найменших квадратів для різних степенів поліному
degrees = [1, 2, 3, 5];
fprintf('\nРезультати апроксимації поліномами різних степенів:\n');

for i = 1:length(degrees)
    degree = degrees(i);
    
    fprintf('\n3.%d Поліном степеня %d\n', i, degree);
    tic;
    [coeffs, y_fit] = least_squares(x, y, degree);
    time_ls = toc;
    
    fprintf('Коефіцієнти поліному:\n');
    disp(coeffs');
    
    % Обчислення похибки апроксимації
    error = norm(y - y_fit) / norm(y);
    fprintf('Відносна похибка апроксимації: %.10e\n', error);
    fprintf('Час виконання: %.6f сек\n', time_ls);
    
    % Візуалізація результатів
    figure;
    plot(x, y, 'bo', 'MarkerSize', 8, 'DisplayName', 'Зашумлені дані');
    hold on;
    
    % Побудова оригінальної кривої
    x_fine = linspace(min(x), max(x), 100);
    y_true_fine = 2*x_fine.^2 - 3*x_fine + 1;
    plot(x_fine, y_true_fine, 'k--', 'LineWidth', 2, 'DisplayName', 'Справжня функція');
    
    % Побудова апроксимуючої кривої
    % Перетворення коефіцієнтів для функції polyval
    p = zeros(1, degree + 1);
    for j = 1:degree + 1
        p(j) = coeffs(degree + 2 - j);
    end
    
    y_fine = polyval(p, x_fine);
    plot(x_fine, y_fine, 'r-', 'LineWidth', 2, 'DisplayName', ['Апроксимація, степінь=' num2str(degree)]);
    
    grid on;
    legend('Location', 'best');
    title(['Апроксимація поліномом степеня ' num2str(degree)]);
    xlabel('x');
    ylabel('y');
end

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 