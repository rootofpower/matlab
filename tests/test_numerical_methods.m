%% ТЕСТУВАННЯ ЧИСЕЛЬНИХ МЕТОДІВ ДЛЯ РОЗВ'ЯЗАННЯ РІВНЯНЬ ТА СИСТЕМ
% Цей скрипт демонструє роботу різних чисельних методів на конкретних прикладах
% та порівнює їх ефективність.

close all;
clear;
clc;

% Додаємо шляхи до директорій з методами і утилітами
addpath('../src/nonlinear');
addpath('../src/linear');
addpath('../src/approximation');
addpath('../src/utils');

fprintf('====== ТЕСТУВАННЯ ЧИСЕЛЬНИХ МЕТОДІВ ======\n\n');

%% ЗАВАНТАЖЕННЯ ВХІДНИХ ДАНИХ З ФАЙЛІВ
% Спробуємо завантажити дані з відповідних файлів, якщо вони існують
nonlinear_data_file = '../data/input/nonlinear_equations_data.txt';
linear_data_file = '../data/input/linear_systems_data.txt';
fixed_point_data_file = '../data/input/fixed_point_data.txt';
least_squares_data_file = '../data/input/least_squares_data.txt';

% Перевіряємо наявність файлів і завантажуємо дані
if exist(nonlinear_data_file, 'file')
    nonlinear_data = read_input_data(nonlinear_data_file);
    fprintf('Завантажено дані для нелінійних рівнянь з файлу.\n');
else
    nonlinear_data = [];
    fprintf('Файл з даними для нелінійних рівнянь не знайдено. Використовуємо значення за замовчуванням.\n');
end

if exist(linear_data_file, 'file')
    linear_data = read_input_data(linear_data_file);
    fprintf('Завантажено дані для систем лінійних рівнянь з файлу.\n');
else
    linear_data = [];
    fprintf('Файл з даними для систем лінійних рівнянь не знайдено. Використовуємо значення за замовчуванням.\n');
end

if exist(fixed_point_data_file, 'file')
    fixed_point_data = read_input_data(fixed_point_data_file);
    fprintf('Завантажено дані для методу простої ітерації з файлу.\n');
else
    fixed_point_data = [];
    fprintf('Файл з даними для методу простої ітерації не знайдено. Використовуємо значення за замовчуванням.\n');
end

if exist(least_squares_data_file, 'file')
    least_squares_data = read_input_data(least_squares_data_file);
    fprintf('Завантажено дані для методу найменших квадратів з файлу.\n');
else
    least_squares_data = [];
    fprintf('Файл з даними для методу найменших квадратів не знайдено. Використовуємо значення за замовчуванням.\n');
end

fprintf('\n');

%% 1. ТЕСТУВАННЯ МЕТОДІВ ДЛЯ НЕЛІНІЙНИХ РІВНЯНЬ
fprintf('1. МЕТОДИ ДЛЯ НЕЛІНІЙНИХ РІВНЯНЬ\n');
fprintf('--------------------------------------\n');

% Визначення тестової функції для методів пошуку коренів
f = @(x) x^3 - x - 2;  % f(x) = x^3 - x - 2 (має корінь x ≈ 1.5214)
df = @(x) 3*x^2 - 1;   % Похідна f'(x) = 3x^2 - 1

% Функція для методу простої ітерації на основі вхідних даних
% За замовчуванням: g = @(x) (x + 2)^(1/3);  % x = (x + 2)^(1/3)
if ~isempty(fixed_point_data) && isfield(fixed_point_data, 'a') && isfield(fixed_point_data, 'b') && ...
   isfield(fixed_point_data, 'c') && isfield(fixed_point_data, 'd') && isfield(fixed_point_data, 'e')
    % Використовуємо перший рядок даних для створення функції
    a = fixed_point_data.a(1);
    b = fixed_point_data.b(1);
    c = fixed_point_data.c(1);
    d = fixed_point_data.d(1);
    e = fixed_point_data.e(1);
    
    % Створюємо функцію виду: ax^4 + bx^3 + cx^2 + dx + e = 0
    % Використовуємо ізольований x: x = f(x)
    g = @(x) nthroot(-(a*x^4 + b*x^3 + c*x^2 + d*x + e), 1);
    fprintf('Використовуємо функцію для методу простої ітерації з файлу: ');
    fprintf('a=%.2f, b=%.2f, c=%.2f, d=%.2f, e=%.2f\n\n', a, b, c, d, e);
else
    % За замовчуванням
    g = @(x) (x + 2)^(1/3);
    fprintf('Використовуємо функцію за замовчуванням для методу простої ітерації: x = (x + 2)^(1/3)\n\n');
end

% Параметри для всіх методів
tol = 1e-6;     % Точність за замовчуванням
max_iter = 100; % Максимальна кількість ітерацій за замовчуванням

% Якщо є дані для точності з файлу, використовуємо їх
if ~isempty(fixed_point_data) && isfield(fixed_point_data, 'eps')
    tol = fixed_point_data.eps(1);
    fprintf('Використовуємо точність з файлу: %.10e\n\n', tol);
end

% 1.1 Метод бісекції
fprintf('\n1.1 Метод бісекції\n');
% Значення за замовчуванням
a = 1; b = 2;  % Відрізок, на якому шукаємо корінь

% Якщо є дані з файлу, використовуємо їх
if ~isempty(nonlinear_data) && isfield(nonlinear_data, 'method') && ...
   isfield(nonlinear_data, 'a') && isfield(nonlinear_data, 'b') && ...
   isfield(nonlinear_data, 'tol') && isfield(nonlinear_data, 'max_iter')
    
    % Знаходимо рядок для методу бісекції
    idx = find(strcmp(nonlinear_data.method, 'bisection'));
    if ~isempty(idx)
        a = nonlinear_data.a(idx);
        b = nonlinear_data.b(idx);
        tol = nonlinear_data.tol(idx);
        max_iter = nonlinear_data.max_iter(idx);
        fprintf('Використовуємо параметри з файлу: a=%.2f, b=%.2f, tol=%.10e, max_iter=%d\n', a, b, tol, max_iter);
    end
end

tic;
[root_bisection, iter_bisection, ~] = bisection_method(f, a, b, tol, max_iter);
time_bisection = toc;
fprintf('Корінь: %.10f\n', root_bisection);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_bisection));
fprintf('Кількість ітерацій: %d\n', iter_bisection);
fprintf('Час виконання: %.6f сек\n', time_bisection);

% 1.2 Метод регула-фалсі
fprintf('\n1.2 Метод регула-фалсі\n');
% Значення за замовчуванням
a = 1; b = 2;

% Якщо є дані з файлу, використовуємо їх
if ~isempty(nonlinear_data) && isfield(nonlinear_data, 'method') && ...
   isfield(nonlinear_data, 'a') && isfield(nonlinear_data, 'b') && ...
   isfield(nonlinear_data, 'tol') && isfield(nonlinear_data, 'max_iter')
    
    % Знаходимо рядок для методу регула-фалсі
    idx = find(strcmp(nonlinear_data.method, 'falsi'));
    if ~isempty(idx)
        a = nonlinear_data.a(idx);
        b = nonlinear_data.b(idx);
        tol = nonlinear_data.tol(idx);
        max_iter = nonlinear_data.max_iter(idx);
        fprintf('Використовуємо параметри з файлу: a=%.2f, b=%.2f, tol=%.10e, max_iter=%d\n', a, b, tol, max_iter);
    end
end

tic;
[root_falsi, iter_falsi, ~] = regula_falsi(f, a, b, tol, max_iter);
time_falsi = toc;
fprintf('Корінь: %.10f\n', root_falsi);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_falsi));
fprintf('Кількість ітерацій: %d\n', iter_falsi);
fprintf('Час виконання: %.6f сек\n', time_falsi);

% 1.3 Метод січних
fprintf('\n1.3 Метод січних\n');
% Значення за замовчуванням
x0 = 1; x1 = 1.5;

% Якщо є дані з файлу, використовуємо їх
if ~isempty(nonlinear_data) && isfield(nonlinear_data, 'method') && ...
   isfield(nonlinear_data, 'x0') && isfield(nonlinear_data, 'x1') && ...
   isfield(nonlinear_data, 'tol') && isfield(nonlinear_data, 'max_iter')
    
    % Знаходимо рядок для методу січних
    idx = find(strcmp(nonlinear_data.method, 'secant'));
    if ~isempty(idx)
        x0 = nonlinear_data.x0(idx);
        x1 = nonlinear_data.x1(idx);
        tol = nonlinear_data.tol(idx);
        max_iter = nonlinear_data.max_iter(idx);
        fprintf('Використовуємо параметри з файлу: x0=%.2f, x1=%.2f, tol=%.10e, max_iter=%d\n', x0, x1, tol, max_iter);
    end
end

tic;
[root_secant, iter_secant, ~] = secant_method(f, x0, x1, tol, max_iter);
time_secant = toc;
fprintf('Корінь: %.10f\n', root_secant);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_secant));
fprintf('Кількість ітерацій: %d\n', iter_secant);
fprintf('Час виконання: %.6f сек\n', time_secant);

% 1.4 Метод Ньютона
fprintf('\n1.4 Метод Ньютона\n');
% Значення за замовчуванням
x0 = 1.5;

% Якщо є дані з файлу, використовуємо їх
if ~isempty(nonlinear_data) && isfield(nonlinear_data, 'method') && ...
   isfield(nonlinear_data, 'x0') && isfield(nonlinear_data, 'tol') && ...
   isfield(nonlinear_data, 'max_iter')
    
    % Знаходимо рядок для методу Ньютона
    idx = find(strcmp(nonlinear_data.method, 'newton'));
    if ~isempty(idx)
        x0 = nonlinear_data.x0(idx);
        tol = nonlinear_data.tol(idx);
        max_iter = nonlinear_data.max_iter(idx);
        fprintf('Використовуємо параметри з файлу: x0=%.2f, tol=%.10e, max_iter=%d\n', x0, tol, max_iter);
    end
end

tic;
[root_newton, iter_newton, ~] = newton_method(f, df, x0, tol, max_iter);
time_newton = toc;
fprintf('Корінь: %.10f\n', root_newton);
fprintf('Значення функції f(x) в корені: %.10e\n', f(root_newton));
fprintf('Кількість ітерацій: %d\n', iter_newton);
fprintf('Час виконання: %.6f сек\n', time_newton);

% 1.5 Метод простої ітерації
fprintf('\n1.5 Метод простої ітерації\n');
% Значення за замовчуванням
x0 = 1.5;

% Якщо є дані з файлу (другий рядок), використовуємо його
if ~isempty(fixed_point_data) && isfield(fixed_point_data, 'a') && isfield(fixed_point_data, 'b') && ...
   isfield(fixed_point_data, 'c') && isfield(fixed_point_data, 'd') && isfield(fixed_point_data, 'e') && ...
   isfield(fixed_point_data, 'eps') && length(fixed_point_data.a) > 1
    
    % Використовуємо другий рядок даних для створення функції
    a = fixed_point_data.a(2);
    b = fixed_point_data.b(2);
    c = fixed_point_data.c(2);
    d = fixed_point_data.d(2);
    e = fixed_point_data.e(2);
    tol = fixed_point_data.eps(2);
    
    % Створюємо функцію виду: ax^4 + bx^3 + cx^2 + dx + e = 0
    g2 = @(x) nthroot(-(a*x^4 + b*x^3 + c*x^2 + d*x + e), 1);
    fprintf('Використовуємо другу функцію для методу простої ітерації з файлу: ');
    fprintf('a=%.2f, b=%.2f, c=%.2f, d=%.2f, e=%.2f, tol=%.10e\n', a, b, c, d, e, tol);
    
    tic;
    [root_fixed2, iter_fixed2, ~] = fixed_point_iteration(g2, x0, tol, max_iter);
    time_fixed2 = toc;
    
    fprintf('Корінь для другої функції: %.10f\n', root_fixed2);
    fprintf('Значення функції в корені: %.10e\n', a*root_fixed2^4 + b*root_fixed2^3 + c*root_fixed2^2 + d*root_fixed2 + e);
    fprintf('Кількість ітерацій: %d\n', iter_fixed2);
    fprintf('Час виконання: %.6f сек\n\n', time_fixed2);
end

% Тестуємо першу функцію методу простої ітерації
tic;
[root_fixed, iter_fixed, ~] = fixed_point_iteration(g, x0, tol, max_iter);
time_fixed = toc;
fprintf('Корінь для першої функції: %.10f\n', root_fixed);
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

% Значення за замовчуванням для тестової системи рівнянь
A = [10, -1, 2; -1, 11, -1; 2, -1, 10];
b = [6; 25; -11];
x0 = zeros(3, 1);

% Параметри для методів
tol = 1e-6;
max_iter = 100;

% Якщо є дані з файлу, використовуємо їх
if ~isempty(linear_data)
    if isfield(linear_data, 'A') && isfield(linear_data, 'b')
        A = linear_data.A;
        b = linear_data.b;
        fprintf('Використовуємо матрицю A і вектор b з файлу.\n');
    end
    
    if isfield(linear_data, 'tol') && isfield(linear_data, 'method')
        % Знаходимо рядок для методу Якобі
        idx = find(strcmp(linear_data.method, 'jacobi'));
        if ~isempty(idx)
            tol_jacobi = linear_data.tol(idx);
            max_iter_jacobi = linear_data.max_iter(idx);
        else
            tol_jacobi = tol;
            max_iter_jacobi = max_iter;
        end
        
        % Знаходимо рядок для методу Гаусса-Зейделя
        idx = find(strcmp(linear_data.method, 'gauss_seidel'));
        if ~isempty(idx)
            tol_gauss = linear_data.tol(idx);
            max_iter_gauss = linear_data.max_iter(idx);
        else
            tol_gauss = tol;
            max_iter_gauss = max_iter;
        end
    else
        tol_jacobi = tol;
        max_iter_jacobi = max_iter;
        tol_gauss = tol;
        max_iter_gauss = max_iter;
    end
else
    tol_jacobi = tol;
    max_iter_jacobi = max_iter;
    tol_gauss = tol;
    max_iter_gauss = max_iter;
end

% 2.1 Метод Якобі
fprintf('\n2.1 Метод Якобі\n');
fprintf('Використовуємо параметри: tol=%.10e, max_iter=%d\n', tol_jacobi, max_iter_jacobi);
tic;
[x_jacobi, iter_jacobi, ~] = jacobi_method(A, b, x0, tol_jacobi, max_iter_jacobi);
time_jacobi = toc;
fprintf('Розв''язок:\n');
disp(x_jacobi);
fprintf('Нев''язка ||Ax-b||: %.10e\n', norm(A*x_jacobi - b));
fprintf('Кількість ітерацій: %d\n', iter_jacobi);
fprintf('Час виконання: %.6f сек\n', time_jacobi);

% 2.2 Метод Гаусса-Зейделя
fprintf('\n2.2 Метод Гаусса-Зейделя\n');
fprintf('Використовуємо параметри: tol=%.10e, max_iter=%d\n', tol_gauss, max_iter_gauss);
tic;
[x_gauss_seidel, iter_gauss_seidel, ~] = gauss_seidel(A, b, x0, tol_gauss, max_iter_gauss);
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

% Створення тестових даних за замовчуванням
n = 20;
x_default = linspace(0, 10, n)';
y_true_default = 2*x_default.^2 - 3*x_default + 1;
rng(42); % Для відтворюваності результатів
noise = randn(size(x_default)) * 5;
y_default = y_true_default + noise;

% Якщо є дані з файлу, використовуємо їх
if ~isempty(least_squares_data) && isfield(least_squares_data, 'x') && isfield(least_squares_data, 'y')
    x = least_squares_data.x;
    y = least_squares_data.y;
    fprintf('Використовуємо дані для апроксимації з файлу.\n');
else
    x = x_default;
    y = y_default;
    fprintf('Використовуємо згенеровані тестові дані для апроксимації.\n');
end

% Тестування методу найменших квадратів для різних степенів поліному
if ~isempty(least_squares_data) && isfield(least_squares_data, 'degree')
    degrees = least_squares_data.degree;
    fprintf('Використовуємо степені поліномів з файлу: ');
    fprintf('%d ', degrees);
    fprintf('\n');
else
    degrees = [1, 2, 3, 5];
    fprintf('Використовуємо степені поліномів за замовчуванням: 1, 2, 3, 5\n');
end

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
    plot(x, y, 'bo', 'MarkerSize', 8, 'DisplayName', 'Вхідні дані');
    hold on;
    
    % Побудова апроксимуючої кривої
    % Перетворення коефіцієнтів для функції polyval
    p = zeros(1, degree + 1);
    for j = 1:degree + 1
        p(j) = coeffs(degree + 2 - j);
    end
    
    x_fine = linspace(min(x), max(x), 100);
    y_fine = polyval(p, x_fine);
    plot(x_fine, y_fine, 'r-', 'LineWidth', 2, 'DisplayName', ['Апроксимація, степінь=' num2str(degree)]);
    
    grid on;
    legend('Location', 'best');
    title(['Апроксимація поліномом степеня ' num2str(degree)]);
    xlabel('x');
    ylabel('y');
end

% Видаляємо шляхи до директорій після завершення
rmpath('../src/nonlinear');
rmpath('../src/linear');
rmpath('../src/approximation');
rmpath('../src/utils');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 