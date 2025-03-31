%% ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ
% Цей скрипт демонструє роботу різних методів для чисельного обчислення інтегралів
% та порівнює їх ефективність.

close all;
clear;
clc;

% Додаємо шлях до директорії з методами інтегрування
addpath('../src/integration');

fprintf('====== ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ ======\n\n');

%% Тестові функції для інтегрування
% Визначимо ряд функцій з відомими аналітичними інтегралами для порівняння точності

% 1. Поліном: f(x) = x^2
f1 = @(x) x.^2;
F1 = @(a, b) (b^3 - a^3) / 3;  % аналітичний інтеграл
name1 = 'x^2';

% 2. Тригонометрична: f(x) = sin(x)
f2 = @(x) sin(x);
F2 = @(a, b) cos(a) - cos(b);  % аналітичний інтеграл
name2 = 'sin(x)';

% 3. Експоненційна: f(x) = e^x
f3 = @(x) exp(x);
F3 = @(a, b) exp(b) - exp(a);  % аналітичний інтеграл
name3 = 'e^x';

% 4. Складніша функція: f(x) = x*sin(x^2)
f4 = @(x) x .* sin(x.^2);
F4 = @(a, b) (1 - cos(b^2))/2 - (1 - cos(a^2))/2;  % аналітичний інтеграл
name4 = 'x*sin(x^2)';

% Масиви для зберігання тестових функцій та їх аналітичних інтегралів
test_functions = {f1, f2, f3, f4};
exact_integrals = {F1, F2, F3, F4};
function_names = {name1, name2, name3, name4};

%% Параметри для тестування
a = 0;      % Нижня межа інтегрування
b = 2;      % Верхня межа інтегрування
n_values = [10, 20, 50, 100, 200, 500, 1000]; % Різні кількості відрізків для тестування

%% Ініціалізація масивів для зберігання результатів
n_funcs = length(test_functions);
n_n = length(n_values);

errors_rectangle = zeros(n_funcs, n_n);
errors_trapz = zeros(n_funcs, n_n);
errors_simpson = zeros(n_funcs, n_n);

time_rectangle = zeros(n_funcs, n_n);
time_trapz = zeros(n_funcs, n_n);
time_simpson = zeros(n_funcs, n_n);

%% Обчислення інтегралів та порівняння
for f_idx = 1:n_funcs
    f = test_functions{f_idx};
    F = exact_integrals{f_idx};
    exact_value = F(a, b);
    
    fprintf('\n\nФункція %d: %s\n', f_idx, function_names{f_idx});
    fprintf('Точне значення інтегралу від %g до %g: %.10f\n', a, b, exact_value);
    fprintf('------------------------------------------------------------\n');
    fprintf('Метод    | Кількість | Значення     | Абс. похибка | Час (сек)\n');
    fprintf('         | відрізків | інтегралу    |              | \n');
    fprintf('------------------------------------------------------------\n');
    
    for n_idx = 1:n_n
        n = n_values(n_idx);
        
        % Метод прямокутників
        tic;
        result_rect = rectangle_method(f, a, b, n);
        time_rectangle(f_idx, n_idx) = toc;
        errors_rectangle(f_idx, n_idx) = abs(result_rect - exact_value);
        
        fprintf('Прямокут.| %9d | %.10f | %.10e | %.6f\n', ...
            n, result_rect, errors_rectangle(f_idx, n_idx), time_rectangle(f_idx, n_idx));
        
        % Метод трапецій
        tic;
        result_trapz = trapezoidal_method(f, a, b, n);
        time_trapz(f_idx, n_idx) = toc;
        errors_trapz(f_idx, n_idx) = abs(result_trapz - exact_value);
        
        fprintf('Трапецій | %9d | %.10f | %.10e | %.6f\n', ...
            n, result_trapz, errors_trapz(f_idx, n_idx), time_trapz(f_idx, n_idx));
        
        % Метод Сімпсона
        tic;
        result_simpson = simpson_method(f, a, b, n);
        time_simpson(f_idx, n_idx) = toc;
        errors_simpson(f_idx, n_idx) = abs(result_simpson - exact_value);
        
        fprintf('Сімпсона | %9d | %.10f | %.10e | %.6f\n', ...
            n, result_simpson, errors_simpson(f_idx, n_idx), time_simpson(f_idx, n_idx));
        
        fprintf('------------------------------------------------------------\n');
    end
    
    % Побудова графіків залежності помилки від кількості відрізків
    figure;
    loglog(n_values, errors_rectangle(f_idx, :), 'r.-', 'LineWidth', 2, 'MarkerSize', 15);
    hold on;
    loglog(n_values, errors_trapz(f_idx, :), 'g.-', 'LineWidth', 2, 'MarkerSize', 15);
    loglog(n_values, errors_simpson(f_idx, :), 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
    grid on;
    xlabel('Кількість відрізків (n)');
    ylabel('Абсолютна похибка');
    title(['Залежність похибки від кількості відрізків для функції ', function_names{f_idx}]);
    legend('Метод прямокутників', 'Метод трапецій', 'Метод Сімпсона');
    set(gca, 'XTick', n_values);
end

%% Порівняння швидкості збіжності різних методів
% Для кожної функції побудуємо графік порядку збіжності (log-log графік)
figure;
f_idx = 1;  % Можна змінити на іншу функцію
loglog(n_values, errors_rectangle(f_idx, :), 'r.-', 'LineWidth', 2, 'MarkerSize', 15);
hold on;
loglog(n_values, errors_trapz(f_idx, :), 'g.-', 'LineWidth', 2, 'MarkerSize', 15);
loglog(n_values, errors_simpson(f_idx, :), 'b.-', 'LineWidth', 2, 'MarkerSize', 15);

% Для порівняння: теоретичний порядок збіжності
c = 0.1;  % Константа для зсуву графіка
loglog(n_values, c ./ n_values.^2, 'r--', 'LineWidth', 1);
loglog(n_values, c ./ n_values.^2, 'g--', 'LineWidth', 1);
loglog(n_values, c ./ n_values.^4, 'b--', 'LineWidth', 1);

grid on;
xlabel('Кількість відрізків (n)');
ylabel('Абсолютна похибка');
title(['Порядок збіжності різних методів для функції ', function_names{f_idx}]);
legend('Метод прямокутників', 'Метод трапецій', 'Метод Сімпсона', ...
       'O(h^2)', 'O(h^2)', 'O(h^4)', 'Location', 'southwest');
set(gca, 'XTick', n_values);

% Видаляємо шлях до директорії з методами інтегрування після завершення
rmpath('../src/integration');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 