%% ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ
% Цей скрипт демонструє роботу різних методів для чисельного обчислення інтегралів
% та порівнює їх ефективність.

close all;
clear;
clc;

% Додаємо шляхи до директорій з методами інтегрування та утилітами
addpath('../src/integration');
addpath('../src/utils');

fprintf('====== ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ ======\n\n');

%% ЗАВАНТАЖЕННЯ ВХІДНИХ ДАНИХ З ФАЙЛУ
integration_data_file = '../data/input/integration_data.txt';

% Перевіряємо наявність файлу і завантажуємо дані
if exist(integration_data_file, 'file')
    integration_data = read_input_data(integration_data_file);
    fprintf('Завантажено дані для інтегрування з файлу.\n');
else
    integration_data = [];
    fprintf('Файл з даними для інтегрування не знайдено. Використовуємо значення за замовчуванням.\n');
end

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
% Значення за замовчуванням
a = 0;      % Нижня межа інтегрування
b = 2;      % Верхня межа інтегрування
n = 100;    % Кількість відрізків за замовчуванням

% Якщо є дані з файлу, використовуємо їх
if ~isempty(integration_data) && isfield(integration_data, 'method') && ...
   isfield(integration_data, 'a') && isfield(integration_data, 'b') && ...
   isfield(integration_data, 'n')
    
    % Використовуємо перший рядок даних
    a = integration_data.a(1);
    b = integration_data.b(1);
    n = integration_data.n(1);
    fprintf('Використовуємо параметри з файлу: a=%.2f, b=%.2f, n=%d\n', a, b, n);
end

% Якщо є дані про функцію для інтегрування, використовуємо її
if ~isempty(integration_data) && isfield(integration_data, 'function_id')
    function_id = integration_data.function_id(1);
    if function_id >= 1 && function_id <= 4
        fprintf('Використовуємо функцію %d з файлу: %s\n', function_id, function_names{function_id});
        test_functions = {test_functions{function_id}};
        exact_integrals = {exact_integrals{function_id}};
        function_names = {function_names{function_id}};
    end
end

%% Ініціалізація масивів для зберігання результатів
n_funcs = length(test_functions);

errors_rectangle = zeros(1, n_funcs);
errors_trapz = zeros(1, n_funcs);
errors_simpson = zeros(1, n_funcs);

time_rectangle = zeros(1, n_funcs);
time_trapz = zeros(1, n_funcs);
time_simpson = zeros(1, n_funcs);

%% Обчислення інтегралів та порівняння
for f_idx = 1:n_funcs
    f = test_functions{f_idx};
    F = exact_integrals{f_idx};
    exact_value = F(a, b);
    
    fprintf('\n\nФункція %d: %s\n', f_idx, function_names{f_idx});
    fprintf('Точне значення інтегралу від %g до %g: %.10f\n', a, b, exact_value);
    fprintf('------------------------------------------------------------\n');
    fprintf('Метод    | Значення     | Абс. похибка | Час (сек)\n');
    fprintf('         | інтегралу    |              | \n');
    fprintf('------------------------------------------------------------\n');
    
    % Метод прямокутників
    tic;
    result_rect = rectangle_method(f, a, b, n);
    time_rectangle(f_idx) = toc;
    errors_rectangle(f_idx) = abs(result_rect - exact_value);
    
    fprintf('Прямокут.| %.10f | %.10e | %.6f\n', ...
        result_rect, errors_rectangle(f_idx), time_rectangle(f_idx));

    % Метод трапецій
    tic;
    result_trapz = trapezoidal_method(f, a, b, n);
    time_trapz(f_idx) = toc;
    errors_trapz(f_idx) = abs(result_trapz - exact_value);
    
    fprintf('Трапецій | %.10f | %.10e | %.6f\n', ...
        result_trapz, errors_trapz(f_idx), time_trapz(f_idx));

    % Метод Сімпсона
    tic;
    result_simpson = simpson_method(f, a, b, n);
    time_simpson(f_idx) = toc;
    errors_simpson(f_idx) = abs(result_simpson - exact_value);
    
    fprintf('Сімпсона | %.10f | %.10e | %.6f\n', ...
        result_simpson, errors_simpson(f_idx), time_simpson(f_idx));
    
    fprintf('------------------------------------------------------------\n');
    
    % Візуалізація результатів
    figure;
    x = linspace(a, b, 1000);
    y = f(x);
    plot(x, y, 'b-', 'LineWidth', 2, 'DisplayName', 'Функція');
    hold on;
    
    % Побудова точок для різних методів
    x_points = linspace(a, b, n+1);
    y_points = f(x_points);
    
    % Метод прямокутників
    for i = 1:n
        rectangle('Position', [x_points(i), 0, (b-a)/n, y_points(i)], ...
                 'FaceColor', [1 0 0 0.2], 'EdgeColor', 'r');
    end
    
    % Метод трапецій
    plot(x_points, y_points, 'g-o', 'MarkerSize', 5, 'DisplayName', 'Точки трапецій');
    
    % Метод Сімпсона
    plot(x_points, y_points, 'm-o', 'MarkerSize', 5, 'DisplayName', 'Точки Сімпсона');
    
    grid on;
    legend('Location', 'best');
    title(['Візуалізація методів інтегрування для функції ', function_names{f_idx}]);
    xlabel('x');
    ylabel('y');
end

% Видаляємо шляхи до директорій після завершення
rmpath('../src/integration');
rmpath('../src/utils');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 