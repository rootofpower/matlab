%% ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА
% Цей скрипт демонструє роботу методу інтерполяції Лагранжа
% та порівнює його з вихідною функцією

close all;
clear;
clc;

% Додаємо шлях до директорії з методами інтерполяції
addpath('../src/interpolation');

fprintf('====== ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА ======\n\n');

%% Тестові функції
% 1. Поліном: f(x) = x^2
f1 = @(x) x.^2;
name1 = 'x^2';

% 2. Тригонометрична: f(x) = sin(x)
f2 = @(x) sin(x);
name2 = 'sin(x)';

% 3. Експоненційна: f(x) = e^x
f3 = @(x) exp(x);
name3 = 'e^x';

% 4. Складніша функція: f(x) = x*sin(x^2)
f4 = @(x) x .* sin(x.^2);
name4 = 'x*sin(x^2)';

% Масиви для зберігання тестових функцій
test_functions = {f1, f2, f3, f4};
function_names = {name1, name2, name3, name4};

%% Параметри для тестування
a = 0;      % Нижня межа інтервалу
b = 2;      % Верхня межа інтервалу
n_points = 5;  % Кількість вузлів інтерполяції
n_interp = 100;  % Кількість точок для побудови графіка

%% Тестування для кожної функції
for f_idx = 1:length(test_functions)
    f = test_functions{f_idx};
    
    % Генеруємо вузли інтерполяції
    x = linspace(a, b, n_points);
    y = f(x);
    
    % Генеруємо точки для побудови графіка
    x_interp = linspace(a, b, n_interp);
    
    % Обчислюємо значення інтерполяційного поліному
    tic;
    [y_interp, L] = lagrange_interpolation(x, y, x_interp);
    time = toc;
    
    % Обчислюємо значення вихідної функції
    y_exact = f(x_interp);
    
    % Обчислюємо похибку
    error = max(abs(y_interp - y_exact));
    
    % Виводимо результати
    fprintf('\nФункція %d: %s\n', f_idx, function_names{f_idx});
    fprintf('Кількість вузлів: %d\n', n_points);
    fprintf('Максимальна похибка: %.10e\n', error);
    fprintf('Час обчислення: %.6f сек\n', time);
    
    % Виводимо коефіцієнти поліному
    fprintf('\nКоефіцієнти поліному Лагранжа (від старшого до молодшого степеня):\n');
    format rational;
    disp(L);
    format;
    
    % Візуалізація результатів
    figure;
    plot(x_interp, y_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Вихідна функція');
    hold on;
    plot(x_interp, y_interp, 'r--', 'LineWidth', 2, 'DisplayName', 'Поліном Лагранжа');
    plot(x, y, 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли інтерполяції');
    
    grid on;
    legend('Location', 'best');
    title(['Інтерполяція Лагранжа для функції ', function_names{f_idx}]);
    xlabel('x');
    ylabel('y');
    
    % Візуалізація базисних поліномів
    figure;
    for i = 1:n_points
        plot(x_interp, L(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('L_%d', i));
        hold on;
    end
    plot(x, ones(1,n_points), 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли');
    
    grid on;
    legend('Location', 'best');
    title(['Базисні поліноми Лагранжа для функції ', function_names{f_idx}]);
    xlabel('x');
    ylabel('y');
end

% Видаляємо шлях до директорії після завершення
rmpath('../src/interpolation');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 