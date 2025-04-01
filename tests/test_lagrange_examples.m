%% ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА НА КОНКРЕТНИХ ПРИКЛАДАХ
% Цей скрипт перевіряє роботу методу Лагранжа на конкретних прикладах

close all;
clear;
clc;

% Додаємо шлях до директорії з методами інтерполяції
addpath('../src/interpolation');

fprintf('====== ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА ======\n\n');

%% Приклад 1: Ваш приклад
fprintf('Приклад 1:\n');
x1 = [0, 1, 2, 3];
y1 = [0.5, 2.1, 8.2, 25.6];
fprintf('Вузли інтерполяції:\n');
fprintf('x = [%s]\n', num2str(x1));
fprintf('y = [%s]\n', num2str(y1));

% Обчислюємо поліном Лагранжа
[y_interp1, L1] = lagrange_interpolation(x1, y1, x1);

% Виводимо коефіцієнти
fprintf('\nКоефіцієнти поліному Лагранжа (від старшого до молодшого степеня):\n');
format rational;
disp(L1);
format;

% Візуалізація
figure;
plot(x1, y1, 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли інтерполяції');
hold on;
x_interp = linspace(min(x1), max(x1), 100);
y_interp = polyval(L1, x_interp);
plot(x_interp, y_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Поліном Лагранжа');
grid on;
legend('Location', 'best');
title('Інтерполяція Лагранжа для Прикладу 1');
xlabel('x');
ylabel('y');

%% Приклад 2: Квадратична функція
fprintf('\nПриклад 2:\n');
x2 = [-1, 0, 1, 2];
y2 = [1, 0, 1, 4];  % f(x) = x^2
fprintf('Вузли інтерполяції:\n');
fprintf('x = [%s]\n', num2str(x2));
fprintf('y = [%s]\n', num2str(y2));

% Обчислюємо поліном Лагранжа
[y_interp2, L2] = lagrange_interpolation(x2, y2, x2);

% Виводимо коефіцієнти
fprintf('\nКоефіцієнти поліному Лагранжа (від старшого до молодшого степеня):\n');
format rational;
disp(L2);
format;

% Візуалізація
figure;
plot(x2, y2, 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли інтерполяції');
hold on;
x_interp = linspace(min(x2), max(x2), 100);
y_interp = polyval(L2, x_interp);
plot(x_interp, y_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Поліном Лагранжа');
grid on;
legend('Location', 'best');
title('Інтерполяція Лагранжа для Прикладу 2 (f(x) = x^2)');
xlabel('x');
ylabel('y');

%% Приклад 3: Тригонометрична функція
fprintf('\nПриклад 3:\n');
x3 = [0, pi/4, pi/2, 3*pi/4];
y3 = [0, 1/sqrt(2), 1, 1/sqrt(2)];  % f(x) = sin(x)
fprintf('Вузли інтерполяції:\n');
fprintf('x = [%s]\n', num2str(x3));
fprintf('y = [%s]\n', num2str(y3));

% Обчислюємо поліном Лагранжа
[y_interp3, L3] = lagrange_interpolation(x3, y3, x3);

% Виводимо коефіцієнти
fprintf('\nКоефіцієнти поліному Лагранжа (від старшого до молодшого степеня):\n');
format rational;
disp(L3);
format;

% Візуалізація
figure;
plot(x3, y3, 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли інтерполяції');
hold on;
x_interp = linspace(min(x3), max(x3), 100);
y_interp = polyval(L3, x_interp);
plot(x_interp, y_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Поліном Лагранжа');
grid on;
legend('Location', 'best');
title('Інтерполяція Лагранжа для Прикладу 3 (f(x) = sin(x))');
xlabel('x');
ylabel('y');

% Видаляємо шлях до директорії після завершення
rmpath('../src/interpolation');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 