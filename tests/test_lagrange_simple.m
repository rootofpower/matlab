%% ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА
% Цей скрипт перевіряє роботу методу Лагранжа на простому прикладі

close all;
clear;
clc;

% Додаємо шлях до директорії з методами інтерполяції
addpath('../src/interpolation');

fprintf('====== ТЕСТУВАННЯ ІНТЕРПОЛЯЦІЙНОГО ПОЛІНОМУ ЛАГРАНЖА ======\n\n');

%% Тестовий приклад
x = [0, 1, 2, 3];
y = [0.5, 2.1, 8.2, 25.6];

fprintf('Вузли інтерполяції:\n');
fprintf('x = [%s]\n', num2str(x));
fprintf('y = [%s]\n', num2str(y));

% Обчислюємо поліном Лагранжа
[y_interp, L] = lagrange_interpolation(x, y, x);

% Виводимо коефіцієнти
fprintf('\nКоефіцієнти поліному Лагранжа (від старшого до молодшого степеня):\n');
format rational;
disp(L);
format;

% Перевіряємо значення в вузлах
fprintf('\nЗначення поліному в вузлах:\n');
disp(y_interp);

% Візуалізація
figure;
plot(x, y, 'ko', 'MarkerSize', 10, 'DisplayName', 'Вузли інтерполяції');
hold on;
x_interp = linspace(min(x), max(x), 100);
y_interp = polyval(L, x_interp);
plot(x_interp, y_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Поліном Лагранжа');
grid on;
legend('Location', 'best');
title('Інтерполяція Лагранжа');
xlabel('x');
ylabel('y');

% Видаляємо шлях до директорії після завершення
rmpath('../src/interpolation');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n'); 