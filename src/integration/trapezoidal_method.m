function [integral, x_points, y_points] = trapezoidal_method(f, a, b, n)
% TRAPEZOIDAL_METHOD Обчислення визначеного інтегралу методом трапецій
%
% Опис методу:
% Метод трапецій (lichobežníková metóda) - це метод чисельного інтегрування,
% який обчислює наближене значення визначеного інтегралу шляхом розбиття
% інтервалу інтегрування на n однакових частин і апроксимації підінтегральної
% функції лінійною функцією на кожному підінтервалі. Площа під кривою 
% апроксимується сумою площ трапецій. 
% Цей метод є точнішим за метод прямокутників для більшості функцій.
%
% Вхідні параметри:
%   f     - функція, яку потрібно проінтегрувати
%   a, b  - межі інтегрування
%   n     - кількість підінтервалів (трапецій)
%
% Вихідні параметри:
%   integral - обчислене значення інтегралу
%   x_points - точки розбиття (для візуалізації)
%   y_points - значення функції в точках розбиття (для візуалізації)

% Визначення ширини підінтервалу
h = (b - a) / n;

% Створення точок розбиття
x_points = a:h:b;
y_points = zeros(1, n+1);

% Обчислення значень функції у всіх точках розбиття
for i = 1:n+1
    y_points(i) = f(x_points(i));
end

% Формула методу трапецій:
% I ≈ h/2 * [f(a) + 2*f(x_1) + 2*f(x_2) + ... + 2*f(x_{n-1}) + f(b)]
integral = h/2 * (y_points(1) + y_points(n+1) + 2*sum(y_points(2:n)));

end

% Приклад використання:
% % Обчислення інтегралу sin(x) від 0 до pi (точне значення = 2)
% f = @(x) sin(x);
% a = 0;
% b = pi;
% 
% % Обчислення для різної кількості підінтервалів
% n_values = [10, 100, 1000];
% 
% for i = 1:length(n_values)
%     n = n_values(i);
%     [I_trap, x_points, y_points] = trapezoidal_method(f, a, b, n);
%     fprintf('n = %d: I ≈ %.10f, Похибка: %.10e\n', n, I_trap, abs(2-I_trap));
% end
% 
% % Візуалізація
% n = 10;  % Мала кількість для кращої візуалізації
% [I_trap, x_points, y_points] = trapezoidal_method(f, a, b, n);
% 
% figure;
% x_fine = linspace(a, b, 1000);
% y_fine = f(x_fine);
% 
% % Графік функції
% plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
% hold on;
% 
% % Побудова трапецій
% for i = 1:n
%     % Координати вершин трапеції
%     x_vertices = [x_points(i), x_points(i), x_points(i+1), x_points(i+1)];
%     y_vertices = [0, y_points(i), y_points(i+1), 0];
%     
%     % Малюємо трапецію
%     fill(x_vertices, y_vertices, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
% end
% 
% % Точки розбиття
% plot(x_points, y_points, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
% 
% title(['Метод трапецій, n = ', num2str(n), ', I ≈ ', num2str(I_trap, '%.6f')]);
% xlabel('x');
% ylabel('f(x)');
% grid on; 