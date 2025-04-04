function [integral, x_points, y_points] = simpson_method(f, a, b, n)
% SIMPSON_METHOD Обчислення визначеного інтегралу методом Сімпсона
%
% Опис методу:
% Метод Сімпсона (Simpsonova metóda) - це метод чисельного інтегрування,
% який обчислює наближене значення визначеного інтегралу шляхом розбиття
% інтервалу інтегрування на парну кількість n підінтервалів і апроксимації
% підінтегральної функції поліномом другого степеня (параболою) на кожній
% парі сусідніх підінтервалів.
% Цей метод, як правило, є точнішим за методи прямокутників і трапецій,
% особливо для функцій, які мають плавну поведінку.
%
% Вхідні параметри:
%   f     - функція, яку потрібно проінтегрувати
%   a, b  - межі інтегрування
%   n     - кількість підінтервалів (має бути парним числом)
%
% Вихідні параметри:
%   integral - обчислене значення інтегралу
%   x_points - точки розбиття (для візуалізації)
%   y_points - значення функції в точках розбиття (для візуалізації)

% Перевірка на парність кількості підінтервалів
if mod(n, 2) ~= 0
    error('Кількість підінтервалів n має бути парним числом.');
end

% Визначення ширини підінтервалу
h = (b - a) / n;

% Створення точок розбиття
x_points = a:h:b;
y_points = zeros(1, n+1);

% Обчислення значень функції у всіх точках розбиття
for i = 1:n+1
    y_points(i) = f(x_points(i));
end

% Формула методу Сімпсона:
% I ≈ h/3 * [f(a) + 4*f(x_1) + 2*f(x_2) + 4*f(x_3) + ... + 2*f(x_{n-2}) + 4*f(x_{n-1}) + f(b)]
% де x_i - це точки розбиття, що відповідають парним та непарним індексам

% Обчислення суми коефіцієнтів для внутрішніх точок
% Коефіцієнти 4 для точок з непарними індексами (x_1, x_3, ...)
% Коефіцієнти 2 для точок з парними індексами (x_2, x_4, ...)
sum_interior = 4 * sum(y_points(2:2:n)) + 2 * sum(y_points(3:2:n-1));

% Повна сума з крайніми точками
integral = h/3 * (y_points(1) + sum_interior + y_points(n+1));

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
%     if mod(n, 2) ~= 0
%         n = n + 1;  % Забезпечуємо парність n
%     end
%     [I_simp, x_points, y_points] = simpson_method(f, a, b, n);
%     fprintf('n = %d: I ≈ %.10f, Похибка: %.10e\n', n, I_simp, abs(2-I_simp));
% end
% 
% % Візуалізація
% n = 10;  % Мала кількість для кращої візуалізації
% if mod(n, 2) ~= 0
%     n = n + 1;
% end
% [I_simp, x_points, y_points] = simpson_method(f, a, b, n);
% 
% figure;
% x_fine = linspace(a, b, 1000);
% y_fine = f(x_fine);
% 
% % Графік функції
% plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
% hold on;
% 
% % Побудова параболічних сегментів
% for i = 1:2:n-1
%     % Створення точок для побудови параболи
%     x_segment = linspace(x_points(i), x_points(i+2), 100);
%     
%     % Коефіцієнти інтерполяційного полінома (параболи)
%     x0 = x_points(i);
%     x1 = x_points(i+1);
%     x2 = x_points(i+2);
%     y0 = y_points(i);
%     y1 = y_points(i+1);
%     y2 = y_points(i+2);
%     
%     % Лагранжева інтерполяція для трьох точок
%     L0 = @(x) ((x-x1).*(x-x2))/((x0-x1)*(x0-x2));
%     L1 = @(x) ((x-x0).*(x-x2))/((x1-x0)*(x1-x2));
%     L2 = @(x) ((x-x0).*(x-x1))/((x2-x0)*(x2-x1));
%     
%     % Інтерполяційний поліном
%     p = @(x) y0*L0(x) + y1*L1(x) + y2*L2(x);
%     
%     % Значення полінома у точках x_segment
%     y_segment = p(x_segment);
%     
%     % Малюємо параболічний сегмент
%     plot(x_segment, y_segment, 'r-', 'LineWidth', 1.5);
%     
%     % Заповнення площі під параболою
%     x_fill = [x_segment, fliplr(x_segment)];
%     y_fill = [y_segment, zeros(size(y_segment))];
%     fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% end
% 
% % Точки розбиття
% plot(x_points, y_points, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
% 
% title(['Метод Сімпсона, n = ', num2str(n), ', I ≈ ', num2str(I_simp, '%.6f')]);
% xlabel('x');
% ylabel('f(x)');
% grid on; 