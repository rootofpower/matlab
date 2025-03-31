function [integral, x_points, y_points] = rectangle_method(f, a, b, n, rule)
% RECTANGLE_METHOD Обчислення визначеного інтегралу методом прямокутників
%
% Опис методу:
% Метод прямокутників (обдĺžniková metóda) - це простий метод чисельного 
% інтегрування, який обчислює наближене значення визначеного інтегралу шляхом
% розбиття інтервалу інтегрування на n однакових частин і обчислення суми площ
% прямокутників. Залежно від вибору точки обчислення підінтегральної функції,
% розрізняють три правила:
% - 'left' - метод лівих прямокутників (значення функції обчислюються в лівих кінцях)
% - 'right' - метод правих прямокутників (значення функції обчислюються в правих кінцях)
% - 'midpoint' - метод середніх прямокутників (значення функції обчислюються в середині)
%
% Вхідні параметри:
%   f     - функція, яку потрібно проінтегрувати
%   a, b  - межі інтегрування
%   n     - кількість підінтервалів (прямокутників)
%   rule  - правило обчислення ('left', 'right', 'midpoint'), за замовчуванням 'midpoint'
%
% Вихідні параметри:
%   integral - обчислене значення інтегралу
%   x_points - точки розбиття (для візуалізації)
%   y_points - значення функції в точках розбиття (для візуалізації)

% Перевірка вхідних параметрів
if nargin < 5
    rule = 'midpoint';
end

% Визначення ширини підінтервалу
h = (b - a) / n;

% Ініціалізація масивів для точок розбиття та значень функції
x_points = zeros(1, n);
y_points = zeros(1, n);

% Обчислення інтегралу в залежності від обраного правила
switch lower(rule)
    case 'left'
        % Метод лівих прямокутників
        for i = 1:n
            x_i = a + (i-1) * h;
            x_points(i) = x_i;
            y_points(i) = f(x_i);
        end
        integral = h * sum(y_points);
        
    case 'right'
        % Метод правих прямокутників
        for i = 1:n
            x_i = a + i * h;
            x_points(i) = x_i;
            y_points(i) = f(x_i);
        end
        integral = h * sum(y_points);
        
    case 'midpoint'
        % Метод середніх прямокутників
        for i = 1:n
            x_i = a + (i-0.5) * h;
            x_points(i) = x_i;
            y_points(i) = f(x_i);
        end
        integral = h * sum(y_points);
        
    otherwise
        error('Невідоме правило інтегрування. Використовуйте "left", "right" або "midpoint".');
end

end

% Приклад використання:
% % Обчислення інтегралу sin(x) від 0 до pi (точне значення = 2)
% f = @(x) sin(x);
% a = 0;
% b = pi;
% n = 100;
% 
% % Обчислення за різними правилами
% [I_left, x_left, y_left] = rectangle_method(f, a, b, n, 'left');
% [I_right, x_right, y_right] = rectangle_method(f, a, b, n, 'right');
% [I_mid, x_mid, y_mid] = rectangle_method(f, a, b, n, 'midpoint');
% 
% fprintf('Метод лівих прямокутників: %.10f, похибка: %.10f\n', I_left, abs(2-I_left));
% fprintf('Метод правих прямокутників: %.10f, похибка: %.10f\n', I_right, abs(2-I_right));
% fprintf('Метод середніх прямокутників: %.10f, похибка: %.10f\n', I_mid, abs(2-I_mid));
% 
% % Візуалізація
% x_fine = linspace(a, b, 1000);
% y_fine = f(x_fine);
% 
% figure;
% plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
% hold on;
% 
% % Побудова прямокутників (наприклад, для методу середніх точок)
% for i = 1:n
%     x_left = a + (i-1) * (b-a)/n;
%     x_right = a + i * (b-a)/n;
%     y_mid = f((x_left + x_right)/2);
%     
%     % Малюємо прямокутник
%     rectangle('Position', [x_left, 0, (b-a)/n, y_mid], 'EdgeColor', 'r', 'LineWidth', 0.5);
% end
% 
% title('Метод прямокутників');
% xlabel('x');
% ylabel('f(x)');
% grid on; 