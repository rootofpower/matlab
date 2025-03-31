%% ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ
% Цей скрипт демонструє роботу різних методів чисельного інтегрування
% та порівнює їх ефективність.

close all;
clear;
clc;

fprintf('====== ТЕСТУВАННЯ МЕТОДІВ ЧИСЕЛЬНОГО ІНТЕГРУВАННЯ ======\n\n');

%% 1. ВИЗНАЧЕННЯ ТЕСТОВИХ ФУНКЦІЙ
fprintf('Тестові функції:\n');
fprintf('--------------------------------------\n');

% Тестова функція 1: f(x) = sin(x), інтеграл від 0 до pi = 2
f1 = @(x) sin(x);
a1 = 0;
b1 = pi;
exact1 = 2;
fprintf('1) f(x) = sin(x), [a,b] = [0,pi], точне значення = %g\n', exact1);

% Тестова функція 2: f(x) = x^2, інтеграл від 0 до 1 = 1/3
f2 = @(x) x.^2;
a2 = 0;
b2 = 1;
exact2 = 1/3;
fprintf('2) f(x) = x^2, [a,b] = [0,1], точне значення = %g\n', exact2);

% Тестова функція 3: f(x) = 1/x, інтеграл від 1 до 2 = ln(2)
f3 = @(x) 1./x;
a3 = 1;
b3 = 2;
exact3 = log(2);
fprintf('3) f(x) = 1/x, [a,b] = [1,2], точне значення = %g\n\n', exact3);

%% 2. ТЕСТУВАННЯ НА РІЗНІЙ КІЛЬКОСТІ ПІДІНТЕРВАЛІВ
% Кількість підінтервалів для тестування
n_values = [10, 100, 1000];

% Функція для перевірки методів на всіх тестових функціях
function test_all_functions(method_name, method_function, f1, a1, b1, exact1, f2, a2, b2, exact2, f3, a3, b3, exact3, n_values)
    fprintf('\n%s:\n', method_name);
    fprintf('-------------------------------------------------------------------\n');
    fprintf('Функція | n      | Наближене значення | Точне значення | Похибка\n');
    fprintf('-------------------------------------------------------------------\n');
    
    % Тестова функція 1
    for i = 1:length(n_values)
        n = n_values(i);
        
        % Для методу Сімпсона n має бути парним
        if strcmp(method_name, 'Метод Сімпсона') && mod(n, 2) ~= 0
            n = n + 1;
        end
        
        if strcmp(method_name, 'Метод прямокутників')
            % Для методу прямокутників тестуємо всі три правила
            [I_left, ~, ~] = rectangle_method(f1, a1, b1, n, 'left');
            [I_right, ~, ~] = rectangle_method(f1, a1, b1, n, 'right');
            [I_mid, ~, ~] = rectangle_method(f1, a1, b1, n, 'midpoint');
            
            fprintf('f1(x)  | %5d  | %.10f (left)  | %.10f     | %.3e\n', n, I_left, exact1, abs(exact1-I_left));
            fprintf('f1(x)  | %5d  | %.10f (right) | %.10f     | %.3e\n', n, I_right, exact1, abs(exact1-I_right));
            fprintf('f1(x)  | %5d  | %.10f (mid)   | %.10f     | %.3e\n', n, I_mid, exact1, abs(exact1-I_mid));
        else
            % Для інших методів
            [I, ~, ~] = method_function(f1, a1, b1, n);
            fprintf('f1(x)  | %5d  | %.10f        | %.10f     | %.3e\n', n, I, exact1, abs(exact1-I));
        end
    end
    
    % Тестова функція 2
    for i = 1:length(n_values)
        n = n_values(i);
        
        % Для методу Сімпсона n має бути парним
        if strcmp(method_name, 'Метод Сімпсона') && mod(n, 2) ~= 0
            n = n + 1;
        end
        
        if strcmp(method_name, 'Метод прямокутників')
            % Для методу прямокутників тестуємо всі три правила
            [I_left, ~, ~] = rectangle_method(f2, a2, b2, n, 'left');
            [I_right, ~, ~] = rectangle_method(f2, a2, b2, n, 'right');
            [I_mid, ~, ~] = rectangle_method(f2, a2, b2, n, 'midpoint');
            
            fprintf('f2(x)  | %5d  | %.10f (left)  | %.10f     | %.3e\n', n, I_left, exact2, abs(exact2-I_left));
            fprintf('f2(x)  | %5d  | %.10f (right) | %.10f     | %.3e\n', n, I_right, exact2, abs(exact2-I_right));
            fprintf('f2(x)  | %5d  | %.10f (mid)   | %.10f     | %.3e\n', n, I_mid, exact2, abs(exact2-I_mid));
        else
            % Для інших методів
            [I, ~, ~] = method_function(f2, a2, b2, n);
            fprintf('f2(x)  | %5d  | %.10f        | %.10f     | %.3e\n', n, I, exact2, abs(exact2-I));
        end
    end
    
    % Тестова функція 3
    for i = 1:length(n_values)
        n = n_values(i);
        
        % Для методу Сімпсона n має бути парним
        if strcmp(method_name, 'Метод Сімпсона') && mod(n, 2) ~= 0
            n = n + 1;
        end
        
        if strcmp(method_name, 'Метод прямокутників')
            % Для методу прямокутників тестуємо всі три правила
            [I_left, ~, ~] = rectangle_method(f3, a3, b3, n, 'left');
            [I_right, ~, ~] = rectangle_method(f3, a3, b3, n, 'right');
            [I_mid, ~, ~] = rectangle_method(f3, a3, b3, n, 'midpoint');
            
            fprintf('f3(x)  | %5d  | %.10f (left)  | %.10f     | %.3e\n', n, I_left, exact3, abs(exact3-I_left));
            fprintf('f3(x)  | %5d  | %.10f (right) | %.10f     | %.3e\n', n, I_right, exact3, abs(exact3-I_right));
            fprintf('f3(x)  | %5d  | %.10f (mid)   | %.10f     | %.3e\n', n, I_mid, exact3, abs(exact3-I_mid));
        else
            % Для інших методів
            [I, ~, ~] = method_function(f3, a3, b3, n);
            fprintf('f3(x)  | %5d  | %.10f        | %.10f     | %.3e\n', n, I, exact3, abs(exact3-I));
        end
    end
end

% Тестування всіх методів
% Метод прямокутників
test_all_functions('Метод прямокутників', @rectangle_method, f1, a1, b1, exact1, f2, a2, b2, exact2, f3, a3, b3, exact3, n_values);

% Метод трапецій
test_all_functions('Метод трапецій', @trapezoidal_method, f1, a1, b1, exact1, f2, a2, b2, exact2, f3, a3, b3, exact3, n_values);

% Метод Сімпсона
test_all_functions('Метод Сімпсона', @simpson_method, f1, a1, b1, exact1, f2, a2, b2, exact2, f3, a3, b3, exact3, n_values);

%% 3. ПОРІВНЯННЯ ШВИДКОСТІ ЗБІЖНОСТІ МЕТОДІВ
fprintf('\n\nПорівняння швидкості збіжності методів:\n');
fprintf('------------------------------------------\n');

% Функція для тесту збіжності
f_test = @(x) sin(x);
a_test = 0;
b_test = pi;
exact_test = 2;

% Кількість підінтервалів для тестування збіжності
n_converge = [4, 8, 16, 32, 64, 128, 256, 512, 1024];
errors_rect_mid = zeros(size(n_converge));
errors_trap = zeros(size(n_converge));
errors_simp = zeros(size(n_converge));

% Обчислення похибок для різної кількості підінтервалів
for i = 1:length(n_converge)
    n = n_converge(i);
    
    % Забезпечуємо парність n для методу Сімпсона
    if mod(n, 2) ~= 0
        n = n + 1;
    end
    
    % Метод прямокутників (середні точки)
    [I_rect, ~, ~] = rectangle_method(f_test, a_test, b_test, n, 'midpoint');
    errors_rect_mid(i) = abs(exact_test - I_rect);
    
    % Метод трапецій
    [I_trap, ~, ~] = trapezoidal_method(f_test, a_test, b_test, n);
    errors_trap(i) = abs(exact_test - I_trap);
    
    % Метод Сімпсона
    [I_simp, ~, ~] = simpson_method(f_test, a_test, b_test, n);
    errors_simp(i) = abs(exact_test - I_simp);
end

% Візуалізація збіжності
figure('Name', 'Порівняння швидкості збіжності методів', 'NumberTitle', 'off');
loglog(n_converge, errors_rect_mid, 'r-o', 'LineWidth', 2, 'DisplayName', 'Прямокутники (середні)');
hold on;
loglog(n_converge, errors_trap, 'g-s', 'LineWidth', 2, 'DisplayName', 'Трапеції');
loglog(n_converge, errors_simp, 'b-d', 'LineWidth', 2, 'DisplayName', 'Сімпсон');

% Теоретичні лінії збіжності
ref_lines = n_converge;
loglog(n_converge, 1./(ref_lines.^2), 'k--', 'LineWidth', 1, 'DisplayName', 'O(h^2)');
loglog(n_converge, 1./(ref_lines.^4), 'k-.', 'LineWidth', 1, 'DisplayName', 'O(h^4)');

grid on;
xlabel('Кількість підінтервалів (n)');
ylabel('Абсолютна похибка');
title('Швидкість збіжності різних методів чисельного інтегрування');
legend('Location', 'southwest');

% Таблиця з порівнянням швидкості збіжності
fprintf('Порядок збіжності для різних методів:\n');
fprintf('----------------------------------------------\n');
fprintf('Метод                 | Теоретичний порядок  \n');
fprintf('----------------------------------------------\n');
fprintf('Прямокутники (ліві)   | O(h)   = O(1/n)      \n');
fprintf('Прямокутники (праві)  | O(h)   = O(1/n)      \n');
fprintf('Прямокутники (середні)| O(h^2) = O(1/n^2)    \n');
fprintf('Трапеції              | O(h^2) = O(1/n^2)    \n');
fprintf('Сімпсона              | O(h^4) = O(1/n^4)    \n');
fprintf('----------------------------------------------\n');

%% 4. ВІЗУАЛІЗАЦІЯ МЕТОДІВ
fprintf('\nВізуалізація методів для функції f(x) = sin(x) на [0,pi]:\n');

n_viz = 10; % Кількість підінтервалів для візуалізації

% Візуалізація методу прямокутників
figure('Name', 'Метод прямокутників', 'NumberTitle', 'off');
subplot(1, 3, 1);
visualize_rectangle(f1, a1, b1, n_viz, 'left');
title('Метод лівих прямокутників');

subplot(1, 3, 2);
visualize_rectangle(f1, a1, b1, n_viz, 'right');
title('Метод правих прямокутників');

subplot(1, 3, 3);
visualize_rectangle(f1, a1, b1, n_viz, 'midpoint');
title('Метод середніх прямокутників');

% Візуалізація методу трапецій
figure('Name', 'Метод трапецій', 'NumberTitle', 'off');
visualize_trapezoidal(f1, a1, b1, n_viz);

% Візуалізація методу Сімпсона
figure('Name', 'Метод Сімпсона', 'NumberTitle', 'off');
if mod(n_viz, 2) ~= 0
    n_viz = n_viz + 1; % Забезпечуємо парність n для методу Сімпсона
end
visualize_simpson(f1, a1, b1, n_viz);

% Порівняння всіх методів
figure('Name', 'Порівняння методів (n=' + string(n_viz) + ')', 'NumberTitle', 'off');
x = linspace(a1, b1, 1000);
y = f1(x);
plot(x, y, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x) = sin(x)');
hold on;

% Обчислення значень інтегралів
[I_rect_mid, ~, ~] = rectangle_method(f1, a1, b1, n_viz, 'midpoint');
[I_trap, ~, ~] = trapezoidal_method(f1, a1, b1, n_viz);
[I_simp, ~, ~] = simpson_method(f1, a1, b1, n_viz);

% Додавання тексту з результатами
text_x = a1 + 0.7*(b1-a1);
text_y = 0.5;
str = {
    ['Точне значення = ' num2str(exact1, '%.10f')]
    ['Прямокутники (середні) = ' num2str(I_rect_mid, '%.10f') ', похибка = ' num2str(abs(exact1-I_rect_mid), '%.3e')]
    ['Трапеції = ' num2str(I_trap, '%.10f') ', похибка = ' num2str(abs(exact1-I_trap), '%.3e')]
    ['Сімпсон = ' num2str(I_simp, '%.10f') ', похибка = ' num2str(abs(exact1-I_simp), '%.3e')]
};
text(text_x, text_y, str, 'FontSize', 10, 'VerticalAlignment', 'middle');

grid on;
title('Порівняння всіх методів чисельного інтегрування');
xlabel('x');
ylabel('f(x)');

fprintf('\n====== ЗАВЕРШЕННЯ ТЕСТУВАННЯ ======\n');

%% Допоміжні функції для візуалізації
function visualize_rectangle(f, a, b, n, rule)
    [integral, x_points, y_points] = rectangle_method(f, a, b, n, rule);
    
    % Створення точного графіка функції
    x_fine = linspace(a, b, 1000);
    y_fine = f(x_fine);
    
    % Графік функції
    plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
    hold on;
    
    % Визначення точок для побудови прямокутників
    h = (b - a) / n;
    
    for i = 1:n
        x_left = a + (i-1) * h;
        x_right = a + i * h;
        
        % Визначення висоти прямокутника в залежності від правила
        switch lower(rule)
            case 'left'
                height = f(x_left);
            case 'right'
                height = f(x_right);
            case 'midpoint'
                height = f((x_left + x_right)/2);
            otherwise
                error('Невідоме правило');
        end
        
        % Малюємо прямокутник
        rectangle('Position', [x_left, 0, h, height], 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', [1 0 0 0.3]);
    end
    
    % Додавання підписів та сітки
    grid on;
    xlabel('x');
    ylabel('f(x)');
    title_str = sprintf('Метод прямокутників (%s), n = %d, I ≈ %.6f', rule, n, integral);
    title(title_str);
end

function visualize_trapezoidal(f, a, b, n)
    [integral, x_points, y_points] = trapezoidal_method(f, a, b, n);
    
    % Створення точного графіка функції
    x_fine = linspace(a, b, 1000);
    y_fine = f(x_fine);
    
    % Графік функції
    plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
    hold on;
    
    % Побудова трапецій
    for i = 1:n
        % Координати вершин трапеції
        x_vertices = [x_points(i), x_points(i), x_points(i+1), x_points(i+1)];
        y_vertices = [0, y_points(i), y_points(i+1), 0];
        
        % Малюємо трапецію
        fill(x_vertices, y_vertices, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
    end
    
    % Точки розбиття
    plot(x_points, y_points, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    % Додавання підписів та сітки
    grid on;
    xlabel('x');
    ylabel('f(x)');
    title_str = sprintf('Метод трапецій, n = %d, I ≈ %.6f', n, integral);
    title(title_str);
end

function visualize_simpson(f, a, b, n)
    [integral, x_points, y_points] = simpson_method(f, a, b, n);
    
    % Створення точного графіка функції
    x_fine = linspace(a, b, 1000);
    y_fine = f(x_fine);
    
    % Графік функції
    plot(x_fine, y_fine, 'b-', 'LineWidth', 2);
    hold on;
    
    % Побудова параболічних сегментів
    for i = 1:2:n-1
        % Створення точок для побудови параболи
        x_segment = linspace(x_points(i), x_points(i+2), 100);
        
        % Коефіцієнти інтерполяційного полінома (параболи)
        x0 = x_points(i);
        x1 = x_points(i+1);
        x2 = x_points(i+2);
        y0 = y_points(i);
        y1 = y_points(i+1);
        y2 = y_points(i+2);
        
        % Лагранжева інтерполяція для трьох точок
        L0 = @(x) ((x-x1).*(x-x2))/((x0-x1)*(x0-x2));
        L1 = @(x) ((x-x0).*(x-x2))/((x1-x0)*(x1-x2));
        L2 = @(x) ((x-x0).*(x-x1))/((x2-x0)*(x2-x1));
        
        % Інтерполяційний поліном
        p = @(x) y0*L0(x) + y1*L1(x) + y2*L2(x);
        
        % Значення полінома у точках x_segment
        y_segment = p(x_segment);
        
        % Малюємо параболічний сегмент
        plot(x_segment, y_segment, 'r-', 'LineWidth', 1.5);
        
        % Заповнення площі під параболою
        x_fill = [x_segment, fliplr(x_segment)];
        y_fill = [y_segment, zeros(size(y_segment))];
        fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Точки розбиття
    plot(x_points, y_points, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    % Додавання підписів та сітки
    grid on;
    xlabel('x');
    ylabel('f(x)');
    title_str = sprintf('Метод Сімпсона, n = %d, I ≈ %.6f', n, integral);
    title(title_str);
end 