function [y_interp, L] = lagrange_interpolation(x, y, x_interp)
    % LAGRANGE_INTERPOLATION Обчислює значення інтерполяційного поліному Лагранжа
    % 
    % Вхідні параметри:
    % x - масив x-координат вузлів інтерполяції
    % y - масив y-координат вузлів інтерполяції
    % x_interp - масив точок, в яких потрібно обчислити значення поліному
    %
    % Вихідні параметри:
    % y_interp - масив значень поліному в точках x_interp
    % L - матриця базисних поліномів Лагранжа (опціонально)
    
    % Перевірка вхідних даних
    if length(x) ~= length(y)
        error('Масиви x та y повинні мати однакову довжину');
    end
    
    n = length(x);
    m = length(x_interp);
    
    % Ініціалізація масиву для результатів
    y_interp = zeros(1, m);
    
    % Обчислення значень поліному в кожній точці x_interp
    for i = 1:m
        xi = x_interp(i);
        sum = 0;
        
        % Обчислення суми базисних поліномів
        for j = 1:n
            term = y(j);
            
            % Обчислення базисного поліному
            for k = 1:n
                if k ~= j
                    term = term * (xi - x(k)) / (x(j) - x(k));
                end
            end
            
            sum = sum + term;
        end
        
        y_interp(i) = sum;
    end
    
    % Якщо потрібно, обчислюємо базисні поліноми
    if nargout > 1
        L = zeros(n, m);
        for j = 1:n
            for i = 1:m
                xi = x_interp(i);
                term = 1;
                
                for k = 1:n
                    if k ~= j
                        term = term * (xi - x(k)) / (x(j) - x(k));
                    end
                end
                
                L(j, i) = term;
            end
        end
    end
end 