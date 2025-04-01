function data = read_input_data(filename)
% READ_INPUT_DATA зчитує вхідні дані з текстового файлу
%   data = READ_INPUT_DATA(filename) - зчитує дані з файлу filename і
%   повертає структуру data з полями, які відповідають заголовкам у файлі.
%
%   Формат файлу:
%   - Першій рядок містить заголовки параметрів, розділені комами
%   - Наступні рядки містять значення параметрів
%   - Рядки, що починаються з #, ігноруються (коментарі)
%   - Пусті рядки ігноруються
%   - Спеціальні секції (як-от матриці) починаються з одного параметра в рядку
%     і продовжуються наступними рядками з даними
%
%   Приклад:
%   a, b, c, eps
%   1, 2, 3, 0.001
%   4, 5, 6, 0.0001
%   
%   # Матриця A
%   A
%   1, 2, 3
%   4, 5, 6
%   7, 8, 9

    % Відкриваємо файл для читання
    fid = fopen(filename, 'r');
    if fid == -1
        error('Не вдалося відкрити файл: %s', filename);
    end
    
    % Ініціалізуємо структуру для зберігання даних
    data = struct();
    
    % Зчитуємо файл рядок за рядком
    line = fgetl(fid);
    if ~ischar(line)
        fclose(fid);
        error('Файл порожній або невірний формат');
    end
    
    % Перевіряємо, чи є заголовки у першому рядку
    if ~startsWith(line, '#') && ~isempty(strtrim(line))
        % Розбираємо заголовки
        headers = strtrim(strsplit(line, ','));
        for i = 1:length(headers)
            headers{i} = strtrim(headers{i});
        end
        
        % Зчитуємо дані для цих заголовків
        while ischar(line)
            line = fgetl(fid);
            if ~ischar(line)
                break;
            end
            
            % Пропускаємо коментарі та порожні рядки
            if startsWith(line, '#') || isempty(strtrim(line))
                continue;
            end
            
            % Перевіряємо, чи це новий розділ (один параметр у рядку)
            parts = strtrim(strsplit(line, ','));
            if length(parts) == 1 && ~contains(line, ',')
                % Це новий розділ, зчитуємо дані для нього
                section_name = strtrim(parts{1});
                section_data = [];
                
                % Зчитуємо дані для цього розділу
                section_line = fgetl(fid);
                while ischar(section_line) && ~startsWith(section_line, '#') && ~isempty(strtrim(section_line))
                    % Розбираємо рядок на числа
                    section_values = strtrim(strsplit(section_line, ','));
                    section_values = cellfun(@(x) str2double(strtrim(x)), section_values);
                    
                    % Додаємо рядок до матриці даних
                    section_data = [section_data; section_values];
                    
                    % Зчитуємо наступний рядок
                    section_line = fgetl(fid);
                end
                
                % Зберігаємо дані розділу у структурі
                data.(section_name) = section_data;
                
                % Якщо досягли кінця файлу, виходимо
                if ~ischar(section_line)
                    break;
                end
                
                % Повертаємо вказівник на рядок назад, якщо це коментар або нова секція
                if startsWith(section_line, '#') || ~isempty(strtrim(section_line))
                    line = section_line;
                    continue;
                end
            else
                % Це рядок з даними для заголовків
                values = strtrim(strsplit(line, ','));
                if length(values) ~= length(headers)
                    warning('Кількість значень не відповідає кількості заголовків у рядку: %s', line);
                    continue;
                end
                
                % Перетворюємо рядки на числа, якщо можливо
                for i = 1:length(values)
                    value = strtrim(values{i});
                    if isStringNumeric(value)
                        values{i} = str2double(value);
                    end
                end
                
                % Додаємо значення до відповідних полів структури
                for i = 1:length(headers)
                    if isfield(data, headers{i})
                        % Якщо поле вже існує, додаємо до масиву
                        if iscell(data.(headers{i}))
                            data.(headers{i}){end+1} = values{i};
                        else
                            data.(headers{i}) = [data.(headers{i}); values{i}];
                        end
                    else
                        % Створюємо нове поле
                        data.(headers{i}) = values{i};
                    end
                end
            end
        end
    end
    
    % Закриваємо файл
    fclose(fid);
end

function result = isStringNumeric(str)
    % Перевіряє, чи є рядок числом
    result = ~isempty(str) && ~isnan(str2double(str));
end 