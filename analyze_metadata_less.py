import pandas as pd
import re
from io import StringIO

# Загрузка данных
data = pd.read_csv("sequences (1).csv")

print(f"Всего строк в датафрейме: {len(data)}")

# Функция для извлечения основных типов вирусов
def extract_virus_type(name):
    if pd.isna(name):
        return "Unknown"
    
    name_str = str(name).lower()
    
    # Основные типы вирусов
    patterns = {
        'Picornavirus': r'picornavirus|picoV|pico-virus',
        'Hepatovirus': r'hepatovirus',
        'Kobuvirus': r'kobuvirus',
        'Mischivirus': r'mischivirus',
        'Shanbavirus': r'shanbavirus',
        'Sapelovirus': r'sapelovirus',
        'Parechovirus': r'parechovirus',
        'Teschovirus': r'teschovirus',
        'Cardiovirus': r'cardiovirus'
    }
    
    for virus_type, pattern in patterns.items():
        if re.search(pattern, name_str):
            return virus_type
    
    return "Other"

# Применяем функцию к данным
data['Virus_Type'] = data['Organism_Name'].apply(extract_virus_type)

# Подсчитываем распределение
virus_distribution = data['Virus_Type'].value_counts()

# Рассчитываем проценты
total_sequences = len(data)
virus_percentages = (virus_distribution / total_sequences * 100).round(2)

# Создаем итоговую таблицу
result_df = pd.DataFrame({
    'Count': virus_distribution,
    'Percentage': virus_percentages
})

print("\nРаспределение вирусов по типам:")
print("=" * 40)
print(f"Всего последовательностей: {total_sequences}")
print("=" * 40)

for virus_type, row in result_df.iterrows():
    print(f"{virus_type:<15} {row['Count']:>3} ({row['Percentage']:>5}%)")

print("=" * 40)
print(f"Сумма всех Count: {result_df['Count'].sum()}")
print(f"Всего типов вирусов: {len(result_df)}")