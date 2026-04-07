import pandas as pd

# читаем с индексом
hmm_info = pd.read_csv("hmm_info_updated_kobu.csv", index_col=0)

# гарантируем, что индекс — int
hmm_info.index = hmm_info.index.astype(int)

# берём следующий индекс
new_index = hmm_info.index.max() + 1

# новая строка
new_row = pd.DataFrame([{
    "Model_ID": "vHMM_21553",
    "Taxon": "Parechovirus",
    "Threshold": 12.0,
    "Positive terms": "polyprotein_custom"
}], index=[new_index])

# добавляем
hmm_info = pd.concat([hmm_info, new_row])

# сохраняем БЕЗ лишних преобразований
hmm_info.to_csv("hmm_info_updated_parecho.csv")