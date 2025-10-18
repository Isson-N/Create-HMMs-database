import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, RocCurveDisplay
import matplotlib.pyplot as plt
import numpy as np


# Функция, чтобы парсить таблицы, полученные из hmmsearch
def parser(file_path):
    col_names = ['target_name', 'target_accession', 'query_name', 'query_accession', 
                 'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                 'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
    df = pd.read_csv(file_path, comment='#', sep=r'\s+', names=col_names, header=None)
    return df

# Применяем, а ещё добавляем инфу о том, вирусный или нет
viral_df = parser("Viral_results.txt")
nonviral_df = parser("NonViral_results.txt")
combined_df = pd.concat([viral_df, nonviral_df], ignore_index=True)

