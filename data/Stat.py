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
    df = pd.read_csv(file_path, comment='#', sep='\s+' , names=col_names, header=None)
    return df

# Применяем, а ещё добавляем инфу о том, вирусный или нет
viral_df = parser("Viral_results.txt")
nonviral_df = parser("NonViral_results.txt")
viral_df["is_viral"] = 1
nonviral_df["is_viral"] = 0


# Зададим график
figure, axis = plt.subplots(2, 2)

# Строим гистограмму для full_score
axis[0, 0].hist([viral_df["full_score"], nonviral_df["full_score"]], 
        bins=100,
        alpha=0.5,
        label=['Viral', 'NonViral'],
        color=['0bb4ff', '50e991])
axis[0, 0].yscale('log')
axis[0, 0].title('Distribution of full_score in Viral sequences and NonViral sequences')
axis[0, 0].legend(['Viral', 'NonViral'])
axis[0, 0].ylabel('Number of sequences')
axis[0, 0].xlabel('Full score')
axis[0, 0].grid(True, alpha=0.3)


# Строим гистограмму для full_evalue
axis[0, 1].hist([viral_df["full_evalue"], nonviral_df["full_evalue"]], 
        bins=100,
        alpha=0.5,
        label=['Viral', 'NonViral'],
        color=['#0bb4ff', '#50e991])
axis[0, 1].yscale('log')
axis[0, 1].title('Distribution of full_evalue in Viral sequences and NonViral sequences')
axis[0, 1].legend(['Viral', 'NonViral'])
axis[0, 1].set_xscale('log')
axis[0, 1].ylabel('Number of sequences')
axis[0, 1].xlabel('Full evalue')
axis[0, 1].grid(True, alpha=0.3)


# Объединим датафреймы
combined_df = pd.concat([viral_df, nonviral_df], ignore_index=True)

# Строим ROC-кривую для full_score
scores = combined_df['full_score']
labels = combined_df['is_viral']
fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

axis[1, 0].plot(fpr, tpr, color='#0bb4ff', linewidth=2, 
                label=f'full_score (AUC = {roc_auc:.3f})')
axis[1, 0].plot([0, 1], [0, 1], 'k--', alpha=0.5)
axis[1, 0].xlabel('False Positive Rate')
axis[1, 0].ylabel('True Positive Rate')
axis[1, 0].title('ROC Curve for Full Score')
axis[1, 0].legend()
axis[1, 0].grid(True, alpha=0.3)

#Строим ROC-кривую для e_value
scores = -np.log10(combined_df['full_evalue'] + 1e-300) #Чтобы рост e_value соотствовал росту вероятности
labels = combined_df['is_viral']
fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

axis[1, 1].plot(fpr, tpr, color='#0bb4ff', linewidth=2, 
                label=f'full_evalue (AUC = {roc_auc:.3f})')
axis[1, 1].plot([0, 1], [0, 1], 'k--', alpha=0.5)
axis[1, 1].xlabel('False Positive Rate')
axis[1, 1].ylabel('True Positive Rate')
axis[1, 1].title('ROC Curve for Full Score')
axis[1, 1].legend()
axis[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
plt.savefig('Results.png')