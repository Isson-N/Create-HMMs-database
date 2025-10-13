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
viral_df["is_viral"] = 1
nonviral_df["is_viral"] = 0


# Строим гистограмму для full_score
plt.figure(figsize=(10, 6))
plt.hist([viral_df["full_score"], nonviral_df["full_score"]],
        bins=100,
        alpha=0.5,
        label=['Viral', 'NonViral'],
        color=['#50e991', '#EE204D'])
plt.yscale('log')
plt.title('Distribution of full_score in Viral sequences and NonViral sequences')
plt.legend(['Viral', 'NonViral'])
plt.ylabel('Number of sequences')
plt.xlabel('Full score')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig('full_score_distribution.png')

# Строим гистограмму для full_evalue
plt.figure(figsize=(10, 6))
plt.hist([viral_df["full_evalue"], nonviral_df["full_evalue"]], 
        bins=100,
        alpha=0.5,
        label=['Viral', 'NonViral'],
        color=['#0bb4ff', '#EE204D'])
plt.yscale('log')
plt.title('Distribution of full_evalue in Viral sequences and NonViral sequences')
plt.legend(['Viral', 'NonViral'])
plt.ylabel('Number of sequences')
plt.xlabel('Full evalue')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig('full_evalue_distribution.png')


# Объединим датафреймы
combined_df = pd.concat([viral_df, nonviral_df], ignore_index=True)

# Строим ROC-кривую для full_score
plt.figure(figsize=(8, 8))
scores = combined_df['full_score']
labels = combined_df['is_viral']
fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

plt.plot(fpr, tpr, color='#0bb4ff', linewidth=2, 
                label=f'full_score (AUC = {roc_auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)

optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]
optimal_fpr = fpr[optimal_idx]
optimal_tpr = tpr[optimal_idx]
plt.scatter(optimal_fpr, optimal_tpr, color='red', s=100, zorder=5,
           label=f'Optimal threshold = {optimal_threshold:.3f}\n(TPR = {optimal_tpr:.3f}, FPR = {optimal_fpr:.3f})')

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Full Score')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig('full_score_ROC.png')


#Строим ROC-кривую для e_value
plt.figure(figsize=(8, 8))
scores = 1-combined_df['full_evalue'] #Чтобы рост e_value соотствовал росту вероятности
labels = combined_df['is_viral']
fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

plt.plot(fpr, tpr, color='#0bb4ff', linewidth=2, 
                label=f'full_evalue (AUC = {roc_auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)

optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]
optimal_fpr = fpr[optimal_idx]
optimal_tpr = tpr[optimal_idx]
plt.scatter(optimal_fpr, optimal_tpr, color='red', s=100, zorder=5,
           label=f'Optimal threshold = {optimal_threshold:.3f}\n(TPR = {optimal_tpr:.3f}, FPR = {optimal_fpr:.3f})')

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for (1 - Full_Evalue)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig('full_evalue_ROC.png')


