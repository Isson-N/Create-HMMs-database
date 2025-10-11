import pandas as pd
import matplotlib.pyplot as plt


def parser(file_path):
    col_names = ['target_name', 'target_accession', 'query_name', 'query_accession', 
                 'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                 'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
    df = pd.read_csv(file_path, comment='#', delim_whitespace=True, names=col_names, header=None)
    return df

viral_df = parser("Viral_results.txt")
nonviral_df = parser("NonViral_results.txt")

viral_df["is_viral"] = 1
nonviral_df["is_viral"] = 0

plt.boxplot([viral_df["full_score"], nonviral_df["full_score"]], 
            labels=['Viral', 'NonViral'],
            patch_artist=True,
            boxprops=dict(facecolor='green', color='green'),
            medianprops=dict(color='black'))

plt.title('Distribution of HMMER Scores: Viral vs Non-Viral Proteins')
plt.ylabel('Full Score')
plt.xlabel('Protein Type')
plt.grid(True, alpha=0.3)

# Сохранение и показ
plt.savefig('score_distribution.png', bbox_inches='tight', dpi=150)
plt.show()