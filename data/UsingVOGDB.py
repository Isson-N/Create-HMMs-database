import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import argparse

# В VOGах были обнаружены повторяющиеся последовательности.
# Нужно по идее узнать откуда они там, потому что в статье говорится, что быть не должно.

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--api_key", required=True, help="Api_Key from NCBI")
parser.add_argument("-f", "--file", required=True, help="File to analyze")
args = parser.parse_args()

Entrez.api_key = args.api_key
file_name = args.file

# Функция, чтобы парсить странные айдишники из VOGDB. Тут удаляются дупликаты
def ids(name):
    try:
        new_name = f'{name.split('.')[1]}.{name.split('.')[2]}'
    except:
        new_name = ""
    return new_name
# Список со всеми айдишниками в нормальном виде
ident = []
seqs = []
for record in SeqIO.parse(file_name, "fasta"):
    if ids(record.id) != "":
        ident.append(ids(record.id))
        seqs.append(str(record.seq))
print(f'Извлёк ids {"\n".join(ident)}')

# Создаём список с фамилиями, чтобы потом закинуть в датафрейм
families = []
with Entrez.efetch(db = 'protein', rettype = 'gb', id = ident) as handle:
    gb = SeqIO.parse(handle, 'genbank')
    for record in gb:
        try:
            family = record.annotations.get('taxonomy')[6]
        except IndexError:
            family = 'Unrecognized'

        families.append(family)

df = pd.DataFrame({'Id': ident, 'Family': families, 'Sequence': seqs})
df.drop_duplicates(keep='first', inplace=True) # Удаление дубликатов, которых хз откуда

# Записываем инфу в файлы. Для каждого семейства должно быть минимум 5 прочтений
for family in df['Family'].unique():
    if len(df[df['Family'] == family]) >= 5:
        print(f'Создаю файл для семейства {family}')
        with open(f'{family.replace(" ", "_").replace(r"\", "")}_{file_name.replace('.faa', '')}.fasta', 'w') as file:
            for iter, row in df[df['Family'] == family].iterrows():
                file.write(f'>{row['Family']}_{row['Id']}')
                file.write('\n')
                file.write(row['Sequence'])
                file.write('\n')
    else:
        print(f'Пропускаю семейство {family}')
