import pandas as pd
import matplotlib.pyplot as plt
import sys

file_path = sys.argv[1]
output_file_name = sys.argv[2]
df = pd.read_csv(file_path, sep='\t')

#revisar
categorias_relevantes = ['EnterovirusA', 'EnterovirusB', 'EnterovirusC', 'EnterovirusD', 'EnterovirusE']

#filtro
df['Genus'] = df['Genus'].apply(lambda x: x if x in categorias_relevantes else 'Others')
grouped_data_detailed = df.groupby(['Genus', 'Specie'])['Count'].sum().unstack(fill_value=0)
ax = grouped_data_detailed.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab10')

#no-legend
ax.legend().remove()

plt.xlabel('Genus')
plt.ylabel('Count')
plt.title('Counts of Different Species within Each Genus')
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(output_file_name)
#optional
#plt.show()
