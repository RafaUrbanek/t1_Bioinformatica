import os
import sys
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Definindo os argumentos de linha de código
parser = argparse.ArgumentParser(description="Script para processamentos dos dados resultantes do alinhamento do blastn.")

# Define o argumento 'input'
parser.add_argument("input", help="Caminho para o arquivo .csv de entrada")
# Define o argumento 'output'
parser.add_argument("output", help="Caminho para os arquivos de saída")

args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_PATH = args.output

# Criando o diretório de resultados caso não exista
if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)

# Carregando os dados
cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
df = pd.read_csv(args.input, names=cols)

# Filtrando bitscore maior ou igual a 700
df = df[df["bitscore"] >= 700]

# Tratando o DataFrame

# Pegando os dados do id
split_qseqid = df['qseqid'].str.split('|', expand=True)
split_sseqid = df['sseqid'].str.split('|', expand=True)
# Adiconando novas colunas no Dataframe para o query
df['q_acc'] = split_qseqid[0] # Acesso
df['q_year'] = (split_qseqid[1].str.split('-').str[0]).astype(int)
df['q_variant'] = split_qseqid[2]
# Adiconando novas colunas no Dataframe para o subject
df['s_acc'] = split_sseqid[0]
df['s_year'] = (split_sseqid[1].str.split('-').str[0]).astype(int)
df['s_variant'] = split_sseqid[2]
# Retirando auto-comparações
df = df[df["q_acc"] != df["s_acc"]]
# Garantindo o incio e fim real 
df["real_sstart"] = df[["sstart","send"]].min(axis=1)
df["real_send"] = df[["sstart","send"]].max(axis=1)
df["real_qstart"] = df[["qstart","qend"]].min(axis=1)
df["real_qend"] = df[["qstart","qend"]].max(axis=1)
# Criando um ponto médio dos alinhamentos
df["smidpoint"] = ((df["real_sstart"] + df["real_send"]) / 2).astype(int)
df["qmidpoint"] = ((df["real_qstart"] + df["real_qend"]) / 2).astype(int)

# --- SIMILARIDADE ---

# Mapa de disperção
temp_pos_identity = df.groupby(['s_year', 'smidpoint'])['pident'].mean().reset_index()
plt.figure(figsize=(16, 6))
# Plotando uma série para cada ano
for year in sorted(temp_pos_identity['s_year'].unique()):
    subset = temp_pos_identity[temp_pos_identity['s_year'] == year]
    plt.scatter(subset['smidpoint'], subset['pident'], s=10, alpha=0.4, label=f'{year}')
plt.title('Similaridade (% Identidade) ao Longo do Genoma por Ano')
plt.xlabel('Posição genômica (nt)')
plt.ylabel('% Identidade Média')
# Ajustando escala do eixo X para cada 2500 nt e limites
plt.xticks(range(0, 30001, 2500))
plt.xlim(0, 30000)
plt.ylim(97, 100.1) # Foco na alta similaridade definida no BLAST
plt.legend(title="Ano Alvo", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Heatmap ano x Posição Genoma
pident_heatmap_df = df.groupby(['s_year', 'smidpoint'])['pident'].mean().unstack(level='s_year').T.fillna(0)
plt.figure(figsize=(18, 8))
# Aplicando o heatmap
sns.heatmap(pident_heatmap_df, cmap='viridis', cbar_kws={'label': 'Média de % Identidade'}, vmin=97, vmax=100) # Ajustar vmin/vmax para o intervalo de identidade
plt.title('Heatmap de Média de % Identidade por Ano e Posição Genômica')
plt.xlabel('Posição Genômica (nt)')
plt.ylabel('Ano Alvo')
# Ajustar os ticks do eixo X para melhor legibilidade
intervalo = 5000
colunas = pident_heatmap_df.columns.values
labels_ticks = np.arange(0, colunas.max() + intervalo, intervalo).astype(int)
posicoes_ticks = [np.argmin(np.abs(colunas - x)) for x in labels_ticks]
plt.xticks(posicoes_ticks, labels_ticks, rotation=45)
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Heatmap ano x ano
# Get all unique years from both query_year and target_year, filter out 'UNK', and sort numerically
all_years = sorted(list(set(df['q_year'].unique()) | set(df['s_year'].unique())))
all_years = [int(y) for y in all_years]
all_years.sort()
heatmap_data = df.pivot_table(
    values="pident",
    index="q_year",
    columns="s_year",
    aggfunc="mean"
)
# Reindex the heatmap_data to include all years on both axes, filling missing data with NaN
heatmap_data = heatmap_data.reindex(index=all_years, columns=all_years)
plt.figure(figsize=(10, 8)) # Aumentei levemente para acomodar a barra de cores lateral
sns.heatmap(
    heatmap_data, 
    annot=True, 
    fmt=".3f", 
    cmap='RdYlBu',
    vmin=99.5, 
    vmax=100, 
    linewidths=.5,   # Adiciona uma linha fina entre os quadrados para clareza
    cbar_kws={'label': '% Identidade Média'}
)
plt.title("Similaridade Temporal SARS-CoV-2", fontsize=14, pad=20)
plt.xlabel("Ano alvo")
plt.ylabel("Ano query")
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

plt.figure(figsize=(10, 8)) # Aumentei levemente para acomodar a barra de cores lateral
sns.heatmap(
    heatmap_data, 
    annot=True, 
    fmt=".3f", 
    cmap='RdYlBu_r', # Mais próximo do "Azul -> Vermelho" do seu título
    vmin=99.5, 
    vmax=100, 
    linewidths=.5,   # Adiciona uma linha fina entre os quadrados para clareza
    cbar_kws={'label': '% Identidade Média'}
)

# --- MISMATCHES ---

# Gráfico em Barras
temp_pos_mismatch = df.groupby(['s_year', 'smidpoint'])['mismatch'].mean().reset_index()
plt.figure(figsize=(16, 6))
for year in sorted(temp_pos_mismatch['s_year'].unique()):
    subset = temp_pos_mismatch[temp_pos_mismatch['s_year'] == year]
    plt.plot(subset['smidpoint'], subset['mismatch'], label=f'{year}', alpha=0.7)
plt.xlabel("Posição genômica (nt)")
plt.ylabel("Média de Mismatches")
plt.title("Comparativo de Mismatches para cada ano ao longo do genoma")
plt.xlim(0, 30000)
plt.xticks(range(0, 30001, 2500))
plt.legend(title="Ano", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Heatmap
plt.figure(figsize=(10, 7))
mismatch_pivot = df.pivot_table(values='mismatch', index='q_year', columns='s_year', aggfunc='mean')
sns.heatmap(mismatch_pivot, annot=True, fmt='.2f', cmap='viridis_r')
plt.title('Média de Mismatches (Temporal)')
plt.xlabel('Ano Alvo')
plt.ylabel('Ano Query')
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# --- GAPS ---

# Gaps por posição
temp_pos_gaps = df.groupby(['s_year', 'smidpoint'])['gaps'].mean().reset_index()
plt.figure(figsize=(16, 6))
for year in sorted(temp_pos_gaps['s_year'].unique()):
    subset = temp_pos_gaps[temp_pos_gaps['s_year'] == year]
    plt.scatter(subset['smidpoint'], subset['gaps'], s=10, alpha=0.5, label=f'{year}')
plt.title('GAPs ao Longo do Genoma')
plt.xlabel('Posição genômica (nt)')
plt.ylabel('Média de GAPs')
# Ajustando escala do eixo X para cada 2500 nt
plt.xticks(range(0, 30001, 2500))
plt.xlim(0, 30000)
plt.legend(title="Ano Alvo", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Gaps por ano
gaps_density_by_year = df.groupby(['s_year', 'smidpoint'])['gaps'].mean().reset_index()
plt.figure(figsize=(16, 6))
for year in sorted(gaps_density_by_year['s_year'].unique()):
    subset = gaps_density_by_year[gaps_density_by_year['s_year'] == year]
    plt.plot(subset['smidpoint'], subset['gaps'], label=f'{year}', alpha=0.7)
plt.xlabel("Posição genômica (nt)")
plt.ylabel("Média de GAPs")
plt.title("Comparativo entre as medias de GAPs para cada ano ao longo do genoma")
plt.xlim(0, 30000)
plt.xticks(range(0, 30001, 2500)) # Ajustando escala do eixo X para cada 2500 nt
plt.legend(title="Ano", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Heatmap temporal
plt.figure(figsize=(10, 7))
gaps_pivot = df.pivot_table(values='gaps', index='q_year', columns='s_year', aggfunc='mean')
sns.heatmap(gaps_pivot, annot=True, fmt='.3f', cmap='coolwarm')
plt.title('Média de GAPs (Temporal)')
plt.xlabel('Ano Alvo')
plt.ylabel('Ano Query')
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Inserções
bp_insertions = df[df['gaps'] > 0].groupby('smidpoint')['gaps'].count()
plt.figure(figsize=(12, 4))
plt.bar(bp_insertions.index, bp_insertions.values, color='teal', alpha=0.6)
plt.xlabel("Posição genômica (nt)")
plt.ylabel("Frequência de Inserções")
plt.title("Distribuição de Inserções por Posição")
plt.xlim(0, 30000)
plt.grid(True, axis='y', alpha=0.3)
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Deleções
gaps_density = df.groupby('smidpoint')['gaps'].mean().reset_index()
plt.figure(figsize=(12, 4))
plt.bar(gaps_density['smidpoint'], gaps_density['gaps'], color='salmon', alpha=0.7, label='Média de GAPs (Deleções)', width=50)
plt.xlabel("Posição genômica (nt)")
plt.ylabel("Média de GAPs")
plt.title("Distribuição de GAPs (Deleções) ao longo do Genoma")
plt.xlim(0, 30000)
plt.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.legend()
plt.tight_layout()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# --- DIVERGÊNCIA GENÉTICA ---

# Divergência temporal
distance_table = df.groupby(["q_year","s_year"])["mismatch"].mean().unstack()
div_table = df.groupby(["q_year","s_year"])["mismatch"].mean().reset_index()
div_table["div_kb"] = (div_table["mismatch"] / 250) * 1000
div_table["time_gap"] = abs(div_table["q_year"] - div_table["s_year"])
time_curve = div_table.groupby("time_gap")["div_kb"].mean().reset_index()
years = sorted(div_table["s_year"].unique())
plt.figure(figsize=(10,6))
for qyear in sorted(div_table["q_year"].unique()):

    subset = div_table[div_table["q_year"] == qyear]
    subset = subset.sort_values("s_year")

    plt.plot(subset["s_year"], subset["div_kb"], marker="o", label=f"{qyear}")
plt.xlabel("Ano alvo")
plt.ylabel("Substituições médias por kb")
plt.title("Tendência temporal da divergência genética SARS-CoV-2")
plt.legend()
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Acúmulo de divergência temporal
x = time_curve["time_gap"].astype(float)
y = time_curve["div_kb"].astype(float)
coef = np.polyfit(x, y, 1)
poly1d_fn = np.poly1d(coef)
y_pred = poly1d_fn(x)
plt.figure(figsize=(8,5))
plt.plot(x, y, marker='o', label='Dados observados (Acúmulo de Divergência)')
plt.plot(x, y_pred, linestyle='--', label=f'Regressão Linear (y = {coef[0]:.3f}x + {coef[1]:.3f})')
plt.xlabel('Distância temporal entre isolados (anos)')
plt.ylabel('Substituições médias por kb')
plt.title('Relógio Molecular: Acúmulo de Divergência e Correlação Temporal')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# --- Proteína Spike ---

# Heatmap de similaridade
spike = df[(df["real_sstart"] >= 21563) & (df["real_send"] <= 25384)]
spike_heatmap = spike.pivot_table(
    values="pident",
    index="q_year",
    columns="s_year",
    aggfunc="mean"
)
plt.figure(figsize=(10, 8))
sns.heatmap(spike_heatmap, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': '% Identidade Média'})
plt.title("Similaridade Temporal na Proteína Spike (nt 21563-25384)")
plt.xlabel("Ano Alvo")
plt.ylabel("Ano Query")
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')

# Mismatches
spike_mismatch_density = spike.groupby("smidpoint")["mismatch"].mean()
spike_start = 21563
spike_end = 25384
plt.figure(figsize=(12, 5))
plt.plot(spike_mismatch_density.index, spike_mismatch_density.values)
plt.xlabel("Posição genômica (nt)")
plt.ylabel("Média de mismatches")
plt.title("Densidade de Mismatches na Proteína Spike")
plt.xlim(spike_start, spike_end) # Limitar o eixo X à região da Spike
plt.grid(True, linestyle='--', alpha=0.7)
texto_do_titulo = plt.gca().get_title()
caminho_completo = os.path.join(OUTPUT_PATH, f"{texto_do_titulo}.png")
plt.savefig(caminho_completo, dpi=300, bbox_inches='tight')
