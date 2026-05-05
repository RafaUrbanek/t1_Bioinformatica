#  Trabalho Semestral de Bioinformática

Este projeto automatiza a análise de sequências de nucleotídeos do vírus SARS-CoV-2 (COVID-19) coletadas no Brasil, utilizando scripts Python e ferramentas de alinhamento como o BLAST.

---

##  Requisitos Técnicos

Para executar este projeto, você precisará do **Python 3.11** e da ferramenta **blastn** instalada no seu sistema (via BLAST+ do NCBI).

### Bibliotecas Python Necessárias:
O projeto utiliza bibliotecas para manipulação de dados biológicos e visualização gráfica:
*   **Bioinformática:** `biopython`
*   **Processamento de Dados:** `numpy`, `pandas`
*   **Visualização:** `matplotlib`, `seaborn`

##  Modo de Uso

### 1. Preparação do Banco de Dados
Para garantir a compatibilidade com os scripts, siga rigorosamente as etapas de coleta no **[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049)**:

1. **Filtros Aplicados:**
   - **Região Geográfica:** `Brazil`
   - **Tamanho do Genoma:** Mínimo `29000`
   - **Data de Coleta:** `05/04/2020` a `05/04/2026`
   - **Tipo de Sequência:** `Nucleotídeos`
2. **Seleção de Amostras:** Escolha aleatoriamente **10 valores** (ou o valor `x` desejado) para cada ano de coleta.
3. **Exportação de Colunas:** Ao baixar, selecione exatamente estas três colunas nesta ordem:
   - `Accession` | `Collection_Date` | `Pangolin`
4. **Download:** Salve o arquivo no diretório `./dataset/sequences` do projeto.

### 2. Execução
Antes de rodar o processo, certifique-se de que possui as bibliotecas Python necessárias instaladas em seu ambiente virtual.

Execute o script principal passando o nome do arquivo baixado como argumento:

```bash
# Dar permissão de execução (se necessário)
chmod +x runblast.sh

# Rodar o script utilizando o nome do arquivo sem extensão
./runblast.sh <nome_do_arquivo_baixado>
```

### 3. Como Ver os Resultados
Após a execução, o script processará os dados e gerará saídas organizadas. Veja como encontrá-las e o que elas significam:

📁 Pasta ``/results``
Todos os arquivos gerados são salvos automaticamente nesta pasta:

Arquivos CSV (results_*.csv): Contém a tabela bruta com os resultados do alinhamento blastn, incluindo o E-value, identidade de sequência e as linhagens Pangolin associadas.

Gráficos de Distribuição(``./graphs/<nome da base de dados>``): Gráficos gerados via seaborn e matplotlib que mostram a evolução das variantes ao longo do período selecionado (2020-2026).
