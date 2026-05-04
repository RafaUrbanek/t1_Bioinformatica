import os
import argparse

# Função para formatar o cabeçalho das sequências no arquivo .fasta
def format_fasta(input_file, output_file):
    try:
        with open(input_file, 'r') as input, open(output_file, 'w') as output:
            for line in input:
                # Verifica se a linha é um cabeçalho (começa com >)
                if line.startswith('>'):
                    # Remove espaços em branco extras
                    clean_header = line.replace(" |", "|").strip()
                    output.write(clean_header + "\n")
                else:
                    # Se for a sequência de DNA/Proteína, apenas escreve como está
                    output.write(line)
    
    except FileNotFoundError:
        print("Erro: O arquivo de entrada não foi encontrado.")
        
# Definindo os argumentos de linha de código
parser = argparse.ArgumentParser(description="Script para formatação do arquivo .fasta.")

# Define o argumento 'input'
parser.add_argument("input", help="Caminho para o arquivo .fasta de entrada")
# Define o argumento 'output'
parser.add_argument("output", help="Caminho para o arquivo .fasta de saída")

args = parser.parse_args()

format_fasta(args.input,args.output)