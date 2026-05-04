import argparse
from Bio import SeqIO

# Definição do intervalo fixo de fragmentação
SSEQ = 21563
ESEQ = 25348
# Definição do intervalo deslizante de fragmentação
WINDOW = 600
STEP = 200

# Função para fragmentar os genes
def extract_fragment(input_file, output_file, window_size, step_size):
    fragments = []

    # Percorrendo todas as sequências
    for register in SeqIO.parse(input_file, "fasta"):

        # Tamanho da sequência
        seq_len = len(register)
        
        # O range vai de 0 até o tamanho da sequência descartanto os finais menores que a janela, pulando de STEP em STEP.
        for start_idx in range(0, seq_len - (window_size//2) + 1, step_size):
            end_idx = start_idx + window_size

            fragment = register[start_idx:end_idx]

            # Descarta fragmentos com muitas bases desconhecidas
            if fragment.seq.count("N") / len(fragment) > 0.05:
                
                # Refragmentando para tentar obter maior aproveitamento
                for sub_start in range(0, len(fragment) - (window_size//2) + 1, step_size//2):
                    sub_end = sub_start + window_size//2
                    
                    sub_fragment = fragment[sub_start:sub_end]
                    
                    # Verifica se o sub-fragmento atingiu qualidade suficiente
                    if sub_fragment.seq.count("N") / len(sub_fragment) <= 0.05:
                        
                        # Calcula a posição absoluta na sequência original para o ID
                        abs_start = start_idx + sub_start
                        abs_end = start_idx + sub_end
                        
                        sub_fragment.id = f"{register.id}|pos({abs_start}-{abs_end})"
                        sub_fragment.description = ""
                        fragments.append(sub_fragment)
                
            else:
                # Fragmento original passou direto no controle de qualidade
                fragment.id = f"{register.id}|pos({start_idx}-{end_idx})"
                fragment.description = ""
                fragments.append(fragment)
    
    SeqIO.write(fragments, output_file, "fasta")
    print(f" Sucesso! {len(fragments)} fragmentos gerados usando janelas de {window_size} com saltos de {step_size} salvos em {output_file}.")

# Definindo os argumentos de linha de código
parser = argparse.ArgumentParser(description="Script para fragmentação de sequências FASTA.")

# Define o argumento 'input'
parser.add_argument("input", help="Caminho para o arquivo .fasta de entrada")
# Define o argumento 'output'
parser.add_argument("output", help="Caminho para o arquivo .fasta de saída")

args = parser.parse_args()

extract_fragment(args.input,args.output,WINDOW,STEP)