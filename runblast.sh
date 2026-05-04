#!/bin/bash

# --- Paths Gerais ---
SEQ="$1"
DATASET="./dataset"
SCRIPTS="./scripts"
FRAGMENTS="$DATASET/fragments"
SEQUENCES="$DATASET/sequences"
BLASTDB="$DATASET/blastdb"
OUT_DIR="./results"
OUT_FILE="$OUT_DIR/results_$SEQ.csv"

# --- Paths Formatação ---
FORIN="$SEQUENCES/$SEQ.fasta"
FOROUT="$SEQUENCES/formated_$SEQ.fasta"

# --- Paths Fragmentação ---
FRAGIN="$FOROUT"
FRAGOUT="$FRAGMENTS/fragmented_$SEQ.fasta"

# --- Paths Banco de Dados blastn ---
DBIN="$FOROUT"
DBOUT="$BLASTDB/$SEQ/db_$SEQ"

# --- Paths Gráficos ---
GOUT="$OUT_DIR/graphs/$SEQ"

# --- Lógica de limpeza ---
if [[ "$1" == "clean" ]]; then
    echo "Removendo resultados..."
    rm -rf "$OUT_DIR"/*
    exit 0
fi

# --- Lógica de limpeza profunda---
if [[ "$1" == "clean-all" ]]; then
    echo "Removendo foramtações, fragmentações, banco de dados e resultados..."
    rm -rf "$FRAGMENTS"/*
    rm -rf "$BLASTDB"/*
    rm -rf "$OUT_DIR"/*
    rm -rf "$SEQUENCES"/formated*
    exit 0
fi

# --- Preparação ---
echo "Iniciando preparação..."

# --- Formatando o arquivo .fasta ---
if [ -f "$FOROUT" ]; then
    echo "Arquivo .fasta formatado $FOROUT encontrado!"
else    
    echo "Realizando formatação do arquivo $SEQ.fasta ..."
    if python3 "$SCRIPTS"/formatdata.py "$FORIN" "$FOROUT"; then
        echo "Formatação realizada com sucesso!"
    else
        echo "Erro na Formatação!"
        exit 1
    fi
fi

# --- Realizando fragmentação dos dados ---
if [ -f "$FRAGOUT" ]; then
    echo "Fragmentos $FRAGOUT encontrado!"
else    
    echo "Realizando fragmentação das sequências em $SEQ..."
    if python3 "$SCRIPTS"/fragmentation.py "$FRAGIN" "$FRAGOUT"; then
        echo "Fragmentação realizada com sucesso!"
    else
        echo "Erro na fragmentação!"
        exit 1
    fi
fi

# --- Criando banco de dados do blastn ---
if [ -f "$DBOUT.nsq" ]; then
    echo "Banco de dados blastn $DBOUT encontrado!"
else
    echo "Criando banco de dados blastn das sequências em $SEQ..."
    makeblastdb -in "$DBIN" -dbtype nucl -out "$DBOUT"
fi

# --- Execução ---
if [ -f "$OUT_FILE" ]; then
    echo "Arquivo de saída do blastn $OUT_FILE encontrado!"
else
    echo "Rodando blastn..."

    blastn \
        -query "$FRAGOUT" \
        -db "$DBOUT" \
        -out "$OUT_FILE" \
        -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    if [ $? -eq 0 ]; then
        echo "BLAST finalizado com sucesso! Saída: $OUT_FILE"
    else
        echo "Erro no blastn."
        exit 1
    fi
fi

# --- Extraindo resultados ---
echo "Gerando gráficos..."
if python3 "$SCRIPTS"/graphs.py "$OUT_FILE" "$GOUT"; then
    echo "Gráficos gerados com sucesso!"
else
    echo "Erro na geração dos gráficos!"
    exit 1
fi