#!/bin/bash

# Configuração inicial
OUTDIR="resultados_likwid"                                                             # Diretório para armazenar os resultados
GERA_ENTRADA="./gera_entrada"                                                          # Caminho para o programa gera_entrada
VERSOES=("ajustaPol" "ajustaPolMelhorado")                                                # Versões do programa
N_VALUES=(10 1000)                                                                     # Valores de N
K_VALUES=(64 128 200 256 512 600 800 1024 2000 3000 4096 6000 7000 10000 50000 100000) # Valores de K
K_VALUES_EXTRA=(1000000 10000000 100000000)                                               # Valores extras para N=10

# Metricas
METRICAS=("L3CACHE" "ENERGY" "FLOPS_DP")

CORE=3

# Criar diretório de saída
mkdir -p $OUTDIR

# Função para rodar os testes com LIKWID
run_test() {
    local n=$1
    local k=$2
    local version=$3

    for metrica in "${METRICAS[@]}"; do
        local output_file="$OUTDIR/${version}_N${n}_K${k}_${metrica}.txt"

        echo "Executando: N=$n, K=$k, Versão=$version, Métrica=$metrica"

        # Executar com LIKWID e salvar a saída
        $GERA_ENTRADA $k $n | likwid-perfctr -C $CORE -g $metrica -m -O ./$version >$output_file
    done
}

# Loop principal para os testes
for n in "${N_VALUES[@]}"; do
    for k in "${K_VALUES[@]}"; do
        for version in "${VERSOES[@]}"; do
            run_test $n $k $version
        done
    done

    # Valores extras de K para N=10
    if [[ $n -eq 10 ]]; then
        for k in "${K_VALUES_EXTRA[@]}"; do
            for version in "${VERSOES[@]}"; do
                run_test $n $k $version
            done
        done
    fi
done

echo "Todos os testes foram concluídos. Resultados armazenados no diretório $OUTDIR."
