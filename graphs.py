import os
import matplotlib.pyplot as plt
import numpy as np

# Diretório contendo os arquivos
resultados_dir = "resultados_likwid"

# Métricas a serem extraídas
metrics = [
    ("L3CACHE", "L3 miss ratio", "l3cache_miss_ratio.png", "L3 miss ratio"),
    ("ENERGY", "Energy [J]", "energia.png", "Energy [J]"),
    ("FLOPS_DP", "DP MFLOP/s", "flops_dp.png", "DP MFLOP/s"),
]

# Função para extrair valores de um arquivo
def extract_metric(file_path, metric_name, metric_key):
    k_values = []
    metric_values = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_k = None
        
        for line in lines:
            # Procurar o valor de K no início do arquivo
            if line.strip().isdigit():
                current_k = int(line.strip())
            
            # Localizar a tabela de métricas
            if metric_name in line and "Group 1 Metric" in line:
                for metric_line in lines[lines.index(line):]:
                    if metric_key in metric_line:
                        # Tentar extrair o valor da métrica
                        try:
                            value = float(metric_line.split(',')[-1].strip())
                            if current_k is not None:
                                k_values.append(current_k)
                                metric_values.append(value)
                        except ValueError:
                            print(f"Valor inválido no arquivo {file_path}, linha: {metric_line}")
                        break
    return k_values, metric_values

# Processar os arquivos e gerar gráficos
for metric, ylabel, output_file, metric_key in metrics:
    plt.figure(figsize=(10, 6))
    
    for root, _, files in os.walk(resultados_dir):
        for file in files:
            if metric in file:
                # Identificar a combinação (N1/N2 e V1/V2) pelo nome do arquivo
                label = file.replace(".txt", "").replace("_", " ")
                file_path = os.path.join(root, file)
                
                # Extrair os dados
                k_values, metric_values = extract_metric(file_path, metric, metric_key)
                
                # Plotar os dados
                if k_values and metric_values:
                    plt.plot(k_values, metric_values, label=label)
    
    # Configurar o gráfico
    plt.xscale('log')
    plt.xlabel('Número de pontos (K)')
    plt.ylabel(ylabel)
    plt.title(f'Gráfico para métrica: {ylabel}')
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    
    # Salvar o gráfico
    plt.savefig(output_file)
    plt.close()

