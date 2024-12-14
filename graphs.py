import os
import matplotlib.pyplot as plt
import re

# Diretório contendo os arquivos
resultados_dir = "resultados_likwid"

# Função para carregar métricas usando regex
def load_metrics(input_dir, metric_name, metric_key):
    data = {}
    for file in os.listdir(input_dir):
        if metric_name in file:
            # Extrair informações do nome do arquivo
            match_k = re.search(r'_K(\d+)_', file)
            if not match_k:
                continue
            K = int(match_k.group(1))

            # Abrir e ler o conteúdo do arquivo
            with open(os.path.join(input_dir, file), "r") as f:
                content = f.read()

                # Buscar a métrica desejada no conteúdo
                metric_regex = {
                    "L3CACHE": r"L3 miss ratio[:,\s]+([\d.]+)",
                    "ENERGY": r"Energy \[J\][:,\s]+([\d.]+)",
                    "FLOPS_DP": r"DP MFLOP/s[:,\s]+([\d.]+)",
                    "FLOPS_AVX_DP": r"AVX DP MFLOP/s[:,\s]+([\d.]+)",
                }

                match = re.search(metric_regex.get(metric_name, ""), content)
                if match:
                    value = float(match.group(1))
                    label = file.replace(".txt", "").replace("_", " ")
                    if label not in data:
                        data[label] = {"K": [], "value": []}
                    data[label]["K"].append(K)
                    data[label]["value"].append(value)
    return data

# Função para plotar os gráficos
def plot_results(data, metric, ylabel, output_file):
    plt.figure(figsize=(10, 6))

    for label, results in data.items():
        plt.plot(results["K"], results["value"], label=label, marker="o")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Número de pontos (K)")
    plt.ylabel(ylabel)
    plt.title(f"Desempenho de {metric}")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.savefig(output_file)
    plt.close()

# Processar os arquivos e gerar gráficos
if __name__ == "__main__":
    metrics = [
        ("L3CACHE", "Cache miss ratio", "l3cache_miss_ratio.png"),
        ("ENERGY", "Energia (J)", "energia.png"),
        ("FLOPS_DP", "FLOPS DP (MFLOP/s)", "flops_dp.png"),
        ("FLOPS_AVX_DP", "FLOPS AVX DP (MFLOP/s)", "flops_avx_dp.png"),
    ]

    for metric_name, ylabel, output_file in metrics:
        data = load_metrics(resultados_dir, metric_name, ylabel)
        plot_results(data, metric_name, ylabel, output_file)

    print("Gráficos gerados com sucesso!")

