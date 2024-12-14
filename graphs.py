import os
import matplotlib.pyplot as plt
import re

# Diretório contendo os arquivos
resultados_dir = "resultados_likwid"

# Função para carregar métricas usando regex
def load_metrics(input_dir, metric_name):
    data = {
        "ajustaPol_N10": {"K": [], "value": []},
        "ajustaPol_N1000": {"K": [], "value": []},
        "ajustaPolMelhorado_N10": {"K": [], "value": []},
        "ajustaPolMelhorado_N1000": {"K": [], "value": []},
    }

    for file in os.listdir(input_dir):
        # Identificar a versão e o valor de N a partir do nome do arquivo
        if "ajustaPol" in file:
            if "Melhorado" in file:
                version = "ajustaPolMelhorado"
            else:
                version = "ajustaPol"

            if "_N10_" in file:
                N = "N10"
            elif "_N1000_" in file:
                N = "N1000"
            else:
                continue

            # Extrair o valor de K
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
                    key = f"{version}_{N}"
                    data[key]["K"].append(K)
                    data[key]["value"].append(value)

    return data

# Função para plotar os gráficos
def plot_results(data, metric, ylabel, output_file):
    plt.figure(figsize=(10, 6))

    for label, results in data.items():
        # Ordenar os valores de K e value para formar retas contínuas
        sorted_data = sorted(zip(results["K"], results["value"]))
        sorted_k, sorted_value = zip(*sorted_data)
        plt.plot(sorted_k, sorted_value, label=label, marker="o")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Número de pontos (K)")
    plt.ylabel(ylabel)
    plt.title(f"Desempenho de {metric}")
    plt.legend()  # Adiciona a legenda
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
        data = load_metrics(resultados_dir, metric_name)
        plot_results(data, metric_name, ylabel, output_file)

    print("Gráficos gerados com sucesso!")

