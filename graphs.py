import os
import matplotlib.pyplot as plt
import re

def load_metrics(input_dir, metric_name):
    data = {}
    versions = ["ajustaPol", "ajustaPolMelhorado"]
    N_values = [10, 1000]

    for version in versions:
        for N in N_values:
            key = f"{version}_N{N}"
            data[key] = {"K": [], "value": []}

            for file in os.listdir(input_dir):
                if file.startswith(f"{version}_N{N}") and metric_name in file:
                    K = int(re.search(r'_K(\d+)_', file).group(1))
                    with open(os.path.join(input_dir, file), "r") as f:
                        content = f.read()
                        if metric_name == "L3CACHE":
                            match = re.search(r"L3 cache miss ratio:\s+([\d.]+)", content)
                        elif metric_name == "ENERGY":
                            match = re.search(r"Energy \[J\]:\s+([\d.]+)", content)
                        elif metric_name == "FLOPS_DP":
                            match = re.search(r"DP MFLOP/s:\s+([\d.]+)", content)
                        elif metric_name == "FLOPS_AVX_DP":
                            match = re.search(r"AVX DP MFLOP/s:\s+([\d.]+)", content)
                        else:
                            match = None

                        if match:
                            value = float(match.group(1))
                            data[key]["K"].append(K)
                            data[key]["value"].append(value)

    return data

def plot_results(data, metric, ylabel, output_file):
    plt.figure(figsize=(10, 6))

    for key, results in data.items():
        plt.plot(
            results["K"], results["value"], label=key, marker="o"
        )

    plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Número de pontos (K)")
    plt.ylabel(ylabel)
    plt.title(f"Desempenho de {metric}")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_dir = "resultados_likwid"

    # Carregar dados e gerar gráficos
    metrics = [
        ("L3CACHE", "Cache miss ratio", "l3cache_miss_ratio.png"),
        ("ENERGY", "Energia (J)", "energia.png"),
        ("FLOPS_DP", "FLOPS DP (MFLOP/s)", "flops_dp.png"),
        ("FLOPS_AVX_DP", "FLOPS AVX DP (MFLOP/s)", "flops_avx_dp.png")
    ]

    for metric_name, ylabel, output_file in metrics:
        data = load_metrics(input_dir, metric_name)
        plot_results(data, metric_name, ylabel, output_file)

    print("Gráficos gerados com sucesso!")
