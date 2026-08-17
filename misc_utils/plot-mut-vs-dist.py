import argparse
import numpy as np
import matplotlib.pyplot as plt


def read_phy_as_pair_dict(path):
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    n = int(lines[0].split()[0])

    names = []
    rows = []

    for line in lines[1:]:
        parts = line.split()
        name = parts[0]
        values = list(map(float, parts[1:]))

        names.append(name)
        rows.append(values)

    if len(names) != n:
        raise ValueError(f"{path}: expected {n} rows, found {len(names)}")

    matrix = np.array(rows, dtype=float)

    if matrix.shape != (n, n):
        raise ValueError(f"{path}: expected {n}x{n} matrix, found {matrix.shape}")

    pair_values = {}

    for i in range(n):
        for j in range(i + 1, n):
            key = tuple(sorted((names[i], names[j])))
            pair_values[key] = matrix[i, j]

    return names, pair_values


def main():
    parser = argparse.ArgumentParser(
        description="Plot evolutionary distance vs number of mutations from two PHYLIP matrices."
    )

    parser.add_argument(
        "mutation_phy",
        help="PHYLIP matrix containing number of mutations, e.g. mutation.count.phy"
    )

    parser.add_argument(
        "distance_phy",
        help="PHYLIP matrix containing evolutionary distances"
    )

    parser.add_argument(
        "-o", "--output",
        default="mutation_vs_distance.png",
        help="Output figure name"
    )

    parser.add_argument(
        "--label",
        default="Corallel",
        help="Plot title"
    )

    args = parser.parse_args()

    mut_names, mut_pairs = read_phy_as_pair_dict(args.mutation_phy)
    dist_names, dist_pairs = read_phy_as_pair_dict(args.distance_phy)

    common_pairs = sorted(set(mut_pairs) & set(dist_pairs))

    if not common_pairs:
        raise ValueError("No matching sample pairs found between the two PHYLIP files.")

    x = np.array([mut_pairs[pair] for pair in common_pairs], dtype=float)
    y = np.array([dist_pairs[pair] for pair in common_pairs], dtype=float)

    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept

    r = np.corrcoef(x, y)[0, 1]
    r2 = r ** 2

    order = np.argsort(x)

    # use a cleaner built-in matplotlib style
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(7, 5), dpi=150)

    # scatter plot
    ax.scatter(
        x,
        y,
        s=42,
        alpha=0.75,
        edgecolor="black",
        linewidth=0.4,
        label=args.label
    )

    # trend line
    ax.plot(
        x[order],
        y_pred[order],
        color="red",
        linewidth=2.2,
        label=f"Linear fit, R² = {r2:.3f}"
    )

    # labels
    ax.set_xlabel("Number of Mutations", fontsize=8)
    ax.set_ylabel("Evolutionary Distance", fontsize=8)

    # grid and spines
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        frameon=True,
        fontsize=9,
        loc="best"
    )

    # cleaner tick labels
    ax.tick_params(axis="both", labelsize=10)

    fig.tight_layout()
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"Mutation samples: {len(mut_names)}")
    print(f"Distance samples: {len(dist_names)}")
    print(f"Common pairs used: {len(common_pairs)}")
    print(f"Pearson r: {r:.6f}")
    print(f"R^2: {r2:.6f}")
    print(f"Slope: {slope:.12g}")
    print(f"Intercept: {intercept:.12g}")
    print(f"Saved plot: {args.output}")


if __name__ == "__main__":
    main()