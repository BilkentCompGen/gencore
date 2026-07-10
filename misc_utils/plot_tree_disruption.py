"""
Plot weighted phylogenetic distance across minimum core-filtering thresholds.

Examples:
    python plot_tree_disruption.py distances.txt -o primates_tree_disruption.png
    python plot_tree_disruption.py distances.txt --csv parsed_metrics.csv
    python plot_tree_disruption.py distances.txt -o figure.pdf --title "My Figure Title"
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


METRIC_KEYS = {
    "RF distance": "rf",
    "Normalized RF distance": "norm_rf",
    "Weighted RF distance": "weighted_rf",
    "Quartet distance": "quartet",
    "Total quartets": "total_quartets",
    "Normalized quartet distance": "norm_quartet",
}


# Colorblind-safe Okabe-Ito-inspired palette.
COLORS = {
    "line": "#222222",
    "stable": "#7A7A7A",
    "minor": "#E69F00",
    "major": "#D55E00",
    "grid": "#D0D0D0",
    "text": "#222222",
}


def set_publication_style() -> None:
    """Set a clean, publication-friendly matplotlib style."""
    mpl.rcParams.update(
        {
            # Font/export settings
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",

            # Axes
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.linewidth": 0.8,
            "axes.edgecolor": COLORS["text"],
            "axes.labelcolor": COLORS["text"],

            # Ticks
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.color": COLORS["text"],
            "ytick.color": COLORS["text"],
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,

            # Legend
            "legend.fontsize": 8,
            "legend.frameon": True,
            "legend.framealpha": 0.95,
            "legend.edgecolor": "#BBBBBB",

            # Figure saving
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
        }
    )


def parse_threshold(
    filename: str,
    threshold_token: str | None = None,
    fallback_last_number: bool = True,
) -> int:
    """
    Parse the minimum core-filtering threshold from a filename.

    Examples parsed by default:
        sample.cc0.nj.newick
        sample.kc5.nj.newick
        sample.kmer10.nj.newick
        sample.core15.nj.newick
        sample.mincore20.nj.newick

    If threshold_token is given, it is tried first.
    For example, threshold_token='cc' parses cc0, cc1, cc2, ...
    """
    patterns = []

    if threshold_token:
        token = re.escape(threshold_token)
        patterns.append(rf"(?:^|[._-]){token}(?:[=_])?(?P<x>\d+)(?:[._-]|$)")

    patterns.extend(
        [
            r"(?:^|[._-])cc(?P<x>\d+)(?:[._-]|$)",
            r"(?:^|[._-])kc(?P<x>\d+)(?:[._-]|$)",
            r"(?:^|[._-])kmer(?P<x>\d+)(?:[._-]|$)",
            r"(?:^|[._-])core(?P<x>\d+)(?:[._-]|$)",
            r"(?:^|[._-])min(?:core|kmer)?(?P<x>\d+)(?:[._-]|$)",
        ]
    )

    for pattern in patterns:
        match = re.search(pattern, filename, flags=re.IGNORECASE)
        if match:
            return int(match.group("x"))

    if fallback_last_number:
        numbers = re.findall(r"\d+", filename)
        if numbers:
            return int(numbers[-1])

    raise ValueError(f"Could not parse threshold from filename: {filename!r}")


def parse_lcp_level(filename: str) -> int | None:
    """
    Parse LCP level from filename, if present.

    Examples:
        sample.lvl6.cc0.nj.newick
        sample.lcp6.cc0.nj.newick
        sample.lcp_6.cc0.nj.newick
    """
    pattern = r"(?:^|[._-])(?:lvl|lcp)(?:[=_])?(?P<x>\d+)(?:[._-]|$)"
    match = re.search(pattern, filename, flags=re.IGNORECASE)

    if not match:
        return None

    return int(match.group("x"))


def parse_numeric_value(value: str) -> int | float | str:
    """Parse numeric values while preserving integer-looking values as ints."""
    try:
        number = float(value)
    except ValueError:
        return value

    return int(number) if number.is_integer() else number


def parse_distance_file(
    path: str | Path,
    threshold_token: str | None = None,
) -> pd.DataFrame:
    """Parse block-style tree-distance report into a tidy dataframe."""
    path = Path(path)
    text = path.read_text().strip()

    rows = []
    current: dict[str, object] | None = None

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        # Filename/header line
        if ":" not in line:
            if current:
                rows.append(current)

            current = {
                "filename": line,
                "threshold": parse_threshold(line, threshold_token=threshold_token),
                "lcp_level": parse_lcp_level(line),
            }
            continue

        if current is None:
            raise ValueError(f"Metric line appeared before filename header: {line!r}")

        key, value = [x.strip() for x in line.split(":", 1)]

        if key in METRIC_KEYS:
            current[METRIC_KEYS[key]] = parse_numeric_value(value)

    if current:
        rows.append(current)

    if not rows:
        raise ValueError(f"No records parsed from {path}")

    df = pd.DataFrame(rows).sort_values("threshold").reset_index(drop=True)

    required = [
        "threshold",
        "rf",
        "norm_rf",
        "weighted_rf",
        "quartet",
        "norm_quartet",
    ]

    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required metrics in parsed file: {missing}")

    # Signed deltas. These can be positive, negative, or zero.
    df["topology_changed"] = df["rf"] > 0
    df["rf_delta_from_previous"] = df["rf"].diff().fillna(0)
    df["weighted_rf_delta_from_previous"] = df["weighted_rf"].diff().fillna(0)
    df["quartet_delta_from_previous"] = df["quartet"].diff().fillna(0)

    return df


def first_threshold_where(df: pd.DataFrame, condition_col: str) -> int | None:
    hits = df.loc[df[condition_col]]

    if hits.empty:
        return None

    return int(hits.iloc[0]["threshold"])


def rf_category(rf: float) -> str:
    """
    Categorize RF disruption severity.

    RF = 0 means topology is preserved.
    Larger RF values indicate stronger topological disruption.
    """
    if rf == 0:
        return "stable"

    if rf < 4:
        return "minor"

    return "major"


def default_title(df: pd.DataFrame) -> str:
    lcp_levels = sorted(x for x in df["lcp_level"].dropna().unique())

    if len(lcp_levels) == 1:
        return (
            "Weighted Phylogenetic Distance Across Minimum Core Filtering "
            f"Thresholds (LCP Level {int(lcp_levels[0])})"
        )

    return "Weighted Phylogenetic Distance Across Minimum Core Filtering Thresholds"


def add_breakpoint_annotations(
    ax: plt.Axes,
    df: pd.DataFrame,
    max_labels: int = 10,
) -> None:
    """Mark thresholds where RF distance changes."""
    breakpoints = df.loc[df["rf_delta_from_previous"] != 0]

    if breakpoints.empty:
        return

    ymin, ymax = ax.get_ylim()
    label_y = ymin + 0.96 * (ymax - ymin)

    for _, row in breakpoints.iterrows():
        threshold = row["threshold"]

        ax.axvline(
            threshold,
            color=COLORS["grid"],
            linestyle=(0, (2, 2)),
            linewidth=0.9,
            zorder=0,
        )

    if len(breakpoints) > max_labels:
        return

    for _, row in breakpoints.iterrows():
        ax.text(
            row["threshold"],
            label_y,
            f"RF {int(row['rf'])}",
            rotation=90,
            ha="right",
            va="top",
            fontsize=7,
            color=COLORS["text"],
            alpha=0.85,
        )


def plot_disruption(
    df: pd.DataFrame,
    output: str | Path,
    title: str | None = None,
    label_breakpoints: bool = True,
) -> None:
    """Create a publication-ready disruption plot."""
    set_publication_style()

    output = Path(output)

    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)

    x = df["threshold"]
    y = df["weighted_rf"]

    plot_df = df.copy()
    plot_df["rf_category"] = plot_df["rf"].apply(rf_category)

    # Main connecting line.
    ax.plot(
        x,
        y,
        color=COLORS["line"],
        linewidth=1.8,
        alpha=0.82,
        zorder=1,
    )

    category_style = {
        "stable": {
            "label": "RF = 0",
            "color": COLORS["stable"],
            "marker": "o",
        },
        "minor": {
            "label": "RF = 2",
            "color": COLORS["minor"],
            "marker": "s",
        },
        "major": {
            "label": "RF = 4",
            "color": COLORS["major"],
            "marker": "^",
        },
    }

    # Draw points by category so the legend is clean and dynamic.
    for category, style in category_style.items():
        subset = plot_df.loc[plot_df["rf_category"] == category]

        if subset.empty:
            continue

        ax.scatter(
            subset["threshold"],
            subset["weighted_rf"],
            s=58,
            marker=style["marker"],
            facecolor=style["color"],
            edgecolor="#222222",
            linewidth=0.6,
            zorder=3,
            label=style["label"],
        )

    # Y-axis padding that also works when all values are equal.
    ymin = float(y.min())
    ymax = float(y.max())
    yrange = ymax - ymin

    if yrange == 0:
        pad = max(abs(ymax) * 0.15, 0.01)
    else:
        pad = yrange * 0.14

    ax.set_ylim(ymin - pad, ymax + pad)

    if label_breakpoints:
        add_breakpoint_annotations(ax, plot_df)
    
    # Highlight default threshold.
    default_threshold = 32

    ax.axvline(
        default_threshold,
        color="red",
        linestyle="--",
        linewidth=1.2,
        alpha=0.9,
        zorder=2,
    )

    ax.text(
        default_threshold,
        ax.get_ylim()[1],
        "Default",
        rotation=90,
        ha="right",
        va="top",
        fontsize=8,
        color="red",
        alpha=0.5,
    )

    # Labels
    ax.set_xlabel("Minimum LCP core-filtering threshold")
    ax.set_ylabel("Weighted RF distance")

    # if title is None:
    #     title = default_title(df)

    ax.set_title(title, pad=10)

    # Ticks
    if len(df) <= 16:
        ax.set_xticks(df["threshold"])
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))

    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))

    # Grid: keep it subtle and mostly horizontal.
    ax.grid(
        True,
        axis="y",
        color=COLORS["grid"],
        linewidth=0.7,
        linestyle="-",
        alpha=0.75,
    )
    ax.grid(
        True,
        axis="x",
        color=COLORS["grid"],
        linewidth=0.5,
        linestyle="-",
        alpha=0.28,
    )

    # Cleaner spines.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend includes line meaning plus observed RF categories.
    handles, labels = ax.get_legend_handles_labels()

    line_handle = Line2D(
        [0],
        [0],
        color=COLORS["line"],
        linewidth=1.8,
        alpha=0.82,
        label="Weighted RF distance",
    )

    ax.legend(
        handles=[line_handle, *handles],
        loc="best",
        borderpad=0.6,
        handlelength=1.8,
    )

    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote figure: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot weighted tree disruption across filtering thresholds."
    )

    parser.add_argument(
        "input",
        help="Block-style distance report.",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="tree_disruption.png",
        help="Output figure path. Use .pdf or .svg for vector output.",
    )

    parser.add_argument(
        "--csv",
        default=None,
        help="Optional output CSV of parsed metrics.",
    )

    parser.add_argument(
        "--title",
        default=None,
        help="Optional figure title.",
    )

    parser.add_argument(
        "--threshold-token",
        default=None,
        help=(
            "Optional filename token used to parse thresholds. "
            "For example, '--threshold-token cc' parses cc0, cc1, cc2, ..."
        ),
    )

    parser.add_argument(
        "--no-breakpoint-labels",
        action="store_true",
        help="Do not label RF-change breakpoints on the plot.",
    )

    args = parser.parse_args()

    df = parse_distance_file(args.input, threshold_token=args.threshold_token)

    first_topology_change = first_threshold_where(df, "topology_changed")

    print(f"Parsed {len(df)} thresholds: {df['threshold'].min()}..{df['threshold'].max()}")

    lcp_levels = sorted(x for x in df["lcp_level"].dropna().unique())
    if len(lcp_levels) == 1:
        print(f"LCP level: {int(lcp_levels[0])}")
    elif len(lcp_levels) > 1:
        print(f"LCP levels detected: {', '.join(map(lambda x: str(int(x)), lcp_levels))}")

    if first_topology_change is None:
        print("No topology disruption detected: RF distance stays 0.")
    else:
        print(f"First topology disruption: threshold={first_topology_change}")

    breakpoint_rows = df.loc[
        df["rf_delta_from_previous"] != 0,
        [
            "threshold",
            "rf",
            "rf_delta_from_previous",
            "norm_rf",
            "quartet",
            "quartet_delta_from_previous",
            "norm_quartet",
            "weighted_rf",
            "weighted_rf_delta_from_previous",
        ],
    ]

    if not breakpoint_rows.empty:
        print("\nTopology breakpoints:")
        print(breakpoint_rows.to_string(index=False))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"Wrote parsed metrics: {args.csv}")

    plot_disruption(
        df,
        args.output,
        title=args.title,
        label_breakpoints=not args.no_breakpoint_labels,
    )


if __name__ == "__main__":
    main()