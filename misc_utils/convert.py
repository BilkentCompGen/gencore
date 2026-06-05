import sys
import argparse


def generate_phy_matrix(input_file, output_file, tool):
    """
    Convert Mash or Sourmash pairwise distance output to PHYLIP .phy matrix.

    Usage:
        python3 convert.py <input_file> <output_file> <mash|sourmash>
    """

    distances = {}
    nodes = set()

    tool = tool.lower()
    if tool not in {"mash", "sourmash"}:
        raise ValueError("tool must be either 'mash' or 'sourmash'")

    with open(input_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            parts = line.split()

            # Mash output is usually:
            # query reference distance p-value matching-hashes
            if tool == "mash":
                if len(parts) < 3:
                    continue
                node1, node2 = parts[0], parts[1]
                distance = float(parts[2])

            # Sourmash output can vary, but common pairwise CSV/TSV style is:
            # name1 name2 distance / similarity ...
            elif tool == "sourmash":
                if len(parts) < 3:
                    continue
                node1, node2 = parts[0], parts[1]
                distance = float(parts[2])

            nodes.add(node1)
            nodes.add(node2)

            distances.setdefault(node1, {})
            distances.setdefault(node2, {})

            distances[node1][node2] = distance
            distances[node2][node1] = distance

    nodes = sorted(nodes)

    with open(output_file, "w") as outfile:
        outfile.write(f"{len(nodes)}\n")

        for node1 in nodes:
            row = [
                f"{distances.get(node1, {}).get(node2, 0.0):.10f}"
                for node2 in nodes
            ]
            formatted_name = f"{node1[:10]:<10}"
            outfile.write(f"{formatted_name} {' '.join(row)}\n")

    print(f".phy file successfully written to: {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert Mash or Sourmash pairwise distance output to PHYLIP .phy matrix."
    )

    parser.add_argument("input_file", help="Input distance file")
    parser.add_argument("output_file", help="Output .phy file")

    tool_group = parser.add_mutually_exclusive_group(required=True)
    tool_group.add_argument("--mash", action="store_true", help="Input is Mash output")
    tool_group.add_argument("--sourmash", action="store_true", help="Input is Sourmash output")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tool = "sourmash" if args.sourmash else "mash"

    try:
        generate_phy_matrix(args.input_file, args.output_file, tool)
    except Exception as e:
        print(f"Error: {e}")
        raise