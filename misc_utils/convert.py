import os
import argparse

def clean_name(name):
    name = os.path.basename(name)

    for ext in [".fasta", ".fa", ".fna", ".fas"]:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break

    return name

def generate_phy_matrix(input_file, output_file, tool):
    """
    Convert Mash, Sourmash, or Dashing2 pairwise output to PHYLIP .phy matrix.

    Usage:
        python3 convert.py <input_file> <output_file> --mash
        python3 convert.py <input_file> <output_file> --sourmash
        python3 convert.py <input_file> <output_file> --dashing2
    """

    distances = {}
    nodes = set()

    tool = tool.lower()
    if tool not in {"mash", "sourmash", "dashing2"}:
        raise ValueError("tool must be either 'mash', 'sourmash', or 'dashing2'")

    if tool == "dashing2":
        names = []
        rows = []

        with open(input_file, "r") as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("#Sources"):
                    names = [clean_name(x) for x in line.split()[1:]]
                    continue

                if line.startswith("#"):
                    continue

                parts = line.split()
                rows.append((clean_name(parts[0]), parts[1:]))

        if not names:
            raise ValueError("No #Sources line found in Dashing2 output")

        for name in names:
            nodes.add(name)

        for i, node1 in enumerate(names):
            distances.setdefault(node1, {})
            distances[node1][node1] = 0.0

        for i, (node1, values) in enumerate(rows):
            for j, value in enumerate(values):
                if value == "-":
                    continue

                # Dashing2 gives similarity, convert to distance
                distance = 1.0 - float(value)

                node2 = names[j]

                distances.setdefault(node1, {})
                distances.setdefault(node2, {})

                distances[node1][node2] = distance
                distances[node2][node1] = distance

    else:
        with open(input_file, "r") as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()

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
        description="Convert Mash, Sourmash, or Dashing2 pairwise output to PHYLIP .phy matrix."
    )

    parser.add_argument("input_file", help="Input distance file")
    parser.add_argument("output_file", help="Output .phy file")

    tool_group = parser.add_mutually_exclusive_group(required=True)
    tool_group.add_argument("--mash", action="store_true", help="Input is Mash output")
    tool_group.add_argument("--sourmash", action="store_true", help="Input is Sourmash output")
    tool_group.add_argument("--dashing2", action="store_true", help="Input is Dashing2 output")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tool = "mash"
    if args.sourmash:
        tool = "sourmash"
    elif args.dashing2:
        tool = "dashing2"

    try:
        generate_phy_matrix(args.input_file, args.output_file, tool)
    except Exception as e:
        print(f"Error: {e}")
        raise