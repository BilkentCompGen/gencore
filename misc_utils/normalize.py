import re

# Step 1: Read ms_tree.txt and extract leaf numbers
tree_file = "ms_tree.txt"
with open(tree_file, "r") as f:
    tree = f.read().strip()

# Find numbers before colon
numbers = re.findall(r"(\d+):", tree)
unique_numbers = []
for n in numbers:
    if n not in unique_numbers:
        unique_numbers.append(n)

# Create mapping {original_number : id}
num_to_id = {int(n): i + 1 for i, n in enumerate(unique_numbers)}

# Step 2: Replace numbers with ids in tree
def replace_number_with_id(match):
    num = int(match.group(1))
    return f"tumor{num_to_id[num]}:"

new_tree = re.sub(r"(\d+):", replace_number_with_id, tree)

with open("ms_tree.newick", "w") as f:
    f.write(new_tree + "\n")

# Step 3: Read tumor.tipnode.map
map_file = "map/tumor.tipnode.map"
tip_map = {}

with open(map_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tip_node, _, cell_str = parts
            try:
                cell = int(cell_str)
                tip_map[cell] = tip_node
            except ValueError:
                continue

# Step 4: Write tip ids according to tumor.tipnode.map order
with open("tumor_fa/shortnames.txt", "w") as f:
    for cell, tip_node in tip_map.items():
        if cell in num_to_id:
            f.write(f"tumor{num_to_id[cell]}\n")
        else:
            sys.exit(1)