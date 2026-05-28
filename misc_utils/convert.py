import sys

"""
Helper function to convert Mash output file to .phy
"""

def generate_phy_matrix(input_file, output_file):
    """
    Generates a .phy file in matrix format based on the input distances.
    
    Parameters:
    - input_file (str): Path to the input file (tab-delimited: node1, node2, distance).
    - output_file (str): Path to the output .phy file.
    """
    try:
        # Read the input file and store distances in a dictionary
        distances = {}
        nodes = set()

        with open(input_file, 'r') as infile:
            for line in infile:
                # Skip empty lines or lines with insufficient columns
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                
                node1, node2, distance = parts[0], parts[1], float(parts[2])
                nodes.add(node1)
                nodes.add(node2)
                
                # Store distances in a dictionary
                if node1 not in distances:
                    distances[node1] = {}
                if node2 not in distances:
                    distances[node2] = {}
                distances[node1][node2] = distance
                distances[node2][node1] = distance

        # Write the output file
        with open(output_file, 'w') as outfile:
            # Write the number of nodes
            outfile.write(f"{len(nodes)}\n")
            
            # Write the distance matrix
            for node1 in nodes:
                row = [f"{distances.get(node1, {}).get(node2, 0):.10f}" for node2 in nodes]
                formatted_name = f"{node1[:10]:<10}"
                outfile.write(f"{formatted_name} {' '.join(row)}\n")
        
        print(f".phy file successfully written to: {output_file}")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 convert.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        generate_phy_matrix(input_file, output_file)
