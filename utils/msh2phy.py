import sys
import numpy as np

def parse(infile):
    distances = {}
    ids = set()

    with open(infile, 'r') as f:
        for line in f:
            reference, query, distance, pvalue, shared_hashes = line.strip().split()
            distance = float(distance)
            
            ids.add(reference)
            ids.add(query)
            
            distances[(reference, query)] = distance
            # distances[(query, reference)] = distance

    return list(ids), distances

def write(outfile, ids, distances):
    n = len(ids)
    
    with open(outfile, 'w') as f:
        f.write(f"{n}\n")
        
        for i in range(n):
            f.write(f"{ids[i]:<10}")
            for j in range(n):
                if i == j:
                    f.write(f"{0:.6f} ")
                else:
                    dist = distances.get((ids[i], ids[j]), 0.0)
                    f.write(f"{dist:.6f} ")
            f.write("\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 msh2phy.py <input_file> <output_file>")
        sys.exit(1)
    
    infile = sys.argv[1]
    outfile = sys.argv[2]

    ids, distances = parse(infile)

    write(outfile, ids, distances)

if __name__ == "__main__":
    main()
