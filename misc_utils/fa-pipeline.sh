# HG001_GRCh38_1_22_v4.2.1_all.vcf.gz: 
#   /giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/SupplementaryFiles/HG001_GRCh38_1_22_v4.2.1_all.vcf.gz
# Cosmic_NonCodingVariants_v100_GRCh38.vcf.gz: 
#   https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/noncodingvariantsvcf (login required)
#
# Required packages:
#   biopython (pip3)
#   ete3 (pip3)
#   six (pip3)
#   PyQt5 (pip3)
# Required scripts and commands:
#   gencore
#   mash
#   phylowizard.py
#   psite.py
#   convert.py (.dist to .phy)
# Required files:
#   hg38.fa
#   HG001_GRCh38_1_22_v4.2.1_all.vcf.gz
#   cfg_template_female.yaml
#   fa-primates
#       ├── hg002v1.1.fasta
#       ├── mGorGor1.dip.cur.20231122.fasta
#       ├── mPanPan1.dip.cur.20231122.fasta
#       ├── mPanTro3.dip.cur.20231122.fasta
#       ├── mPonAbe1.dip.cur.20231205.fasta
#       ├── mPonPyg2.dip.cur.20231122.fasta
#       ├── mSymSyn1.analysis-dip.20240514.fasta
#       ├── input.txt
#       └── shortnames.txt
#   
# If the chromosome names in FASTA file is given in chrN format, convert it into numberic style.
#   sed 's/^>chr/>/' human_38.fasta > hg38.fa


mkdir fa-tumor-human
cd fa-tumor-human

# Simulate coalescent tree of 20 tumor cells
ms 20 1 -T -G 1 | tail -n1 > ms_tree.txt

# Filter the raw variants to get the phased SNPs as the germline variants of the sample to simulate.
zcat ../HG001_GRCh38_1_22_v4.2.1_all.vcf.gz | \
    awk '/^#/ || ($NF~/^[01]\|[01]/ && length($4)==1 && length($5)==1)' | \
    gzip -c > HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz

mv HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz hg001.vcf.gz

source ~/scripts/activate_conda

# Simulate the (normal) germline genomes of a female individual, by integrating germline SNPs into the human reference genome. (~3 min.)
/bin/time -v python3 ../psite.py vcf2fa \
    -r ../hg38.fa \
    -v hg001.vcf.gz \
    --autosomes 1..22 \
    --sex_chr X,X \
    -o normal_fa > fa-tumor-pipeline.out 2>&1

# Simulate somatic variants with tumor chain (~30 sec.)
/bin/time -v python3 ../psite.py phylovar \
    -t ms_tree.txt \
    --config ../cfg_template_female.yaml \
    --purity 0.8 \
    --sex_chr X,X \
    --trunk_length 2.0 \
    --prune 0.05 \
    --chain tumor_chain \
    --map map >> fa-tumor-pipeline.out 2>&1

# Build the genomes of tumor cells (~40 min.)
/bin/time -v python3 ../psite.py chain2fa \
    -c tumor_chain \
    -n normal_fa/normal.parental_0.fa,normal_fa/normal.parental_1.fa \
    -o tumor_fa \
    --cores 16 >> fa-tumor-pipeline.out 2>&1

# Merge diploid files
nodes=$(ls ./tumor_fa | grep -oP 'node\d+' | sort | uniq); \
for node in $nodes; do \
	parental_0="./tumor_fa/${node}.parental_0.fa"; \
	parental_1="./tumor_fa/${node}.parental_1.fa"; \
    cat "$parental_0" "$parental_1" > ./tumor_fa/${node}.fa; \
    echo "Merged $parental_0 and $parental_1 into ./tumor_fa/${node}.fa"; \
    rm "$parental_0" "$parental_1"; \
done

dir=$(pwd)
tot=0; snv=0; amp=0; del=0; oth=0; ind=0;

echo "." > "$dir/summary"
for chain in "$dir/tumor_chain"/*.genome.chain; do \
    file=$(basename "$chain" .genome.chain); \
    tot_loc=$(grep -c -vE '(REF|>)' "$chain"); \
    snv_loc=$(grep -c 'SNV' "$chain"); \
    amp_loc=$(grep -c 'AMP' "$chain"); \
    del_loc=$(grep -c 'DEL' "$chain"); \
    oth_loc=$((tot_loc-snv_loc-amp_loc-del_loc)); \
    echo "$file:" >> "$dir/summary"; \
    echo -e "\tALL: ${tot_loc}, SNV: ${snv_loc}, AMP: ${amp_loc}, DEL: ${del_loc}, OTH: ${oth_loc}" >> "$dir/summary"; \
    tot=$((tot + tot_loc)); \
    snv=$((snv + snv_loc)); \
    amp=$((amp + amp_loc)); \
    del=$((del + del_loc)); \
	oth=$((oth + oth_loc)); \
	ind=$((ind + 1)); \
done

echo "Total:" >> "$dir/summary"
echo -e "\tALL: ${tot}, SNV: ${snv}, AMP: ${amp}, DEL: ${del}, OTH: ${oth}" >> "$dir/summary"

if [ "$ind" -gt 0 ]; then \
    avg_tot=$((tot / ind)); \
    avg_snv=$((snv / ind)); \
    avg_amp=$((amp / ind)); \
    avg_del=$((del / ind)); \
    avg_oth=$((oth / ind)); \
    echo "Average per file:" >> "$dir/summary"; \
    echo -e "\tALL: ${avg_tot}, SNV: ${avg_snv}, AMP: ${avg_amp}, DEL: ${avg_del}, OTH: ${avg_oth}" >> "$dir/summary"; \
fi

echo "Number of genomes in ./fa-tumor-human/tumor_fa: $ind"
echo -e "\tALL: ${tot}, SNV: ${snv}, AMP: ${amp}, DEL: ${del}, OTH: ${oth}"

# Create the reference newick
map_file="map/tumor.tipnode.map"
newick_file="ms_tree.txt"

## Rename nodes (to tumor#)
declare -A cell_to_node

while read -r tip_node cell_count cell; do
    if [[ "$tip_node" == \#* ]]; then
        continue  # skip header
    fi
    cell_to_node["$cell"]="$tip_node"
done < "$map_file"

newick=$(cat "$newick_file")

## Find numbers (cells) to replace - stricter match, remove duplicates, sort descending
numbers=$(echo "$newick" | grep -oP '(?<=\(|,)[0-9]+(?=:)')
numbers=$(echo "$numbers" | sort -nru)

## Replace each cell ID with tip_node. Find all unique numbers from Newick
for num in $numbers; do
    node=${cell_to_node[$num]}
    if [ -n "$node" ]; then
        newick=$(echo "$newick" | sed -E "s/([,(])$num:([0-9])/\1$node:\2/g")
    fi
done

## map nodeXX to tumorY. Find all unique node names
node_names=$(echo "$newick" | grep -oE 'node[0-9]+' | sort -uV)

declare -A node_to_tumor
tumor_id=1
for node in $node_names; do
    node_to_tumor["$node"]="tumor$tumor_id"
    tumor_id=$((tumor_id + 1))
done

for node in "${!node_to_tumor[@]}"; do
    tumor=${node_to_tumor[$node]}
    newick=$(echo "$newick" | sed -E "s/\b$node\b/$tumor/g")
done

echo "$newick" > ms_tree.newick

# Create input.txt and shortnames.txt files
cd tumor_fa

## Get fasta file names
find . -maxdepth 1 -name "node*.fa" | sort -V > input.txt
## Trim the names
cat input.txt | sed 's|^\./||; s|\.fa$||' > input-trimmed.txt
## Create short names
for i in $(seq 1 $(wc -l < input-trimmed.txt)); do \
    echo "tumor$i"; \
done > shortnames.txt

source ~/scripts/activate_conda

rm -f gencore-tumor-fa-out.txt mash-tumor-fa-out.txt

for l in 4 5 6 7 8 9; do \
    /bin/time -v ../../gencore fa \
    -i input.txt \
    -s shortnames.txt \
    -t 4 \
    -l "$l" \
    -p tumor \
    -v >> gencore-tumor-fa-out.txt 2>&1; \
    /bin/time -v python3 ../../phylowizard.py tumor.set.jc.lvl${l}.phy --normalize; \
done

for s in 1000 5000 50000 500000 5000000; do \
    /bin/time -v mash sketch \
        -o tumor.mash.${s}.msh *.fa \
        -s ${s} \
        -p 4 >> mash-tumor-fa-out.txt 2>&1; \
    /bin/time -v mash dist \
         tumor.mash.${s}.msh tumor.mash.${s}.msh > tumor.mash.${s}.dist; \
    python3 ../../convert.py tumor.mash.${s}.dist tumor.mash.${s}.phy; \
    /bin/time -v python3 ../../phylowizard.py tumor.mash.${s}.phy --normalize; \
done

for s in 1000 5000 50000 500000 5000000; do 
    newick=$(cat "tumor.mash.$s.nj.newick" | sed -E "s/.fa//g")
    node_names=$(echo "$newick" | grep -oE 'node[0-9]+' | sort -uV)
    declare -A mash_node_to_tumor
    tumor_id=1
    for node in $node_names; do
        mash_node_to_tumor["$node"]="tumor$tumor_id"
        tumor_id=$((tumor_id + 1))
    done
    for node in "${!mash_node_to_tumor[@]}"; do
        tumor=${mash_node_to_tumor[$node]}
        newick=$(echo "$newick" | sed -E "s/\b$node\b/$tumor/g")
    done
    echo "$newick" > "tumor.mash.$s.nj.renamed.newick"
done


cd ../../

cd fa-primates

# input.txt:
#   find . -maxdepth 1 -name "*.fasta" | sort -V | sed 's|^\./||' > input.txt
# shortnames.txt:
#   Human
#   Gorilla
#   Bonobo
#   Chimpanzee
#   Sum_Orang
#   Bor_Orang
#   Siamang

rm -f gencore-primates-fa-out.txt mash-primates-fa-out.txt

for l in 4 5 6 7; do \
    /bin/time -v ../gencore fa \
    -i input.txt \
    -s shortnames.txt \
    -t 7 \
    -l "$l" \
    -p primates \
    -v >> gencore-primates-fa-out.txt 2>&1; \
    /bin/time -v python3 ../phylowizard.py primates.set.jc.lvl${l}.phy; \
done

/bin/time -v mash sketch \
    -o primates.mash.1000.msh *.fasta \
    -s 1000 \
    -p 7 >> mash-primates-fa-out.txt 2>&1; \
/bin/time -v mash dist primates.mash.1000.msh primates.mash.1000.msh > primates.mash.1000.dist
python3 ../convert.py primates.mash.1000.dist primates.mash.1000.phy
/bin/time -v python3 ../phylowizard.py primates.mash.1000.phy

conda deactivate