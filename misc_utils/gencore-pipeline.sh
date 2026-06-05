# Required files:
#   fa-primates
#       ├── hg002v1.1.fasta (+ .fai)
#       ├── mGorGor1.dip.cur.20231122.fasta (+ .fai)
#       ├── mPanPan1.dip.cur.20231122.fasta (+ .fai)
#       ├── mPanTro3.dip.cur.20231122.fasta (+ .fai)
#       ├── mPonAbe1.dip.cur.20231205.fasta (+ .fai)
#       ├── mPonPyg2.dip.cur.20231122.fasta (+ .fai)
#       ├── mSymSyn1.analysis-dip.20240514.fasta (+ .fai)
#       ├── mSymSyn1.analysis-dip.20240514.fasta.fai
#       ├── filenames.txt
#       └── shortnames.txt
#
#   The content of the input can be filled via following command:
#       filenames.txt:
#           find . -maxdepth 1 -name "*.fasta" | sort -V | sed 's|^\./||' > filenames.txt
#   The content of the shortname is as follows:
#       shortnames.txt:
#           Human
#           Gorilla
#           Bonobo
#           Chimpanzee
#           Sum Orang
#           Bor Orang
#           Siamang
#

RUN_GENCORE_TUMOR="true"
RUN_MASH_TUMOR="true"
RUN_SOURMASH_TUMOR="true"
RUN_ANI="true"
RUN_GENCORE_PRIM="true"
RUN_MASH_PRIM="true"
RUN_SOURMASH_PRIM="true"
RUN_GENCORE_FQ="true"

GENCORE="gencore"
PHYLOWIZARD="phylowizard.py"
PLOT="misc_utils/plot-mut-vs-dist.py"
MASH="mash"
CONVERT="misc_utils/convert.py"
SOURMASH="sourmash"

PRIMATES_FQ_DIR=.

WORK_DIR=.

MUTATION_COUNT_PHY="mutation.count.phy"

CONFIG_FILE="gencore-config.sh"

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
fi

cd $WORK_DIR

# GenCore Tumor
if [ "$RUN_GENCORE_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f gencore-tumor-fa-out.txt

    for l in 4 5 6 7 8 9; do

        /bin/time -v ${GENCORE} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 4 \
            -l "$l" \
            -p tumor \
            -v >> gencore-tumor-fa-out.txt 2>&1

        python3 ${PHYLOWIZARD} tumor.set.evol.lvl${l}.phy --normalize >> gencore-tumor-fa-out.txt 2>&1

    done

    cd ..
fi

# Mash Tumor
if [ "$RUN_MASH_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f mash-tumor-fa-out.txt

    for s in 1000 5000 50000 500000 5000000; do

        /bin/time -v ${MASH} sketch \
            -o tumor.mash.${s}.msh *.fa \
            -s ${s} \
            -p 4 >> mash-tumor-fa-out.txt 2>&1

        /bin/time -v ${MASH} dist \
            tumor.mash.${s}.msh tumor.mash.${s}.msh > tumor.mash.${s}.dist

        python3 ${CONVERT} --mash tumor.mash.${s}.dist tumor.mash.${s}.phy

        if [ -f "tumor.mash.${s}.phy" ]; then
            python3 ${PHYLOWIZARD} tumor.mash.${s}.phy --normalize >> mash-tumor-fa-out.txt 2>&1
        else
            echo "File tumor.mash.${s}.phy not found!"
        fi
    done

    for s in 1000 5000 50000 500000 5000000; do
        if [ -f "tumor.mash.$s.nj.newick" ]; then

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
        fi
    done

    cd ..
fi

# Sourmash tumor
if [ "$RUN_SOURMASH_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f sourmash-tumor-fa-out.txt

    /bin/time -v ${SOURMASH} compute -k 21 *.fa > sourmash-tumor-fa-out.txt 2>&1
    /bin/time -v ${SOURMASH} compare -p 8 *.sig -o sourmash.tumor >> sourmash-tumor-fa-out.txt 2>&1

    cd ..
fi

# GenCore-ANI tumor
if [ "$RUN_ANI" = "true" ]; then

    mkdir -p mut-vs-ani-plots
    cd mut-vs-ani-plots

    if [ -f "${MUTATION_COUNT_PHY}" ]; then
        for l in 4 5 6 7 8 9; do
            python3 ${PLOT} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.jaccard.lvl${l}.phy --label "GenCore (${l})" -o count-gencore-j-${l}.pdf
            python3 ${PLOT} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.evol.lvl${l}.phy --label "GenCore (${l})" -o count-gencore-p-${l}.pdf
        done
    fi

    cd ..
fi

# GenCore primates
if [ "$RUN_GENCORE_PRIM" = "true" ]; then

    cd fa-primates

    rm -f gencore-primates-fa-out.txt

    # GenCore
    for l in 4 5 6 7; do
        /bin/time -v ${GENCORE} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 7 \
            -l "$l" \
            -p primates \
            -v >> gencore-primates-fa-out.txt 2>&1;

        python3 ${PHYLOWIZARD} primates.set.evol.lvl${l}.phy;
    done

    cd ..
fi

# Mash primates
if [ "$RUN_MASH_PRIM" = "true" ]; then

    cd fa-primates

    rm -f mash-primates-fa-out.txt

    /bin/time -v ${MASH} sketch \
        -o primates.mash.1000.msh *.fasta \
        -s 1000 \
        -p 7 >> mash-primates-fa-out.txt 2>&1

    /bin/time -v mash dist primates.mash.1000.msh primates.mash.1000.msh > primates.mash.1000.dist
    python3 ${CONVERT} --mash primates.mash.1000.dist primates.mash.1000.phy
    python3 ${PHYLOWIZARD} primates.mash.1000.phy

    cd ..
fi

# Mash primates
if [ "$RUN_SOURMASH_PRIM" = "true" ]; then

    cd fa-primates

    rm -f sourmash-primates-fa-out.txt

    /bin/time -v ${SOURMASH} compute -k 21 *.fasta > sourmash-primates-fa-out.txt 2>&1
    /bin/time -v ${SOURMASH} compare -p 8 *.sig -o sourmash.primates >> sourmash-primates-fa-out.txt 2>&1

    cd ..
fi

# Gencore primates fq
if [ "$RUN_GENCORE_FQ" = "true" ]; then

    mkdir -p fq-primates

    cd $PRIMATES_FQ_DIR

    rm -f gencore-primates-fq-out.txt

    for l in 4 5 6; do

        /bin/time -v ${GENCORE} fq \
            -i filenames.txt \
            -s shortnames.txt \
            -l "$l" \
            -t 32 \
            -r 6 \
            -p primates.fq \
            --min-cc 32 \
            -v >> gencore-primates-fq-out.txt 2>&1

        python3 ${PHYLOWIZARD} primates.fq.set.evol.lvl${l}.phy >> gencore-primates-fq-out.txt 2>&1

    done

    mv *.phy $WORK_DIR/fq-primates
    mv *.newick $WORK_DIR/fq-primates
    mv gencore-primates-fq-out.txt $WORK_DIR/fq-primates

    cd -
fi

cd $WORK_DIR
