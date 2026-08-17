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

RUN_CORALLEL_TUMOR="true"
RUN_MASH_TUMOR="true"
RUN_SOURMASH_TUMOR="true"
RUN_DASHING2_TUMOR="true"
RUN_ANI="true"
RUN_CORALLEL_PRIM="true"
RUN_MASH_PRIM="true"
RUN_SOURMASH_PRIM="true"
RUN_DASHING2_PRIM="true"
RUN_CORALLEL_FQ="true"
RUN_CORALLEL_FQ_CC="true"

CORALLEL="corallel"
PHYLOWIZARD="phylowizard.py"
PLOT_ANI="misc_utils/plot-mut-vs-dist.py"
MASH="mash"
CONVERT="misc_utils/convert.py"
SOURMASH="sourmash"
DASHING2="dashing2"
CMP_TREE="misc_utils/cmp_trees.py"
PLOT_DIST="misc_utils/plot_tree_disruption.py"

PRIMATES_FQ_DIR=.

WORK_DIR=.

MUTATION_COUNT_PHY="mutation.count.phy"

CONFIG_FILE="corallel-config.sh"

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
fi

cd $WORK_DIR

# Corallel Tumor
if [ "$RUN_CORALLEL_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f corallel-tumor-fa-out.txt

    for l in 4 5 6 7 8 9; do

        /bin/time -v ${CORALLEL} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 4 \
            -l "$l" \
            -p tumor \
            -v >> corallel-tumor-fa-out.txt 2>&1

        python3 ${PHYLOWIZARD} tumor.set.evol.lvl${l}.phy --normalize >> corallel-tumor-fa-out.txt 2>&1

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

# Dashing2 Tumor
if [ "$RUN_DASHING2_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f dashing2-tumor-fa-out.txt

    for s in 1000 5000 50000 500000 5000000; do

        /bin/time -v ${DASHING2} sketch \
            -F filenames.txt \
            -p 4 \
            -S${s} \
            --cmpout "tumor.dashing2.s${s}.phy" >> dashing2-tumor-fa-out.txt 2>&1

        python3 ${CONVERT} --dashing2 tumor.dashing2.s${s}.phy tumor.dashing2.s${s}.norm.phy

        if [ -f "tumor.dashing2.s${s}.norm.phy" ]; then
            python3 ${PHYLOWIZARD} tumor.dashing2.s${s}.norm.phy --normalize >> dashing2-tumor-fa-out.txt 2>&1
        else
            echo "tumor.dashing2.s${s}.norm.phy not found!"
        fi

    done

    cd ..
fi

# Corallel-ANI tumor
if [ "$RUN_ANI" = "true" ]; then

    mkdir -p mut-vs-ani-plots
    cd mut-vs-ani-plots

    if [ -f "${MUTATION_COUNT_PHY}" ]; then
        for l in 4 5 6 7 8 9; do
            python3 ${PLOT_ANI} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.jaccard.lvl${l}.phy --label "Corallel (${l})" -o count-corallel-j-${l}.pdf
            python3 ${PLOT_ANI} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.evol.lvl${l}.phy --label "Corallel (${l})" -o count-corallel-p-${l}.pdf
        done
    fi

    cd ..
fi

# Corallel primates
if [ "$RUN_CORALLEL_PRIM" = "true" ]; then

    cd fa-primates

    rm -f corallel-primates-fa-out.txt

    # Corallel
    for l in 4 5 6 7; do
        /bin/time -v ${CORALLEL} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 7 \
            -l "$l" \
            -p primates \
            -v >> corallel-primates-fa-out.txt 2>&1;

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

# Dashing2 primates
if [ "$RUN_DASHING2_PRIM" = "true" ]; then

    cd fa-primates

    rm -f dashing2-primates-fa-out.txt

    /bin/time -v ${DASHING2} sketch \
        -F filenames.txt \
        -p 6 \
        --cmpout "primates.dashing2.default.phy" >> dashing2-primates-fa-out.txt 2>&1

    python3 ${CONVERT} --dashing2 primates.dashing2.default.phy primates.dashing2.default.norm.phy
    python3 ${PHYLOWIZARD} primates.dashing2.default.norm.phy >> dashing2-primates-fa-out.txt 2>&1

    cd ..
fi

# Corallel primates fq
if [ "$RUN_CORALLEL_FQ" = "true" ]; then

    mkdir -p fq-primates

    cd $PRIMATES_FQ_DIR

    rm -f corallel-primates-fq-out.txt

    for l in 4 5 6; do

        /bin/time -v ${CORALLEL} fq \
            -i filenames.txt \
            -s shortnames.txt \
            -l "$l" \
            -t 32 \
            -r 11 \
            -p primates.fq \
            -v >> corallel-primates-fq-out.txt 2>&1

        python3 ${PHYLOWIZARD} primates.fq.set.evol.lvl${l}.phy >> corallel-primates-fq-out.txt 2>&1

    done

    mv *.phy $WORK_DIR/fq-primates
    mv *.newick $WORK_DIR/fq-primates
    mv corallel-primates-fq-out.txt $WORK_DIR/fq-primates

    cd -

    cd ..
fi

# Corallel primates fq
if [ "$RUN_CORALLEL_FQ_CC" = "true" ]; then

    mkdir -p fq-primates-cc

    cd $PRIMATES_FQ_DIR

    for l in 4 5 6; do
        /bin/time -v ${CORALLEL} fq \
            -i filenames.txt \
            -s shortnames.txt \
            -o binnames${l}.txt \
            -l "$l" \
            -t 32 \
            -r 6 \
            -p primates.fq \
            -v >> corallel-primates-fq-out.txt 2>&1
    done

    for l in 4 5 6; do
        /bin/time -v ${CORALLEL} ld \
            -i binnames${l}.txt \
            -s shortnames.txt \
            -t 6 \
            -r 6 \
            -l "$l" \
            -p primates.fq -v
        
        rm -f distances.lvl${l}.txt

        for cc in {0..200..4}; do
            python3 ${PHYLOWIZARD} primates.fq.cc${cc}.set.evol.lvl${l}.phy >> corallel-primates-fq-out.txt 2>&1
            python3 ${CMP_TREE} ${PRIMATES_GROUND} primates.fq.cc${cc}.set.evol.lvl${l}.nj.newick >> distances.lvl${l}.txt
            echo "" >> distances.lvl${l}.txt
        done

        python  distances.lvl${l}.txt -o tree_disruption.lvl${l}.pdf
    done

    mv *.phy $WORK_DIR/fq-primates-cc
    mv *.newick $WORK_DIR/fq-primates-cc
    mv corallel-primates-fq-out.txt $WORK_DIR/fq-primates-cc
    mv distances.lvl*.txt $WORK_DIR/fq-primates-cc
    mv tree_disruption.lvl*.pdf $WORK_DIR/fq-primates-cc
    
    cd -

    cd ..
fi

cd $WORK_DIR
