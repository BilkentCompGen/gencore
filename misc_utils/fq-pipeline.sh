# Required packages:
#   biopython (pip3)
#   ete3 (pip3)
#   six (pip3)
#   PyQt5 (pip3)
# Required scripts, files, and directories:
#   ├── mGorGor
#       ├── m54329U_210319_174352.hifi_reads.fq.gz
#       ├── m54329U_211102_230231.hifi_reads.fq.gz
#       ├── m54329U_211107_082940.hifi_reads.fq.gz
#       ├── m64076_210208_175256.hifi_reads.fq.gz
#       ├── m64076_210213_010909.hifi_reads.fq.gz
#       ├── m64076_210214_073735.hifi_reads.fq.gz
#       ├── m64076_210215_140546.hifi_reads.fq.gz
#       ├── m64076_210326_192259.hifi_reads.fq.gz
#       ├── m64076_230112_193924-bc2013.hifi_reads.fq.gz
#       ├── m64076_230114_050324-bc2013.hifi_reads.fq.gz
#       └── m64076_230115_145043-bc2013.hifi_reads.fq.gz
#   ├── mPanPan1
#       ├── m54329U_210509_045858.hifi_reads.fq.gz
#       ├── m54329U_210510_223954.hifi_reads.fq.gz
#       ├── m54329U_210512_043859.hifi_reads.fq.gz
#       ├── m54329U_211104_082916.hifi_reads.fq.gz
#       ├── m64076_210509_004533.hifi_reads.fq.gz
#       ├── m64076_210511_192234.hifi_reads.fq.gz
#       ├── m64076_210802_234415.hifi_reads.fq.gz
#       ├── m64076_210804_194759.hifi_reads.fq.gz
#       ├── m64076_210806_064833.hifi_reads.fq.gz
#       └── m64076_210807_174839.hifi_reads.fq.gz
#   ├── mPanTro3
#       ├── m54329U_220226_122930.hifi_reads.fq.gz
#       ├── m54329U_220304_132403.hifi_reads.fq.gz
#       ├── m64076_210810_005444.hifi_reads.fq.gz
#       ├── m64076_210813_021703.hifi_reads.fq.gz
#       ├── m64076_210814_131501.hifi_reads.fq.gz
#       ├── m64076_210816_001611.hifi_reads.fq.gz
#       └── m64076_220303_022914.hifi_reads.fq.gz
#   ├── mPonAbe1
#       ├── m54329U_210228_013345.hifi_reads.fq.gz
#       ├── m54329U_210404_020346.hifi_reads.fq.gz
#       ├── m54329U_220227_232630.hifi_reads.fq.gz
#       ├── m54329U_220313_173323.hifi_reads.fq.gz
#       ├── m64076_210219_011735.hifi_reads.fq.gz
#       ├── m64076_210221_195503.hifi_reads.fq.gz
#       ├── m64076_210223_194521.hifi_reads.fq.gz
#       ├── m64076_210225_020019.hifi_reads.fq.gz
#       ├── m64076_210330_204128.hifi_reads.fq.gz
#       ├── m64076_220312_125202.hifi_reads.fq.gz
#       └── m64076_220313_224245.hifi_reads.fq.gz
#   ├── mPonPyg2
#       ├── m54329U_210506_061718.hifi_reads.fq.gz
#       ├── m54329U_220303_022849.hifi_reads.fq.gz
#       ├── m64076_210430_224715.hifi_reads.fq.gz
#       ├── m64076_210505_002126.hifi_reads.fq.gz
#       ├── m64076_210506_062052.hifi_reads.fq.gz
#       ├── m64076_210507_194603.hifi_reads.fq.gz
#       └── m64076_220226_122604.hifi_reads.fq.gz
#   ├── mSymSyn1
#       ├── m54329U_210828_003431.hifi_reads.fq.gz
#       ├── m54329U_210829_112929.hifi_reads.fq.gz
#       ├── m54329U_210901_104742.hifi_reads.fq.gz
#       ├── m54329U_210902_233521.hifi_reads.fq.gz
#       ├── m54329U_210904_103226.hifi_reads.fq.gz
#       ├── m54329U_210905_202241.hifi_reads.fq.gz
#       ├── m54329U_210910_231056.hifi_reads.fq.gz
#       ├── m64076_211107_004425.hifi_reads.fq.gz
#       ├── m64076_211111_101817.hifi_reads.fq.gz
#       └── m64076_211112_233105.hifi_reads.fq.gz

GENCORE_DIR=
FQ_PRIMATES_DIR=

rm -rf fq-primates

mkdir fq-primates
cd fq-primates

ls -1 "${FQ_PRIMATES_DIR}" | sed "s|^|${FQ_PRIMATES_DIR}/|" > input.txt

echo "Gorilla" > shortnames.txt
echo "Bonobo" >> shortnames.txt
echo "Chimpanzee" >> shortnames.txt
echo "Sum Orang" >> shortnames.txt
echo "Bor Orang" >> shortnames.txt
echo "Siamang" >> shortnames.txt

/bin/time -v ${GENCORE_DIR}/gencore fq \
    -i input.txt \
    -s shortnames.txt \
    -l 4 \
    -t 6 \
    -p primates.fq \
    --min-cc 32 \
    -v > gencore-primates-fq-out.txt 2>&1

python3 ${GENCORE_DIR}/phylowizard.py primates.fq.set.evol.lvl4.phy >> gencore-primates-fq-out.txt 2>&1

/bin/time -v ${GENCORE_DIR}/gencore fq \
    -i input.txt \
    -s shortnames.txt \
    -l 5 \
    -t 6 \
    -p primates.fq \
    --min-cc 32 \
    -v >> gencore-primates-fq-out.txt 2>&1

python3 ${GENCORE_DIR}/phylowizard.py primates.fq.set.evol.lvl5.phy >> gencore-primates-fq-out.txt 2>&1

/bin/time -v ${GENCORE_DIR}/gencore fq \
    -i input.txt \
    -s shortnames.txt \
    -l 6 \
    -t 6 \
    -p primates.fq \
    --min-cc 32 \
    -v >> gencore-primates-fq-out.txt 2>&1

python3 ${GENCORE_DIR}/phylowizard.py primates.fq.set.evol.lvl6.phy >> gencore-primates-fq-out.txt 2>&1

cd ..