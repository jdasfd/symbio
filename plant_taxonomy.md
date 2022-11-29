# RLK identification

## Preparation

Visit [Anaconda](https://www.anaconda.com/) to download install script.

- BUSCO

Visit [BUSCO userguide](https://busco.ezlab.org/busco_userguide.html#conda-package) for database.

```bash
# install anaconda
bash Anaconda3-2022.05-Linux-x86_64.sh
# enter to allow and complete installation
# start linux again, directly activating base

# avoid directly activating base when start linux
conda config --set auto_activate_base false
# conda install busco
conda create -n genome -c conda-forge -c bioconda busco=5.4.2
# activate env
conda activate genome
```

- OrthoFinder

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install orthofinder
conda install orthofinder
# primary_transcript.py to get the proteins from the primary transcripts
```

- HMMER

```bash
mkdir -p ~/data/symbio/HMM/PFAM
cd ~/data/symbio/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

hmmpress Pfam-A.hmm
```

- Other software

```bash
brew install dos2unix
```

## CDS or protein selection

### Genomes selection and CDS renaming

There are something to declare: 

- Multiple transcripts were not excluded here, which means different types of the species transcripts were not specified. Maybe the primary transcripts would be adopted afterwards.
- All genomes were put into directories named by `species`, then moved into directory `GENOMES`. Every directory basically had four files including genome, gff, cds and pep.
- Download directly from URL (no records in database)
  - fagopyrum_esculentum
  - allium_sativum
  - elaeis_guineensis
  - passiflora_edulis
  - daemonorops_jenkinsiana
  - sphagnum_magellanicum
  - adiantum_capillus-veneris (no longest gff)
  - ceratopteris_richardii, marchantia_paleacea
  - carica_papaya
  - taxus_chinensis
  - sphagnum_fallax
- Discard
  - punica_granatum (gff has problem)
  - picea_abies (ftp cannot be accessed)
  - phoenix_dactylifera

```bash
mkdir ~/data/symbio/info
cd ~/data/symbio/info

# manually pick some genomes into taxo_genom.xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f taxo_genom.xlsx |
    tsv-select -d, -f 1 |
    sed 1d |
    perl -nl -e '$_ = lc $_;print $_' \
    > genome.lst

wc -l genome.lst
#148 genome.lst

# process genomes
cd ~/data/symbio/GENOMES

for dir in $(ls)
do
    echo "==>{dir}"
    faops size ${dir}/genome.cds |
        cut -f 1 |
        awk -v spe=${dir} '{print ($0"\t"spe_$0)}' \
        > tmp
    faops replace ${dir}/genome.cds tmp ../CDS/${dir}.cds.fa
    rm tmp
done

# there would have a question: faops size will automatically stop at the first space
# so some cds files would be named incorrectly
# check whether changed successfully
cd ~/data/symbio

cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo "==>{}"
        faops size CDS/{}.cds.fa | cut -f 1 | head -n 1
    '
# manually check those cds done wrong

# manually change those species with problems
# following is an example for dealing with the `.||.||.` type cds
cd ~/data/symbio/GENOMES/<species>

# <species> means you should replace to the specific species
cat genome.cds |
    perl -nle '
        print uc $_ if /^[AGCTN]/i;
        print "><species>_$1" if /^>.+\|(.+?)\|\|-*\d\|\|CDS/;
    ' \
    > ../../CDS/<species>.cds.fa 

# finally check the renamed cds numbers
cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo "==>{}"
        CDS=$(cat CDS/{}.cds.fa | faops size stdin | wc -l)
        GEN=$(cat GENOMES/{}/genome.cds | faops size stdin | wc -l)
        if [[ $CDS == $GEN ]]; then
            echo check ok
        else
            echo check wrong
        fi
    '
# nothing wrong
```

### Protein processing

CDS always failed in grouping orthogroups correctly, so I chose protein files.

```bash
mkdir ~/data/symbio/PROTEIN
cd ~/data/symbio

cat info/genome.lst |
    parallel -j 6 --keep-order '
        if [ ! -f GENOMES/{}/genome.pep ]; then
            echo "==>{} without protein"
        else
            cp GENOMES/{}/genome.pep PROTEIN/{}.pep.fa
        fi
    '
# remove all species without pro seq in the genome.lst

conda activate genome

for file in $(ls PROTEIN)
do
    echo "==> ${file}"
    python /home/j/anaconda3/envs/genome/bin/primary_transcript.py PROTEIN/${file}
    echo
done

cd PROTEIN
ls primary_transcripts | wc -l
#148

mv primary_transcripts/ ..
cd ..

# make sure which file should change seq_id
# check seq_id could not be change
for file in $(ls primary_transcripts)
do
    id=$(cat primary_transcripts/${file} | grep '^>' | head -n 1)
    result=$(echo $id | grep '\s')
    if [[ $result != "" ]]; then
        echo "==>${file} name should change"
        echo $id
        echo
    fi
done
#persea_americana
#lemna_minor
#ceratophyllum_demersum
#cat primary_transcripts/<species>.pep.fa |
#    perl -nle '
#        print ">$1" if /^>.+?\|\|.+?\|\|.+?\|\|.+?\|\|(.+?)\|/;
#        print if /^[A-Z]/i;
#    ' > tmp && mv tmp primary_transcripts/<species>.pep.fa

for file in $(ls primary_transcripts)
do
    id=$(cat primary_transcripts/${file} | grep '^>' | head -n 1)
    result=$(echo $id | grep '\s')
    if [[ $result != "" ]]; then    
        cat primary_transcripts/${file} | perl -nle '
            print ">$1" if /^>(.+?)\s/;
            print if /^[A-Z]/i;
        ' \
        > tmp && mv tmp primary_transcripts/${file}
    fi
done
# run above code for sure again

# add species name before the id
for name in $(ls primary_transcripts | perl -p -e 's/\..+$//')
do
    faops size primary_transcripts/${name}.pep.fa |
        cut -f 1 |
        awk -v name=${name} '{print ($0"\t"name"_"$0)}' \
        > name.tmp
    faops replace primary_transcripts/${name}.pep.fa name.tmp tmp.fa \
    && mv tmp.fa primary_transcripts/${name}.pep.fa
done

# some protein files using . for termination codon
# change them to * for satisfying diamond
for file in $(ls primary_transcripts)
do
    SEQ=$(faops some -l 0 primary_transcripts/${file} \
            <(faops size primary_transcripts/${file}) \
            stdout | grep -v '^>')
    result=$(echo $SEQ | grep '\.')
    if [[ $result != "" ]]; then
        echo "==>${file} end is ."
    fi
done

# manually change them
for file in $(ls primary_transcripts)
do
    cat primary_transcripts/${file} |
        perl -nle 'if(/^>/){print;}else{$_ =~ s/\./\*/g; print $_;}' > tmp \
        && mv tmp primary_transcripts/${file}
done
# run above code for sure again

rsync -avP primary_transcripts/ name@ip:jyq/data/symbio/
```

## BUSCO

Running BUSCO for checking the genomes or assemblies sequencing quality.

```bash
mkdir ~/data/symbio/BUSCO
cd ~/data/symbio

# check busco of every genome
bash scripts/busco.sh
```

## OrthoFinder identifying gene orthogroups

### OrthoFinder processing

```bash
# do following to transmit data to server
#rsync -avP /home/j/data/symbio/CDS name@ip:jyq/data/

# connect to server
# ssh name@ip

cd ~/jyq/data/symbio

# orthofinder in conda base env
conda activate

orthofinder -f ./primary_transcripts -og
# orthofinder [options] -f <dir>
# -og: Stop after inferring orthogroups

# cp results to local
#rsync -avP name@ip:jyq/data/symbio/primary_transcripts/OrthoFinder/Results_Oct20/Orthogroups ~/data/symbio/

# If hard limit, h > r already, then you just need to increase the soft limit:
ulimit -n 22004

# check the limit
ulimit -Sn

# continue after orthogroups inferred
orthofinder -fg Results_Oct20/ -M msa -X
# -fg <dir>: Start OrthoFinder from pre-computed orthogroups in <dir>
# -M <txt>: Method for gene tree inference. Options 'dendroblast' & 'msa' [Default = dendroblast]
# -X: Don't add species names to sequence IDs
```

### OrthoFinder results

Acquire results from the workstation using rsync and put in the dir `Orthogroups` (command not shown).

```bash
cd ~/data/symbio/Orthogroups
mkdir groups

cat Orthogroups.tsv | cut -f 1 | sed 1d > ortho.lst

wc -l ortho.lst
#135923 ortho.lst
# Totally 135923 orthogroups were identified

# split each orthogroup into a tsv
cat ortho.lst |
    parallel -j 4 --keep-order '
        echo "==> {}"
        bash ../scripts/ortho_extract.sh {} ./groups/
    '

mkdir pro3
# extract 3 longest proteins
cat ortho.lst |
    parallel -j 16 --keep-order '
        echo "==> {}"
        bash ../scripts/pro_extract.sh groups/{}.csv ../primary_transcripts/ {}
        if [ -f {}.fa ]; then
            cat {}.fa |
                faops size stdin |
                sort -nk 2,2 -r |
                head -n 3 |
                cut -f 1 \
                > {}.lst
            faops some {}.fa {}.lst ./pro3/{}.fa
            rm {}.lst {}.fa
        else
            echo problem
        fi
    '
```

### Renaming and extracting

After OrthoFinder processing, those pep files contained `:` in their seq_ids will be converted to `_`. Protein extraction will fail if you do not change your seq_ids.

The pep file of ceratophyllum_demersum has predicted sequences generated from six frame translation, remove them manually (code not shown).

```bash
cd ~/data/symbio

rm info/rename_id.lst
# check those files whose seq_ids should be changed
for file in $(ls primary_transcripts)
do
    if grep -i -q ':' primary_transcripts/${file}; then
        echo "${file}" >> info/rename_id.lst
    fi
done
#anthoceros_angustus.pep.fa
#bonia_amplexicaulis.pep.fa
#calamus_simplicifolius.pep.fa
#cocos_nucifera.pep.fa
#daemonorops_jenkinsiana.pep.fa
#mesotaenium_endlicherianum.pep.fa
#olyra_latifolia.pep.fa
#phalaenopsis_equestris.pep.fa
#phyllostachys_edulis.pep.fa
#rhododendron_delavayi.pep.fa
#welwitschia_mirabilis.pep.fa

for file in $(cat info/rename_id.lst)
do
    cat primary_transcripts/${file} |
        perl -p -e 's/:/_/g' \
        > tmp && mv tmp primary_transcripts/${file}
done

# run again to check whether renaming succeed
for file in $(ls primary_transcripts)
do
    if grep -i -q ':' primary_transcripts/${file}; then
        echo "${file}"
    fi
done
# OK
```

## RLK identifying

RLK (receptor-like receptor) structures:

From N-terminal ECD (extracellular domain) - TMD (transmembrane domain) - KD (kinase domain) to C-terminal.

So domains should be formatted correctly for a true RLK.

### Identify domains among all proteins via `hmmscan`

```bash
cd ~/data/symbio
mkdir -p ~/data/symbio/DOMAIN/pfam

cat info/genome.lst |
    parallel -j 6 --keep-order '
        echo {}
        hmmscan --cpu 2 -E 0.1 --domE 0.1 -o DOMAIN/pfam/{}.txt \
            ./HMM/PFAM/Pfam-A.hmm primary_transcripts/{}.pep.fa
    '

ls DOMAIN/pfam/*.txt |
    parallel -j 6 --keep-order '
        echo "==> {/.}"
        perl scripts/hmm_results.pl -i {} \
        > DOMAIN/pfam/{/.}.pfam.tsv
    '

ls DOMAIN/pfam/*.txt |
    parallel -j 12 --keep-order 'echo {/.}' \
    > DOMAIN/species.lst
```

- Extract kinase domain

```bash
cd ~/data/symbio/DOMAIN

rm all_domains.tsv
# all domains were extracted
for file in $(ls pfam/*.pfam.tsv)
do
    cat ${file} |
        tsv-select -H -f Domain |
        sed 1d \
        >> all_domains.tsv
done

# basic info of domains identified
cat all_domains.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_domains.tsv

cat all_domains.tsv | wc -l
#19622
# all domains identified, including those repeated domains from a location
# this will be considered later

# extract all domains of pkinase with cutoff E-value <= 1e-4
# be aware of pk may contained domains not only pkinase
ls pfam/*.pfam.tsv |
    perl -p -e 's/\.pfam\.tsv$//' |
    parallel -j 12 --keep-order '
        echo "==> {/}"
        cat {}.pfam.tsv |
            tsv-filter -H --iregex Domain:pk --le E_value:"1e-4" |
            sed 1d |
            tsv-select -f 1 |
            tsv-uniq \
            > KD/{/}.lst
        cat {}.pfam.tsv |
            sed 1d |
            tsv-join -f KD/{/}.lst -k 1 |
            tsv-select -f 1,2,4,5,3 \
            > KD/{/}.tsv
        rm KD/{/}.lst
    '

rm all_KD.tsv
for file in $(ls KD/*.tsv)
do
    cat ${file} |
        tsv-select -f 2 |
        tsv-filter --iregex 1:pk \
        >> all_KD.tsv
done

cat all_KD.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_KD.tsv

rm all_KD_species.tsv
# all proteins with the pk domain among species
ls KD/*.tsv |
    parallel -j 12 --keep-order '
        cat {} |
            cut -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={/.} '\''{print (spe"\t"$0)}'\'' \
            >> all_KD_species.tsv
    '

# all proteins
cat all_KD_species.tsv |
    cut -f 2 |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#204137

# extract all protein sequences contained the kinase domain
mkdir -p ~/data/symbio/DOMAIN/seq_KD

ls KD/*.tsv |
    parallel -j 12 --keep-order '
        echo "==> {/.}"
        faops some -l 0 \
            ../primary_transcripts/{/.}.pep.fa \
            <(cat {} | cut -f 1 | tsv-uniq) \
            ./seq_KD/{/.}.fa
    '

# check extraction
wc -l ./seq_KD/*.fa | grep 'total' | perl -p -e 's/\s+(\d+).+$/$1\/2/' | bc
#204137
# extraction complete
```

- `TMHMM2` for transmembrane domain identification

Manually upload files to [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) and copy results into a tsv file.

Although `selenium` could deal with this automatically, but there are always problems that I cannot solve. But still, a script will be saved in [scripts](https://github.com/jdasfd/symbio/tree/main/scripts). It will be completed later (maybe). XD

RLK (receptor-like receptor) structures:

From N-terminal ECD (extracellular domain) - TMD (transmembrane domain) - KD (kinase domain) to C-terminal.

So domains should be formatted correctly for a true RLK.

```bash
cd ~/data/symbio/DOMAIN

mkdir -p ~/data/symbio/DOMAIN/TMD
# manually doing this steps:
# upload pep.fa acquired from the previous step to TMHMM 2.0 websites
# copy all results into a tsv

mkdir ~/data/symbio/DOMAIN/RLK

dos2unix TMD/*.tsv

# if there is at least a TMD
# retrieved as potential RLK proteins
ls KD/*.tsv |
    parallel -j 12 --keep-order '
        echo "==> {/.}"
        perl ../scripts/tmhmm_result.pl -i TMD/{/.}.tsv |
            tsv-select -f 1,4,2,3 \
            > RLK/{/.}.tsv
        cat {} |
            tsv-select -f 1,2,3,4 \
            >> RLK/{/.}.tsv
        cat RLK/{/.}.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            tsv-sort > RLK/{/.}.tmp.lst
        if [ -f RLK/{/.}.sort.tsv ]; then
            rm RLK/{/.}.sort.tsv
        fi
        for acc in $(cat RLK/{/.}.tmp.lst)
        do
            cat RLK/{/.}.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-sort -nk 3,3 \
                >> RLK/{/.}.sort.tsv
        done
        rm RLK/{/.}.tmp.lst
    '

ls RLK/*.sort.tsv |
    sed 's/\.sort\.tsv$//' |
    parallel -j 12 --keep-order '
        echo "==> {/}"
        cat {}.sort.tsv |
            cut -f 1 |
            uniq |
            wc -l |
            awk -v spe={/} '\''{print (spe"\t"$0)}'\'' \
            >> all_RLK_pot_species.tsv
    '

# extract RLK to .lst according to RLK structure
ls RLK/*.sort.tsv |
    sed 's/\.sort\.tsv$//' |
    parallel -j 12 --keep-order '
        echo "==> {/}"
        cat {}.sort.tsv |
            perl ../scripts/RLK.pl \
            > {}.RLK.lst
        cat {}.sort.tsv |
            tsv-join -f {}.RLK.lst -k 1 \
            > {}.RLK.tsv
        rm {}.RLK.lst
    '

ls RLK/*.RLK.tsv |
    sed 's/\.RLK\.tsv$//' |
    parallel -j 12 --keep-order '
        echo "==> {/}"
        cat {}.RLK.tsv |
            cut -f 1 |
            uniq |
            wc -l |
            awk -v spe={/} '\''{print (spe"\t"$0)}'\'' \
            >> all_RLK_true_species.tsv
    '
```



```

## RNA-seq of AM symbiosis

```bash
mkdir ~/data/symbio/sra
cd ~/data/symbio/sra


```

```bash
for file in $(ls)
do
    hmmscan --cpu 8 -E 1e-4 --domE 1e-4 -o ${file}.txt \
    ../../symbio/HMM/PFAM/Pfam-A.hmm ${file}
done
```
