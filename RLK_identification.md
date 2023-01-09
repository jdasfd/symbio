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
brew install wang-q/tap/nwr

cd ~/Scripts
git clone https://github.com/wang-q/fig_table.git

# script need python2
sudo apt install python2.7
```

## Genomes selection and processing

There are something to declare: 

- Multiple transcripts were not excluded here, which means different types of the species transcripts were not specified.
- All genomes were put into directories named by their own names of species, then moved into directory `GENOMES`. Most species directory basically had four files including genome, gff, cds and pep (some species).
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
  - All algae
- All inforamation of species is recorded into `info/taxo_genom.xlsx`. Taxonomy (Col2) could be checked through code below.

### Collect genomes information

All genomes info were collected into `taxo_genom.xlsx`.

```bash
mkdir ~/data/symbio/info
cd ~/data/symbio/info

# manually pick some genomes into taxo_genom.xlsx

bash taxonomy.sh <species_name>
# taxonomy.sh relied on nwr
# col2 in taxo_genom.xlsx were acquired from the script

perl ~/Scripts/fig_table/xlsx2csv.pl -f taxo_genom.xlsx |
    tsv-select -d, -f 1 |
    sed 1d |
    perl -nl -e '$_ = lc $_;print $_' \
    > genome.lst

wc -l genome.lst
#177 genome.lst
```

### Find primary transcripts

Using [getLongestProteinFromGFF.py](https://github.com/MWSchmid/Ngou-et-al.-2022/blob/master/scripts/getLongestProteinFromGFF.py) from Ngou acquiring longest transcripts.

```bash
mkdir ~/data/symbio/PROTEINS
cd ~/data/symbio

cat info/genome.lst |
    parallel -j 16 --keep-order '
        echo "==> {}"
        python2.7 scripts/getLongestProteinFromGFF.py \
            GENOMES/{}/genome.pep GENOMES/{}/genome.gff \
            PROTEINS/{}.longest.pep
    '

# check primary results
cat info/genome.lst |
    parallel -j 16 -k '
        echo "==> {}"
        before=$(
            cat GENOMES/{}/genome.pep |
            grep "^>" |
            wc -l
            )
        after=$(
            cat PROTEINS/{}.longest.pep |
            grep "^>" |
            wc -l
        )
        echo -e "{}\t${before}\t${after}"
    '

# Abnormal results:
#==> persea_americana
#persea_americana        24616   1
#==> lemna_minor
#lemna_minor     22382   1
#==> ceratophyllum_demersum
#ceratophyllum_demersum  25373   1

for genome in persea_americana lemna_minor ceratophyllum_demersum
do
    cat GENOMES/${genome}/genome.pep |
        perl -nle '
            if($_ =~ /^>/){
                @name = split /\|\|/, $_;
                print ">$name[4]";
            }
            else{
                print uc $_;
            }
        ' \
        > GENOMES/${genome}/genome.rename.pep
done

for genome in persea_americana lemna_minor ceratophyllum_demersum
do
    python2.7 scripts/getLongestProteinFromGFF.py \
        GENOMES/${genome}/genome.rename.pep \
        GENOMES/${genome}/genome.gff \
        PROTEINS/${genome}.longest.pep
done

# ceratophyllum_demersum has six box translation
# all gene_id recorded in info/rm_ceratophyllum_demersum.lst
# rm them
faops some PROTEINS/ceratophyllum_demersum.longest.pep \
    <(cat PROTEINS/ceratophyllum_demersum.longest.pep |
        faops size stdin |
        cut -f 1 |
        tsv-join -f info/rm_ceratophyllum_demersum.lst -k 1 -e
    ) tmp && mv tmp PROTEINS/ceratophyllum_demersum.longest.pep

ls PROTEINS | wc -l
#177
# All primary transcripts were successfully extracted
```

### Protein processing

Add taxonomy and species name before all proteins and make them satisfying the next step.

```bash
cd ~/data/symbio

perl ~/Scripts/fig_table/xlsx2csv.pl -f info/taxo_genom.xlsx |
    tsv-select -d ',' -f 1,3 |
    sed 1d |
    perl -nla -F',' -e '
        $F[0] =~ s/.*/\L$&/g;
        print "CL_$F[0]" if $F[1] =~ /^Chloro/;
        print "CA_$F[0]" if $F[1] =~ /^Charo/;
        print "BR_$F[0]" if $F[1] =~ /^Bryo/;
        print "LY_$F[0]" if $F[1] =~ /^Lyco/;
        print "FE_$F[0]" if $F[1] =~ /^Fern/;
        print "GY_$F[0]" if $F[1] =~ /^Gymno/;
        print "AN_$F[0]" if $F[1] =~ /^Angio/;
    ' \
    > info/name.lst

# add species name before the id
cat info/genome.lst |
    parallel -j 12 -k '
        new_name=$(cat info/name.lst | grep "{}")
        faops size PROTEINS/{}.longest.pep |
            cut -f 1 |
            awk -v name=${new_name} '\''{print ($0"\t"name"_"$0)}'\'' \
            > PROTEINS/{}.name.tmp
        faops replace PROTEINS/{}.longest.pep \
            PROTEINS/{}.name.tmp \
            PROTEINS/{}.tmp.fa \
            && mv PROTEINS/{}.tmp.fa PROTEINS/{}.longest.pep
        rm PROTEINS/{}.name.tmp
    '

# some protein files using . for termination codon
# change them to * for satisfying diamond
for file in $(ls PROTEINS)
do
    SEQ=$(faops some -l 0 PROTEINS/${file} \
            <(faops size PROTEINS/${file}) \
            stdout | grep -v '^>')
    result=$(echo $SEQ | grep '\.')
    if [[ $result != "" ]]; then
        echo "==>${file} end is ."
    fi
done

# manually change them
cat info/genome.lst |
    parallel -j 16 -k '
        echo "==> {}"
        cat PROTEINS/{}.longest.pep |
            perl -nle '\''
                if (/^>/){
                    print;
                }
                else{
                    $_ =~ s/\./\*/g;
                    print $_;
                }
            '\'' > PROTEINS/{}.tmp \
                && mv PROTEINS/{}.tmp PROTEINS/{}.longest.pep
    '

# run above code for sure again

# upload to the workstation
rsync -avP PROTEINS name@ip:jyq/data/symbio/
```

## BUSCO

Running BUSCO for checking the genomes or assemblies sequencing quality.

```bash
mkdir ~/data/symbio/BUSCO
cd ~/data/symbio

# check busco of every genome
bash scripts/busco.sh
```

## Identifying gene orthogroups via OrthoFinder

### OrthoFinder processing

```bash
cd ~/data/symbio

# orthofinder in conda base env
conda activate

orthofinder -f ./PROTEINS -og
# orthofinder [options] -f <dir>
# -og: Stop after inferring orthogroups

# continue after orthogroups inferred
#orthofinder -fg Results_Oct20/ -M msa -X
# -fg <dir>: Start OrthoFinder from pre-computed orthogroups in <dir>
# -M <txt>: Method for gene tree inference. Options 'dendroblast' & 'msa' [Default = dendroblast]
# -X: Don't add species names to sequence IDs

# orthofinder will make a working directory under the -f, which is named by the date
#cp -r ./PROTEINS/OrthoFinder/Results_Dec22/Orthogroups ./
#rsync -avP name@ip:jyq/data/symbio/Orthogroups ~/data/symbio/
```

### OrthoFinder results

```bash
cd ~/data/symbio/Orthogroups
mkdir groups

cat Orthogroups.tsv | cut -f 1 | sed 1d > ortho.lst

wc -l ortho.lst
#150270 ortho.lst
# Totally 150270 orthogroups were identified

# split each orthogroup into a tsv
cat ortho.lst |
    parallel -j 12 -k '
        echo "==> {}"
        bash ../scripts/ortho_extract.sh {} ./groups/
    '

# combine groups .tsv files into one
# col1: ortho_name, col2: prot_name
cat ortho.lst |
    parallel -j 1 -k '
        echo "==> {}"
        cat groups/{}.csv |
            mlr --icsv --otsv cat |
            cut -f 2 |
            sed 1d |
            tr "," "\n" |
            sed '\''s/^\s//'\'' |
            sed '\''/^$/d'\'' |
            awk -v ortho={} '\''{print ($0"\t"ortho)}'\'' \
            >> all_pro_ortho.tsv
    '

wc -l all_pro_ortho.tsv
#5450693 all_pro_ortho.tsv

dos2unix Orthogroups.GeneCount.tsv

cat Orthogroups.GeneCount.tsv |
    tsv-select -H -f Orthogroup,Total \
    > ortho_gene.tsv

cat ortho_gene.tsv | tsv-summarize -H --sum Total | sed 1d
#5450693

# rm -rf ./groups
```

### Count the length and the number of proteins

There are few things to be declared:

- Length

```bash
cd ~/data/symbio
mkdir -p DOMAIN/protein

# all proteins' length
ls PROTEINS/*.longest.pep |
    sed 's/\.longest\.pep$//' |
    parallel -j 12 -k '
        echo "==> {/}"
        faops size {}.longest.pep \
        > DOMAIN/protein/{/}.length.tsv
    '

cat DOMAIN/protein/*.length.tsv | wc -l
#5800767
```

## Identify RLK

RLK (receptor-like receptor) structures:

From N-terminal ECD (extracellular domain) - TMD (transmembrane domain) - KD (kinase domain) to C-terminal.

So domains should be formatted correctly for a true RLK.

### Identify domains among all proteins via `hmmscan`

```bash
cd ~/data/symbio
mkdir -p DOMAIN/pfam

cat info/genome.lst |
    parallel -j 10 -k '
        echo {}
        hmmscan --cpu 2 -E 0.1 --domE 0.1 -o DOMAIN/pfam/{}.txt \
            ./HMM/PFAM/Pfam-A.hmm PROTEINS/{}.longest.pep
    '

ls DOMAIN/pfam/*.txt | wc -l
#177

# the number of all query sequences
ls DOMAIN/pfam/*.txt |
    parallel -j 12 -k '
        cat {} |
            grep -v '^#' |
            grep '^Query' |
            grep '\]$' |
            wc -l
        ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#5747688

ls DOMAIN/pfam/*.txt |
    parallel -j 12 -k '
        echo "==> {/.}"
        perl scripts/hmm_results.pl -i {} \
        > DOMAIN/pfam/{/.}.pfam.tsv
    '
# screen will show this command:
#Missed this line:    [No individual domains that satisfy reporting thresholds (although complete target did)]
# It means that those domains could not be identified among proteins

# check whether pfam results contained all proteins
ls DOMAIN/pfam/*.pfam.tsv |
    parallel -j 12 --k '
        cat {} |
            tsv-select -H -f QUERY |
            sed 1d |
            tsv-uniq |
            wc -l
    ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#4468845

# Results without domains
ls DOMAIN/protein/*.length.tsv |
    sed 's/\.length\.tsv$//' |
    parallel -j 12 -k '
        echo "==>{/}"
        cat {}.length.tsv |
            tsv-join -f <(
                cat DOMAIN/pfam/{/}.pfam.tsv |
                    sed 1d |
                    cut -f 1
                ) \
                -k 1 -e \
            > DOMAIN/protein/{/}.no_hmm.lst
    '

# The number of all proteins without domains scanned
ls DOMAIN/protein/*.no_hmm.lst |
    parallel -j 12 -k '
        cat {} |
            tsv-select -f 1 |
            tsv-uniq |
            wc -l
    ' |
    tr '\n' '+' |
    sed 's/+$/\n/' |
    bc
#1331922

# echo 4468845+1331922 | bc
#5800767
```

### Extract pkinase domain

- Kinase domains identification
 
RLK domains identification process was referred to the [article](https://www.nature.com/articles/s41477-022-01260-5). The kinase domain was identified as PFAM [PF00069](https://www.ebi.ac.uk/interpro/entry/pfam/PF00069/) named Protein kinase domain (short name *Pkinase*).

Therefore, all domains contained `pk` were collected here. This step was to determine whether there were other kinase domains named not pkinase. The reason why doing this here is that we use `hmmscan` for domains scanning across the whole protein sequences by `PFAM-A` database. So double check them.

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
# all domains identified, including those domains scanned repeatedly from a location
# this will be considered later

# extract all domains contained pk with cutoff E-value <= 1e-5
rm all_pk.tsv
for species in $(cat ../info/genome.lst)
do
    echo "==> ${species}"
    cat pfam/${species}.pfam.tsv |
        tsv-filter -H --or \
            --iregex Domain:pk \
            --iregex Domain:kinase |
        tsv-filter -H --le E_value:"1e-5" |
        sed 1d |
        tsv-select -f 2 \
    >> all_pk.tsv
done

cat all_pk.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_pk.tsv
```

True pkinase domain were kept after searching on [InterPro](https://www.ebi.ac.uk/interpro/)

- True pkinase domains
  - Pkinase
  - PK_Tyr_Ser-Thr
  - Pkinase_fungal
  - Pkinase_C

E-value cutoff <= 1e-5 for all KDs (pkinase domain) - kinase domains are relatively more conserved than ECD.

```bash
cd ~/data/symbio/DOMAIN
mkdir -p ~/data/symbio/DOMAIN/KD

# extract all domains contained pk with cutoff E-value <= 1e-5
cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat pfam/{}.pfam.tsv |
            tsv-filter -H --or \
                --str-eq Domain:Pkinase \
                --str-eq Domain:Pkinase_fungal \
                --str-eq Domain:PK_Tyr_Ser-Thr \
                --str-eq Domain:Pkinase_C |
            tsv-filter -H --le E_value:"1e-5" |
            sed 1d |
            tsv-select -f 1 |
            tsv-uniq \
            > KD/{}.lst
        cat pfam/{}.pfam.tsv |
            sed 1d |
            tsv-join -f KD/{}.lst -k 1 |
            tsv-select -f 1,2,4,5,3 \
            > KD/{}.tsv
    '

wc -l ./KD/*.lst | grep 'total'
# 194150 total

rm all_true_KD.tsv
for species in $(cat ../info/genome.lst)
do
    echo "==> ${species}"
    cat KD/${species}.tsv |
        tsv-filter --or \
            --str-eq 2:Pkinase \
            --str-eq 2:Pkinase_fungal \
            --str-eq 2:PK_Tyr_Ser-Thr \
            --str-eq 2:Pkinase_C |
        tsv-filter -H --le 5:"1e-5" |
        tsv-select -f 2 \
    >> all_true_KD.tsv
done

cat all_true_KD.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_true_KD.tsv

# extract all protein sequences contained pk
mkdir -p ~/data/symbio/DOMAIN/seq_KD

cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        faops some -l 0 \
            ../PROTEINS/{}.longest.pep \
            KD/{}.lst \
            seq_KD/{}.fa
    '

# check extraction
wc -l ./seq_KD/*.fa | grep 'total' | perl -p -e 's/\s+(\d+).+$/$1\/2/' | bc
#194150
# extraction complete
```

### Predict transmembrane domain via `TMHMM2`

Manually upload files to [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0) and copy results into a tsv file.

Although `selenium` could deal with this automatically, but there are always problems that I cannot solve. But still, a script will be saved in [scripts](https://github.com/jdasfd/symbio/tree/main/scripts). It will be completed later (maybe XD).

```bash
cd ~/data/symbio/DOMAIN
mkdir -p TMD

cat ../info/genome.lst |
    parallel -j 12 -k '
        touch TMD/{}.tsv
    '
# manually doing this steps:
# upload pep.fa of each specie acquired from the previous step to TMHMM 2.0 websites
# copy all predicted results into a tsv named by species

# tsv-utils could only deal with those files in unix style (LF), convert them using dos2unix
dos2unix TMD/*.tsv

wc -l TMD/* | grep 'total'
#  194150 total
# alright

# retrieved potential pk proteins if they had at least a TMD
cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        perl ../scripts/tmhmm_result.pl -i TMD/{}.tsv |
            tsv-select -f 1,4,2,3 \
            > TMD/{}.TMD.tsv
        '

cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat TMD/{}.tsv |
            cut -f 1 |
            tsv-join -f TMD/{}.TMD.tsv \
            -k 1 -e \
            > TMD/{}.not_TMD.lst
    '
```

### Remove repeated domains

Because of redundant results acquired from `hmmscan` (using PFAM-A database with cutoff E-value 1e-1), a process should be adopted for reducing repeated domains.

The main principle - filtering according to E-value.

Protein domains of each protein were sorted by E-value (from less to more). The domains with smaller E-value were kept, and every repeated domains with larger E-value were removed.

```bash
cd ~/data/symbio/DOMAIN
mkdir -p UNIQ

# sort according to E_value from less to more
cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        if [ -f UNIQ/{}.sort.tsv ]; then
            rm UNIQ/{}.sort.tsv
        fi
        for acc in $(cat KD/{}.lst)
        do
            cat pfam/{}.pfam.tsv |
                tsv-filter -H --str-eq QUERY:${acc} |
                sed 1d |
                tsv-select -f 1,2,3,4,5 |
                sort -gk 3,3 \
                >> UNIQ/{}.sort_e.tsv
        done
    '
# sort -g: use generic numerical value

# a script called domain_uniq to remove repeated domains
cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat UNIQ/{}.sort_e.tsv |
            perl ../scripts/domain_uniq.pl |
            tsv-select -f 1,4,5,2 \
            > UNIQ/{}.uniq.tsv
    '

# combine TMD to uniq.tsv
cat ../info/genome.lst |
    parallel -j 12 -keep-order '
        echo "==> {}"
        for acc in $(cat KD/{}.lst)
        do
            cat TMD/{}.TMD.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-select -f 1,3,4,2 \
                >> UNIQ/{}.uniq.tsv
        done
    '

# sort by start pos
cat ../info/genome.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat UNIQ/{}.uniq.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            tsv-sort > UNIQ/{}.tmp.lst
        if [ -f UNIQ/{}.sort_start.tsv ]; then
            rm UNIQ/{}.sort_start.tsv
        fi
        for acc in $(cat UNIQ/{}.tmp.lst)
        do
            cat UNIQ/{}.uniq.tsv |
                tsv-filter --str-eq 1:${acc} |
                tsv-sort -nk 2,2 \
                >> UNIQ/{}.sort_start.tsv
        done
        rm UNIQ/{}.tmp.lst
    '

rm all_uniq_KD.tsv
for species in $(cat ../info/genome.lst)
do
    cat UNIQ/${species}.sort_start.tsv |
        tsv-select -f 4 |
        tsv-filter --or \
            --str-eq 1:Pkinase \
            --str-eq 1:Pkinase_fungal \
            --str-eq 1:PK_Tyr_Ser-Thr \
            --str-eq 1:Pkinase_C \
        >> all_uniq_KD.tsv
done

cat all_uniq_KD.tsv |
    tsv-summarize -g 1 --count |
    sort -r -nk 2,2 \
    > tmp && mv tmp all_uniq_KD.tsv
# these are all pkinase unrepeated domains
```

### Pick pkinase domains only

All domains contained `pk` were collected. But RLK only contains the pkinase domain.

Those domains in `all_KD_uniq.tsv` were checked through [InterPro](https://www.ebi.ac.uk/interpro/), and true pkinase domain were kept.

- Pkinase
- PK_Tyr_Ser-Thr
- Pkinase_fungal

```bash
cd ~/data/symbio/DOMAIN
mkdir -p RLK

# all pkinase protein names
cat species.lst |
    parallel -j 12 --keep-order '
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-filter --or --str-eq 4:Pkinase \
                --str-eq 4:Pkinase_fungal \
                --str-eq 4:PK_Tyr_Ser-Thr |
            cut -f 1 |
            tsv-uniq \
            > RLK/{}.RLK.lst
    '

# all RLK that truly contained the pkinase domain
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-join -f RLK/{}.RLK.lst -k 1 \
            > RLK/{}.RLK.tsv
    '

echo -e "species\tRLK_num" > all_RLK_species.tsv
cat species.lst |
    parallel -j 12 --keep-order '
        echo "==> {}"
        cat RLK/{}.RLK.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={} '\''{print (spe"\t"$0)}'\'' \
            >> all_RLK_species.tsv
    '

cat all_PK_species.tsv |
    tsv-join -H -f all_RLK_species.tsv -k species -a RLK_num \
    > all_count_species.tsv

perl ~/Scripts/fig_table/xlsx2csv.pl -f ../info/taxo_genom.xlsx |
    mlr --icsv --otsv cat |
    tsv-select -H -f SpeciesName,ComName |
    sed 1d |
    perl -nle 'print lc $_;' |
    sed '1ispecies\ttaxo' \
    > taxonomy.tsv

cat taxonomy.tsv |
    tsv-join -H -f all_count_species.tsv \
        -k species -a pot_PK,true_PK,RLK_num \
    > Species_result.tsv
```

### Count all ECDs

```bash
mkdir -p ~/data/symbio/RLK
cd ~/data/symbio/RLK

cat ../DOMAIN/RLK/*.RLK.lst > all_RLK.lst
cat all_RLK.lst | wc -l
#63888

cat all_RLK.lst |
    tsv-join -f ../Orthogroups/all_pro_ortho.tsv -k 1 -a 2 \
    > RLK_ortho.tsv

cat RLK_ortho.tsv | wc -l
#63766

cat all_RLK.lst |
    tsv-join -f RLK_ortho.tsv -k 1 -e \
    > RLK_without_ortho.lst
cat RLK_without_ortho.lst | wc -l
#122

cat RLK_ortho.tsv |
    tsv-summarize -g 2 --count |
    sort -nk 2,2 -r |
    sed '1iOrthogroup\tRLK_num' \
    > ortho_count.tsv

cat ortho_count.tsv | sed 1d | wc -l
#1486
# All 63766 RLKs were grouped into 1486 orthogroups

cat ortho_count.tsv |
    tsv-join -H -f ../Orthogroups/ortho_gene.tsv \
        -k Orthogroup -a Total \
        > Orthogroup_result.tsv

cat Orthogroup_result.tsv |
    sed 1d |
    perl -nlae '
        $ratio = $F[1]/$F[2];
        printf "%s\t%.3f\n",$F[0], $ratio;
    ' |
    sort -nk 2,2 -r |
    sed '1iOrthogroup\tRatio'\
    > Orthogroup_ratio.tsv

# deprecated due to the linux system issue
# using inside content in the RStudio for substitution
Rscript -e '
    library(ggplot2)
    args <- commandArgs(T)
    ratio <- read_tsv(args[1])
    p <- ggplot(ratio, aes(x = Orthogroup, y = Ratio)) +
         geom_bar(stat = "identity", position = "dodge")
    ggsave(p, file = "Ortho_ratio.pdf", width = 9, height = 4)
'
```

### Picture all domains

```bash
cd ~/data/symbio/DOMAIN
mkdir -p PICTURE
mkdir -p SEQ

ls RLK/*.RLK.lst |
    sed 's/\.RLK\.lst$//' |
    parallel -j 12 -keep-order '
        echo "==> {/}"
        faops some -l 0 ../primary_transcripts/{/}.pep.fa \
            {}.RLK.lst SEQ/{/}.RLK.fa
        faops size SEQ/{/}.RLK.fa \
            > SEQ/{}.len.tsv
    '
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

```bash
perl ../scripts/tmhmm_result.pl -i TMD/actinidia_chinensis.tsv > actinidia_chinensis.RLK.lst

cat pfam/actinidia_chinensis.pfam.tsv | sed 1d | tsv-join -f actinidia_chinensis.RLK.lst -k 1 > actinidia_chinensis.RLK.tsv


```

unrecognized option '--gdwarf-5'

OG0001905
OG0002289
OG0005768
OG0001821
OG0010386

```bash
cd ~/data/symbio/RLK

cat ortho.lst |
    parallel -j 1 --keep-order '
        cat ../Orthogroups/all_pro_ortho.tsv |
            tsv-filter --str-eq 2:{} |
            cut -f 1 |
            tsv-join -f ../PROTEIN/all_pro_spe.name.tsv -k 1 -a 2 |
            tsv-select -f 2 |
            tsv-join -f <(cat ../DOMAIN/taxonomy.tsv | sed 1d) -k 1 -a 2 |
            tsv-select -f 2 |
            tsv-uniq |
            tr '\''\n'\'' '\'','\'' |
            sed '\''s/,$/\n/'\'' |
            awk -v OG={} '\''{print (OG"\t"$0)}'\'' \
            >> OG_cover/taxo_cover.tsv
    '
```


### Pick pkinase domains only

```bash
cd ~/data/symbio/DOMAIN
mkdir -p RLK

# all pkinase protein names
cat species.lst |
    parallel -j 12 -k '
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-filter --or --str-eq 4:Pkinase \
                --str-eq 4:Pkinase_fungal \
                --str-eq 4:PK_Tyr_Ser-Thr |
            cut -f 1 |
            tsv-uniq \
            > RLK/{}.RLK.lst
    '

# all RLK that truly contained the pkinase domain
cat species.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat UNIQ/{}.sort_uniq.tsv |
            tsv-join -f RLK/{}.RLK.lst -k 1 \
            > RLK/{}.RLK.tsv
    '

echo -e "species\tRLK_num" > all_RLK_species.tsv
cat species.lst |
    parallel -j 12 -k '
        echo "==> {}"
        cat RLK/{}.RLK.tsv |
            tsv-select -f 1 |
            tsv-uniq |
            wc -l |
            awk -v spe={} '\''{print (spe"\t"$0)}'\'' \
            >> all_RLK_species.tsv
    '

cat all_PK_species.tsv |
    tsv-join -H -f all_RLK_species.tsv -k species -a RLK_num \
    > all_count_species.tsv

perl ~/Scripts/fig_table/xlsx2csv.pl -f ../info/taxo_genom.xlsx |
    mlr --icsv --otsv cat |
    tsv-select -H -f SpeciesName,ComName |
    sed 1d |
    perl -nle 'print lc $_;' |
    sed '1ispecies\ttaxo' \
    > taxonomy.tsv

cat taxonomy.tsv |
    tsv-join -H -f all_count_species.tsv \
        -k species -a pot_PK,true_PK,RLK_num \
    > Species_result.tsv
```
