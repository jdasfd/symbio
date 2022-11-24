#!/usr/bin/env sh

usage () { echo "bash pro_extract.sh <.csv> <primary_path> <orthogroup_num>" 1>&2; exit; }
[ $# -lt 3 ] && usage

CSV=$1;
PRI_PATH=$2;
FILE_NAME=$3;

cat ${CSV} |
    perl -nla -F',' -e '
        shift @F; next if @F == 0;
        print join ("\n", @F);
        ' |
    sed 's/^"//g' |
    sed 's/"$//g' |
    sed 's/^\s//g' |
    sed 1d \
    > ${FILE_NAME}.tmp.lst

if [ -f ${FILE_NAME}.fa ]; then
    rm ${FILE_NAME}.fa
fi

for file in $(ls ${PRI_PATH})
do
    faops some ${PRI_PATH}/${file} \
    ${FILE_NAME}.tmp.lst stdout \
    >> ${FILE_NAME}.fa
done

name_num=$(cat ${FILE_NAME}.tmp.lst | wc -l)
pro_num=$(cat ${FILE_NAME}.fa | grep '^>' | wc -l)

if [[ $name_num == $pro_num ]]; then
    echo "==> Extraction complete"
    rm ${FILE_NAME}.tmp.lst
else
    echo "==> Please check name.tmp.lst"
    rm ${FILE_NAME}.fa
fi
