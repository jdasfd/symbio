# 植物分类

- Selaginella为卷柏属，Lycopodium为石松属，之前都归在Pteridophyta（蕨类植物门）之下。后来将蕨类植物中的卷柏分出，作为石松类植物。目前石松类植物只有两个基因组，均为卷柏属下的植物，其中一个并不含有注释文件。目前Lycopodiopsida（石松纲）下有3个部分：Isoetales/quillworts（水韭目），Lycopodiales（石松目），Selaginellales/spike mosses（卷柏目）。
- 蕨类（Polypodiidae/Fern）：目前主要包括3个已测序的目
  - 桫椤目（Cyatheales）
  - 水龙骨目（Polypodiales），两个测序的物种都在凤尾蕨科（Pteridaceae）
  - 槐叶萍目（Salviniales），两个测序的物种都在槐叶萍科（Salviniaceae）
    - Adiantum capillus-veneris（铁线蕨）目前刚刚被测序完成并组装，且提供注释文件，目前尚未被数据库收录。
- 裸子植物（Gymnosperm）：目前主要包括4个主要的lineages
  - 苏铁（Cycads）
  - 银杏（Ginkgo）
  - 松柏类（Conifers）
  - 买麻藤（Gnetophyte）

```bash
cd ~/data/symbio/info

SPECIES=$(
    nwr member -r species Azolla |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v " cf." |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    wc -l
```

## 一些中文名称

- Fern - 蕨类：
  - Alsophila spinulosa - 桫椤
  - Ceratopteris richardii - 美洲水蕨
  - Adiantum capillus - 铁线蕨
  - Azolla filiculoides - 细叶满江红
  - Salvinia cucullata - 勺叶槐叶萍
- Gymnosperm - 裸子植物
  - Cycas panzhihuaensis - 攀枝花苏铁 
  - Conifers II：
    - Taxus wallichiana - 红豆杉
    - Taxus wallichiana var. chinensis (old: Taxus chinensis) - 红豆杉
    - Sequoiadendron giganteum - 巨杉（giant sequoia）
  - Conifers I：
    - Pinus taeda - 火炬松（loblolly pine）
    - Pinus tabuliformis - 油松（southern Chinese pine）
    - Picea abies - 欧洲云杉（Norway spruce）
    - Picea glauca - 白云杉（White spruce）
    - Abies alba - 银冷杉（silver fir）
  - Gnetidae（买麻藤纲）：
    - Welwitschia mirabilis - 百岁兰
    - Gnetum montanum - 买麻藤
- Angiosperm - 被子植物
  - 
