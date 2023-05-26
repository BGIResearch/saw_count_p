- [Bam2Gem程序处理流程](#bam2gem程序处理流程)
  - [程序版本](#程序版本)
  - [数据准备及格式](#数据准备及格式)
    - [bam格式要求](#bam格式要求)
    - [注释文件](#注释文件)
    - [其他](#其他)
  - [multi-map处理](#multi-map处理)
    - [前提条件](#前提条件)
    - [multi-map处理逻辑](#multi-map处理逻辑)
  - [读取bam](#读取bam)
  - [过滤](#过滤)
  - [注释及去重](#注释及去重)
    - [读取基因注释文件](#读取基因注释文件)
    - [注释reads](#注释reads)
    - [去重](#去重)
  - [umi校正及二次去重](#umi校正及二次去重)
    - [umi校正](#umi校正)
    - [二次去重](#二次去重)
  - [表达矩阵](#表达矩阵)
  - [输出](#输出)
    - [输出bam文件](#输出bam文件)
    - [统计](#统计)
    - [饱和度](#饱和度)
    - [基因表达矩阵](#基因表达矩阵)
      - [旧格式](#旧格式)
      - [新格式](#新格式)

# Bam2Gem程序处理流程

读取bam -> 比对质量过滤 -> 注释及去重 -> umi校正及二次去重 -> 输出结果

## 程序版本

Bam2Gem v2.3.2

## 数据准备及格式

输入文件包括bam文件和注释文件

### bam格式要求

* 按照比对位置排序
* 存在tag **Cx:i** 保存x坐标值,值区间为 [0,2^22),即最大坐标值小于4194304
* 存在tag **Cy:i** 保存y坐标值,值区间为 [0,2^22),即最大坐标值小于4194304
* 存在tag **UR:Z** 保存umi,将序列ACGT转换为[0,1,2,3]并使用16进制编码,字符串最大长度为16,即原始umi序列最长支持32bp
* 比对位置值的长度应小于 2^31

### 注释文件

基因注释文件gtf/gff的要求:

* 注释行以 # 开始
* 数据行共九列,以tab分隔;前八列有固定含义,第九列为以分号分隔的属性键值对列表
* gtf的属性键值对以**空格**分隔,gff以**等号**分隔
* 文件后缀名支持 *gtf/gtf.gz gff/gff.gz gff3/gff3.gz*
* start/end的最大值需小于2^31
* 最大有效基因数必须小于2^20,即1048576(程序目前将坐标x,y与gene转换为数值并编码成uint64,因此对坐标范围和基因个数有所限制)
* gtf文件不可以乱序,即同一个gene的transcript/exons需按顺序排列
* gff文件可以乱序,但仍需满足 *gene必须出现在对应mRNA之前,mRNA必须出现在对应的exon之前* 的规则

出现以下任一情况,读取注释文件时会忽略对应基因:

* gtf格式的属性没有gene_name gene_id transcript_name transcript_id(对gene只需要有gene_name gene_id)
* gff格式的属性没有ID Name Parent(对gene无需判断Parent)
* 同一个gene下的数据包含多个gene_id,日志打印 "Multiple gene IDs for gene xxx: id1, id2..."
* 同一个gene下的数据同时包含正反链,日志打印 "Strand disagreement for gene xxx"
* transcript或exon没有transcript_id,日志打印 "Record does not have transcriptID for gene xxx"
* 同一个gene有多条transcript的id相同,日志打印 "Transcript appears more than once for xxx"
* 存在exon的start > end,日志打印 "Exon has 0 or negative extent for xxx"
* 同一个transcript下的exons之间有overlap,日志打印 "Exons overlap for xxx"
* 一个gene没有任何transcript,日志打印 "No transcript for gene xxx"

ps: 一个contig下出现多条gene有相同的的gene_name,合并为一个gene

### 其他

* 最小线程数需不低于8. v2.0.0之后的版本对速度进行优化,对最小线程个数有一定要求

## multi-map处理

针对multi-map reads,以qname分组,每组挑选一条注释结果最可靠的矫正为uniq read

### 前提条件

以单个bam文件为单位处理multi-map reads,因此对每个bam文件有以下要求:

* 一个文件内,第一列qname相同则为同一条read的比对结果
* 程序针对BGI测序仪下机的数据有优化,要求qname格式形如FP200001234L1C003R02200689580
  其可拆分为:
  * FP200001234 芯片号
  * L1 通道号,长度两位,第一位固定为L,第二位为数字
  * C003 列号,第一位固定为C,其他为三位数字
  * R022 行号,第一位固定为R,其他为三位数字
  * 00689580 read序号,长度为8的数字
* 序列长度最大为255

其他格式qname程序同样支持,只是无法按规范进行编码,需要保存字符串,会导致内存消耗较大;程序优先按照BGI的qname格式进行处理,如不符合要求则按字符串处理

### multi-map处理逻辑

单read注释

* 对每条read进行注释,根据与基因组的overlap长度大于等于50%,计算得到其locus为exon/intron/intergenic的其中一个,并记录overlap长度
* 如遇到注释到多个基因的情况,按照优先级exon>intron选择基因
* 如果多个基因同属于exon或intron,选择exon或intron 的overlap长度最长的基因
* 如果overlap长度相同,则认为此read不注释到任何一个基因,归类为intergenic

选择一条read

* 将reads按qname分组
* 对每一组数据,按照优先级exon>intron选择read
* 如果多条reads同属于exon或intron,选择exon或intron 的overlap长度最长的read 
* 如果overlap长度相同,则无法挑选read矫正回uniq read,将所有reads的mapq设置为0,且flag设置为BAM_FSECONDARY
* 否则将overlap长度最长的read作为uniq read,将其mapq矫正为255,将剩余reads的mapq设置为0,且flag设置为BAM_FSECONDARY

## 读取bam

按contigs字典序查找每个contig数据,并按行读取

对于多个输入bam,同时打开所有文件,按比对位置顺序读取,相当于将多个有序列表合并为一个有序列表

## 过滤

按照比对质量过滤,丢掉小于比对质量阈值的reads,默认的比对质量阈值为10

## 注释及去重

### 读取基因注释文件

* 打开基因注释文件,跳过空行和注释行,以tab分割每一行
* 基础属性提取 chr/featureType/start/end/strand 分别对应列数 1/3/4/5/7
* 属性列提取
  * 针对gtf,gene_id/gene_name/transcript_id/transcript_name
  * 针对gff,读取ID/Name/Parent,并转换为gene_id/gene_name/transcript_id/transcript_name
* 将所有记录按照gene_name聚合并构建区间树,用于后续注释

### 注释reads

对每一条read,通过与注释文件的基因区间进行overlap的查找,计算得出注释上的gene_name/strand以及属于EXONIC/INTRONIC/INTERGENIC的哪一个,并对个数进行统计

* 解析cigar信息. 获取read全长以及每个align block的起始位置与长度
* 查找overlap的所有基因
* 确定最终选择哪个基因. 对每一个gene进行计算:
  * 对每一个align block,与当前gene的每一个transcript进行计算,根据与exon/非exon区域overlap的长度,得出多个exon_cnt/intron_cnt,取最大的exon_cnt(优先级为exon>intron)
  * 累加每一个align block的cnts,得出最优的exon_cnt/intron_cnt,如果exon_cnt大于等于reads长度的50%,则标记为EXONIC,否则,如果intron_cnt大于等于reads长度的50%,则标记为INTRONIC,否则标记为INTERGENIC
  * 从多个gene中选择最可靠的一个.
    * 首先得到注释结果最优的gene列表(优先级为EXONIC>INTRONIC>INTERGENIC)
    * 从这些gene中,取最大的overlap的gene作为注释结果
    * 如果多个gene的overlap长度相同,则随机选一个基因(具体挑选规则为选择gene的start与end较小的那个)
* 返回注释结果
  * 如果没有注释到任何一个gene,则直接返回,否则
  * 注释区域为EXONIC/INTRONIC/INTERGENIC三者中的一个
  * 如果为EXONIC/INTRONIC两者中的一个,则同时返回gene_name和gene_strand
  * 将注释结果更新到bam

### 去重

注释之后,按reads顺序构建 `{barcode_gene : {mid : cnt}}` 的嵌套map,即同一个DNB和gene组合下,含有的umi的种类及个数,
当遇到某一个组合下的umi出现次数大于1时,标记当前read为重复

## umi校正及二次去重

### umi校正

将由于测序错误等原因导致的错误umi根据汉明距离校正回正确umi,其过程:

* 数据准备. 形如 `{barcode_gene : {mid : cnt}}` 的嵌套map
* 校正.
  * 设置参数. 最小umi种类数阈值/容错个数阈值/umi长度,默认5/1/10
  * 校正. 对嵌套map的每组数据:
    1. 检查umi种类数,大于等于阈值则继续处理
    2. 按照umi的cnt个数降序排序,得到形如 `[(umi1,cnt1),(umi2,cnt2)...]` 的列表
    3. 逆序遍历排序后的列表,即从cnt最小的umi开始,分别与其他umi计算碱基错误数,如果满足容错个数阈值,则将当前umi校正为cnt较大的umi,并将当前umi的cnt转移给正确umi
    4. 得到校正后的形如 `{barcode_gene : {old_umi : new_umi}}` 的嵌套map
* 更新Bam. 对非重复且umi属于old_umi的reads增加 *UB* tag记录校正后的 new_umi

示例:

给定容错个数阈值为1,假设有umi及cnt数据排序后为:

```cpp
{
"AAA": 5,
"GGA": 4,
"AGA": 3,
"AAT": 2,
"GGG": 1,
"CCC": 1
}
```

校正过程:

* 将"CCC" 分别与 "AAA" - "GGG" 计算容错个数,均大于阈值1
* 将"GGG" 分别与 "AAA" - "AAT" 计算容错个数,发现当遇到"GGA"时容错个数为1,则更新两者的cnt值,并记录校正前后的对应关系
* 将"AAT" 分别与 "AAA" - "AGA" 计算容错个数,发现当遇到"AAA"时容错个数为1,则更新两者的cnt值,并记录校正前后的对应关系
* 将"AGA" 分别与 "AAA" - "GGA" 计算容错个数,发现当遇到"AAA"时容错个数为1,则更新两者的cnt值,并记录校正前后的对应关系
* 将"GGA" 与 "AAA" 计算容错个数,大于阈值1

更新原数据为:

```cpp
{
"AAA": 10,
"GGA": 5,
"AGA": 0,
"AAT": 0,
"GGG": 0,
"CCC": 1
}
```

保存校正前后的映射关系:

```cpp
{
"AGA": "AAA",
"AAT": "AAA",
"GGG": "GGA"
}
```

### 二次去重

将umi被校正的reads标记为重复,并增加tag标记正确的umi

## 表达矩阵

表达矩阵计算逻辑

* 挑选注释到EXON或INTRON的reads
* 过滤掉read方向与注释基因链方向相反的reads
* 按照坐标,gene,MID的顺序进行分组
* 统计每个坐标每个gene的唯一MID的个数,即为表达矩阵

## 输出

### 输出bam文件

根据输入bam进行修改,并输出新的bam

* 去掉未通过比对质量过滤的reads,如果设置参数--save_lq则保存reads,并将flag设置BAM_FQCFAIL
* 去掉重复read,如果设置参数--save_dup则保存reads,并将flag设置BAM_FDUP
* 添加注释信息
  * GE tag 表示基因名称
  * GS tag 表示基因链方向
  * XF tag 表示注释区域(EXONIC:0 INTRONIC:1 INTERGENIC:2)
* 添加umi校正信息,UB tag表示正确的umi

### 统计

输出统计文件,包含两部分

1. 总reads个数以及通过过滤和通过去重的reads个数

| TOTAL_READS | PASS_FILTER | UNIQUE_READS | FAIL_FILTER_RATE | DUPLICATION_RATE |
| ----------- | ----------- | ------------ | ---------------- | ---------------- |
| 9999404     | 6450176     | 4036230      | 35.4944          | 37.4245          |

2. 注释统计结果

| TOTAL_READS | MAP    | EXONIC | INTRONIC | INTERGENIC | TRANSCRIPTOME | ANTISENSE |
| ----------- | ------ | ------ | -------- | ---------- | ------------- | --------- |
| 949906      | 949906 | 855200 | 354      | 42033      | 850051        | 9833      |
| 100.0       | 100.0  | 90.0   | 0.0      | 9.9        | 89.5          | 1.0       |

统计指标描述:

| metrics       | 1.x                                      | 2.x                                                                                                           |
| ------------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| TOTAL_READS   | 进行注释的reads总数                      | 进行注释的reads总数                                                                                           |
| MAP           | 比对质量>=255                            | 比对质量>=255                                                                                                 |
| EXONIC        | 满足MAP且与exon的overlap大于等于50%      | 满足MAP且与<font color=red>链同向gene的</font>exon的overlap大于等于50%                                        |
| INTRONIC      | 满足MAP且不满足EXONIC且与intron有overlap | 满足MAP且不满足EXONIC且与<font color=red>链同向gene的</font>intron的overlap<font color=red>大于等于50%</font> |
| INTERGENIC    | 满足MAP且不满足EXONIC和INTRONIC          | 满足MAP且不满足EXONIC和INTRONIC                                                                               |
| TRANSCRIPTOME | 满足EXONIC且与链同向                     | 满足EXONIC<font color=red>或INTRONIC</font>且与链同向                                                         |
| ANTISENSE     | 注释到的gene与read链方向相反             | <font color=red>存在overlap gene</font>与read链方向相反                                                       |

2.x 相比 1.x 的更新:

* INTRONIC的判断从与intron有交集改为50%阈值
* 链方向的判断. EXONIC/INTRONIC的统计均加入链同向的要求;表达矩阵的数据亦满足链同向要求
* 注释文件gtf/gff文件要求更严格,如gff必须要由ID Name Parent等属性

### 饱和度

测序饱和度中间文件(仅当参数'--sat_file'生效时生成),共五列,以空格分隔

| coorY  | coorX | geneIndex | MIDIndex | readCount |
| ------ | ----- | --------- | -------- | --------- |
| 108565 | 84833 | 10539     | 1018895  | 1         |
| 114929 | 88696 | 10539     | 929853   | 4         |
  
coorY,coorX: dnb坐标
geneIndex: 基因的索引,编码为整型  
MIDIndex: MID序列索引,编码为整型 
readCount: MID的read个数

### 基因表达矩阵

#### 旧格式

旧格式为文本文件,文件共四列,以tab分隔,包含header
 
| geneID | x     | y      | MIDCount |
| ------ | ----- | ------ | -------- |
| Gm1992 | 94950 | 111098 | 1        |
| Gm1992 | 78404 | 113605 | 3        |

#### 新格式

新格式为二进制文件,gef格式,读取方式类似于hdf5文件

* dataset "/geneExp/bin1/expression" 包含一个类型为H5T_COMPOUND的数据列表,
其结构为{x, y, count},分别表示 x坐标,y坐标,MIDCount
* ATTRIBUTE "/geneExp/bin1/expression" 含有六项属性
  * maxExp dataset中最大的count数值
  * minX dataset中x坐标的offset,即坐标加上minX为原始坐标
  * minY dataset中y坐标的offset,即坐标加上minY为原始坐标
  * maxX dataset中原始坐标最大的x
  * maxY dataset中原始坐标最大的y
  * resolution 芯片分辨率
* dataset "/geneExp/bin1/gene" 包含一个类型为H5T_COMPOUND的数据列表,
其结构为{gene, offset, count},分别表示 基因名,此基因在expression中的起始位置,此基因的坐标个数

由两个dataset即可还原回旧格式的表达矩阵,即expression 按照顺序记录每个坐标的位置和MIDCount,gene 记录每个坐标对应的基因名称; 举例,假如有gene记录如 {"Xkr4",17,155}, 表示expression中第17个坐标对应基因为"Xkr4",且从第17个坐标开始的155个坐标均对应基因"Xkr4",注意offset是0-based
