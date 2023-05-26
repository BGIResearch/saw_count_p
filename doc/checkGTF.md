# 检查注释文件格式

## 检查格式

## 修复格式

无论gtf/gff,均通过gene id/transcript id来对数据进行分组,也要按组输出

如gff缺失Parent,按照顺序将前面最近合法记录作为Parent

### 已知错误格式类型

低级错误

* 第七列,正负链的 "-" "_" 混肴

属性缺失

* gtf 缺失gene_id/gene_name/transcript_id/transcript_name
* gtf 缺失transcript的行
* gtf 缺失gene的行
* gff mRNA缺失Parent

### 解决方案

低级错误,如正负链问题,对每一行进行检查并替换

gtf缺失属性,对每一行,使用已存的id/name来补全缺失的id/name,如果两个都缺失则无法弥补

gtf缺失行,需要读取exon的数据,解析gene/transcript id/name等属性,当一组数据结束,则在前面添加缺失的 gene/transcript行

gff缺失属性,ID/Name来补齐,如果确实的Parent,则需要记录上一条记录的ID,作为当前缺失的Parent