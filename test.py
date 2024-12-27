

#统计360+ wells中UMI的分布
import pandas as pd
import json
import re
import matplotlib.pyplot as plt
import numpy as np
import gzip




#R1 reads
# 打开压缩的FASTQ文件
with gzip.open('SRR14730301_1.fastq.gz', 'rt') as f:
    for line in f:
        # FASTQ文件格式是每4行为一个循环：序列名称，序列，加，质量分数
        # 通常，我们只关心序列和质量分数
        sequence_identifier = line.strip()  # 序列标识符，通常是 '@' 开头的行
        sequence_ID=re.sub(" [0-9]*\\/[0-9]*","",sequence_identifier)
        if sequence_identifier.startswith('@'):
            sequence = next(f).strip()  # 序列行
            barcode_id=sequence[:10]
            UMI=sequence[-10:]
            # 跳过加号行，因为通常我们不使用它
            _ = next(f).strip()
            quality_scores = next(f).strip()  # 质量分数行
            # 现在你可以处理序列和质量分数
            #print(sequence_ID, barcode_id, UMI)
            R1_reads[sequence_ID]=barcode_id+"_"+UMI


with open("data.json") as f:
    results=f.read()
    R1_reads=json.loads(results)
    
    
barcode_id=pd.read_csv("barcode_id.csv")
barcodes=list(barcode_id['barcode_ID'])


import re
rs={} #用于存储解析后的barcode_UMI，仅限有效的barcodes
all_barcodes=[] #存储所有的barcodes


for k,v in R1_reads.items():
    r1=re.sub("\\_.*","",v)
    if r1 in barcodes:
        r2=re.sub(".*\\_(.*)","\\1",v)
        if r1 in rs:
            if r2 in rs[r1]:
                rs[r1][r2] +=1
            else:
                rs[r1][r2]=1
        else:
            rs[r1]={r2:1}

#每个barcode下的UMI，转录本数目
barcode_umi_library_size={}
for k,v in rs.items():
    barcode_umi_library_size[k]=len(v)
    
df=pd.DataFrame.from_dict(barcode_umi_library_size, orient='index')
df.columns=['N']
df.to_csv("UMI_counts_in_Barcodes.csv")

plt.hist(df[['N']], bins=30, edgecolor='black')  # bins参数控制直方图的条形数，edgecolor设置条形边缘颜色
# 添加标题和标签
plt.title('N of transcripts in 360 wells')
plt.xlabel('Value')
plt.ylabel('Frequency')
# 显示直方图
plt.savefig("UMI_counts_distribution.png")


#每个wells中，UMI>2的转录本个数
barcode_umi_dupliates={}
for k,v in rs.items():
    barcode_umi_dupliates[k]=sum(1 for value in v.values() if value > 1)/len(v)

df=pd.DataFrame.from_dict(barcode_umi_dupliates, orient='index')
df.columns=['Rate']
df.to_csv("UMI_duplicates.csv")
plt.hist(df[['Rate']], bins=30, edgecolor='black')  # bins参数控制直方图的条形数，edgecolor设置条形边缘颜色
# 添加标题和标签
plt.title('UMI duplicate rates')
plt.xlabel('Value')
plt.ylabel('Frequency')
# 显示直方图
plt.savefig("UMI_duplicates.png")



######################################################
##bam 比对
##
##
##
######################################################




######################################################
##开始对bam文件进行计数，获取各个barcode中基因表达数据
##
##
##
######################################################

reads_genes={}#包括： reads-> genes
with pysam.AlignmentFile("test.Aligned.toTranscriptome.out_sort.bam", "r") as bamfile:
    for read in bamfile.fetch():
        if not read.is_unmapped and read.query_name.startswith('SRR'):
            read_records=read.to_dict()
            if read_records['name'] in reads_genes:
                reads_genes[read_records['name']].append(read_records['ref_name'])
            else:
                reads_genes[read_records['name']]=[read_records['ref_name']]                


tmp={}
for k,v in reads_genes.items():
    r1=re.sub("\\_.*","",R1_reads["@"+k])
    if r1 in barcodes:
        if R1_reads["@"+k] in tmp:
            tmp[R1_reads["@"+k]].update({g:np.round(1/len(v),2) for g in v})
        else:
            tmp[R1_reads["@"+k]]={g:np.round(1/len(v),2) for g in v}


barcode_expression={}
for k,v in tmp.items():
    r1=re.sub("\\_.*","",k)#获取barcode 序列
    for g,c in v.items():
        if r1 in barcode_expression:
            if g in barcode_expression[r1]:
                barcode_expression[r1][g] += c
            else:
                barcode_expression[r1][g]=c
        else:
            barcode_expression[r1] = {g:c}
            
barcode_expression_df=pd.DataFrame.from_dict(barcode_expression)
barcode_expression_df.to_csv("Genes_expression_in_Barcodes.csv")
