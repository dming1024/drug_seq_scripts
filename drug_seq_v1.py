

#drug seq data preprocess
import os
import sys


import pandas as pd
import json
import re
import matplotlib.pyplot as plt
import numpy as np
import gzip
import pysam


######################################################
##1. 读入fastq1，获取barcode、UMI
##2. 对数据进行质控与可视化
######################################################
def get_fastq1(barcode,fastq1):
    R1_reads={}
    with gzip.open(fastq1, 'rt') as f:
        for line in f:
            sequence_identifier = line.strip()  # 序列标识符，通常是 '@' 开头的行
            sequence_ID=re.sub(" [0-9]*\\/[0-9]*","",sequence_identifier)
            if sequence_identifier.startswith('@'):
                sequence = next(f).strip()  # 序列行
                barcode_id=sequence[:10]
                UMI=sequence[-10:]
                _ = next(f).strip()
                quality_scores = next(f).strip()  # 质量分数行
                R1_reads[sequence_ID]=barcode_id+"_"+UMI

    barcodes=get_barcode_library(barcode)
    rs={} #用于存储解析后的barcode_UMI，仅限有效的barcodes
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
    plt.title('N of transcripts in 360 wells')# 添加标题和标签
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.savefig("UMI_counts_distribution.png")# 显示直方图
    plt.close()


    #每个wells中，UMI>1的转录本个数
    barcode_umi_dupliates={}
    for k,v in rs.items():
        barcode_umi_dupliates[k]=sum(1 for value in v.values() if value > 1)/len(v)
    df=pd.DataFrame.from_dict(barcode_umi_dupliates, orient='index')
    df.columns=['Rate']
    df.to_csv("UMI_duplicates.csv")
    
    plt.hist(df[['Rate']], bins=30, edgecolor='black')  # bins参数控制直方图的条形数，edgecolor设置条形边缘颜色
    plt.title('UMI duplicate rates')# 添加标题和标签
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.savefig("UMI_duplicates.png")# 显示直方图
    plt.close()
    
    return(R1_reads)


######################################################
##获取有效barcode，减少不必要的检索
##
######################################################
def get_barcode_library(barcode):
    barcode_id=pd.read_csv(barcode)
    barcodes=list(barcode_id['barcode_ID'])
    return barcodes

######################################################
##将fastq2比对到参考基因组
##
######################################################
def align_to_reference(sample,fastq2):
    _cmd={}
    _cmd['STAR']= '''/OIU/sge_software/STAR/STAR-2.4.2a/bin/Linux_x86_64/STAR \
        --genomeDir /OIU/database_NGS/hg38/reference/STAR_242a_index \
        --sjdbGTFfile /OIU/database_NGS/hg38/GTF/gencode.v45.primary_assembly.annotation.gtf \
        --limitBAMsortRAM 40000000000 \
        --runThreadN 4 \
        --limitIObufferSize 500000000 \
        --outFilterType BySJout \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterMultimapNmax 20 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterIntronMotifs None \
        --outSJfilterReads All \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline ID:{prefix} SM:{prefix} PL:ILLUMINA \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimScoreMin 0 \
        --chimScoreDropMax 20 \
        --chimScoreSeparation 10 \
        --chimScoreJunctionNonGTAG -1 \
        --quantMode TranscriptomeSAM \
        --quantTranscriptomeBan IndelSoftclipSingleend \
        --outReadsUnmapped Fastx \
        --readFilesIn  {R2} \
        --readFilesCommand zcat \
        --outFileNamePrefix  {prefix}.'''.format(prefix=sample,R2=fastq2)
        
    _cmd['sort']='''/OIU/sge_software/samtools/currentVersion/samtools sort {prefix}.Aligned.toTranscriptome.out.bam > {prefix}.Aligned.toTranscriptome.out_sort.bam'''.format(prefix=sample)
    _cmd['index']='''/OIU/sge_software/samtools/currentVersion/samtools index {prefix}.Aligned.toTranscriptome.out_sort.bam'''.format(prefix=sample)
    os.system(_cmd['STAR'])
    os.system(_cmd['sort'])
    os.system(_cmd['index'])


######################################################
##解析bam文件，获取reads和genes
##
######################################################
def get_reads_genes_from_BAM(sample):
    reads_genes={}#包括： reads-> genes
    with pysam.AlignmentFile("{}.Aligned.toTranscriptome.out_sort.bam".format(sample), "r") as bamfile:
        for read in bamfile.fetch():
            if not read.is_unmapped and read.query_name.startswith('SRR'):
                read_records=read.to_dict()
                if read_records['name'] in reads_genes:
                    reads_genes[read_records['name']].append(read_records['ref_name'])
                else:
                    reads_genes[read_records['name']]=[read_records['ref_name']]   

    return reads_genes

######################################################
##联合R1,R2分析结果，将readsID 置换为 barcode_ID 和 UMI
##
######################################################
def exchange_readsID_with_barcodeUMI(barcode,R1_reads,reads_genes):
    tmp={}
    barcodes=get_barcode_library(barcode)
    for k,v in reads_genes.items():
        r1=re.sub("\\_.*","",R1_reads["@"+k])
        if r1 in barcodes:#即为有效barcode
            if R1_reads["@"+k] in tmp:
                #print(tmp[R1_reads["@"+k]])
                tmp[R1_reads["@"+k]].update({g:np.round(1/len(v),2) for g in v})
                #print(tmp[R1_reads["@"+k]])
                #i +=1
                #print(i)#难道只是PCR重复带来的UMI重复？
                #print("--------------------------")
            else:
                tmp[R1_reads["@"+k]]={g:np.round(1/len(v),2) for g in v}
    return(tmp)

######################################################
##定量，在不同barcode中，定量转录本的表达
##
######################################################
def counts_UMI_for_transcripts_in_each_barcodes(tmp):
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
    return barcode_expression_df


######################################################
##定量，在不同barcode中，定量基因水平的表达
##
######################################################



#计算信息
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python drugseq.py barcode fastq1 fastq2 sample")
    else:
        barcode=sys.argv[1]
        fq1=sys.argv[2]
        fq2=sys.argv[3]
        sample=sys.argv[4]
        #step1
        R1_reads=get_fastq1(barcode, fq1) #barcode,fastq1
        #step2
        #align_to_reference(sample, fq2) ##barcode,fastq2, 比对
        #step3
        reads_genes=get_reads_genes_from_BAM(sample)
        #step4
        tmp=exchange_readsID_with_barcodeUMI(barcode,R1_reads,reads_genes)
        count_matrix=counts_UMI_for_transcripts_in_each_barcodes(tmp)        
        count_matrix.to_csv(sample+"_count_matrix.csv")
        
