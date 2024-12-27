

# Drug Seq Analysis test

参考文章： [DRUG-seq Provides Unbiased Biological Activity Readouts for Neuroscience Drug Discovery](https://pubs.acs.org/doi/10.1021/acschembio.1c00920)

+ 数据集： GSE176150

## 数据集下载


```
#使用ftp，匿名登录进行下载：可以借助winscp进行登录
ftp ftp.sra.ebi.ac.uk
Name: anonymous
Password:
ftp> cd vol1/fastq/ERR164/ERR164407
ftp> get ERR164407.fastq.gz


#ascp下载
#download方法：https://www.ebi.ac.uk/ena/browser/downloading-data

era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR147/001/SRR14730301/SRR14730301_1.fastq.gz
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz

http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz	
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz	

#测试数据下载
cli/bin/ascp.exe -T -i asperaweb_id_dsa.openssh --mode=recv --host=fasp.sra.ebi.ac.uk -P 33001 --user=era-fasp vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz .
cli/bin/ascp.exe  -T -i asperaweb_id_dsa.openssh --mode=recv --host=fasp.sra.ebi.ac.uk -P 33001 --user=era-fasp vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz .
cli/bin/ascp.exe -P33001 -O33001 -QT -l 200M -i asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR991/009/SRR9916599/SRR9916599_1.fastq.gz .
cli/bin/ascp.exe -P33001 -O33001 -QT -l 200M -i asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz .

cli/bin/ascp.exe  -P33001 -O33001 -QT -l 200M  -i asperaweb_id_dsa.openssh --mode=recv --host=fasp.sra.ebi.ac.uk  --user=era-fasp vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz .
$ascp -P33001 -O33001 -QT -l 200M  -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh --mode=recv --host=fasp.sra.ebi.ac.uk  --user=era-fasp vol1/fastq/SRR147/001/SRR14730301/SRR14730301_2.fastq.gz .

ls -lhtr
 8.7G Jun 12  2021 SRR14730301_1.fastq.gz
 21G Apr  8  2024 SRR14730301_2.fastq.gz

```


## 数据预处理

参考 `drug_seq_v2.py` ，获取样本的基因表达矩阵，分析结果在：

 + test/



## 下游分析

参考 `drug_seq.ipynb` 文档，主要包括：

+ 相关性
+ 聚类分析
+ 差异比较
+ 功能分析

结果包括在 `result_202404.pptx`


## GSM5357044验证

从GEO下载数据集，进行单次trial的验证

