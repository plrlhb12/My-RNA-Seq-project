#!/bin/bash
#SBATCH -J tophat-cufflink-0816
#SBATCH -p cn-short
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -o tophat-cufflink-0816.out
#SBATCH -e tophat-cufflink-0816.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000cns

## HL60_RNAseq tophat-cufflink pipeline

####—————————————————step 0:set the index gtf and inputfile—————————————————————————————————————————————————————————————————————————————————————————————

BOWTIE2IDX=/lustre1/lch3000_pkuhpc/chenq/05_PublicData/UCSC/hg19/Sequence/Bowtie2Index/genome
GenomeSEQ=/lustre1/lch3000_pkuhpc/chenq/05_PublicData/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
GTF=/lustre1/lch3000_pkuhpc/chenq/05_PublicData/UCSC/hg19/Annotation/Genes/genes.gtf

INPUTDIR=/lustre1/lch3000_pkuhpc/chenq/01_data/pha/HL60-RNA-seq-170727/P101SC17051666-01zjy-88G/data_release/clean_data
OUTPUTDIR=/lustre1/lch3000_pkuhpc/chenq/02_results/pha_RNAseq_HL60_TophatCufflinks

# mkdir TophatOut_HL60-1_H5
# tophat2 -p 8 -o TophatOut_HL60-1_H5 -G  $GTF $BOWTIE2IDX  $INPUTDIR/RNA-HL60-1_H53LHCCXY_L5_1.clean.fq.gz $INPUTDIR/RNA-HL60-1_H53LHCCXY_L5_2.clean.fq.gz

# ####—————————————————step1: reads alignments by tophat2. ——————————————————————————————————————————————————————————————————————————————————————————————
# stime=`date +"%Y-%m-%d %H:%M:%S"`
# echo "[$stime] step1: reads alignments by tophat2"
    
# ## ls |cat |grep "_1.clean.fq.gz"|wc

# cd  $INPUTDIR
# for  i in $(ls |cat | grep "_1.clean.fq.gz")
# do
# 	cd $OUTPUTDIR
# 	mkdir "TophatOut_"${i:4:9}
# 	# format: tophat2 -p thread -o outputdir -G annotation.gtf/gff bowtie2idx_base read1 read2
# 	tophat2 -p 8 -o "TophatOut_"${i:4:9} -G $GTF $BOWTIE2IDX  $INPUTDIR/${i}  $INPUTDIR/${i/_1.clean.fq.gz/_2.clean.fq.gz}   
# 	# echo "${i:4:9} reads alignments completed"
# done

# stime=`date +"%Y-%m-%d %H:%M:%S"`
# echo "[$stime]  step1: reads alignments by tophat2 completed......."

tophat2 -p 8 -o TophatOut_HL60-4h-2 -G $GTF $BOWTIE2IDX  $INPUTDIR/RNA-HL60-4h-2_H53LHCCXY_L5_1.clean.fq.gz  $INPUTDIR/RNA-HL60-4h-2_H53LHCCXY_L5_2.clean.fq.gz
tophat2 -p 8 -o TophatOut_HL60-4h-3 -G $GTF $BOWTIE2IDX  $INPUTDIR/RNA-HL60-4h-3_H53LHCCXY_L3_1.clean.fq.gz  $INPUTDIR/RNA-HL60-4h-3_H53LHCCXY_L3_2.clean.fq.gz


cd $OUTPUTDIR
# ####————————————————— step2: transcripts assembly by cufflinks. ——————————————————————————————————————————————————————————————————————————————————————————————


# stime=`date +"%Y-%m-%d %H:%M:%S"`
# echo "[$stime] step2: transcripts assembly by cufflinks"
    
# for  i in $(ls |cat | grep "TophatOut")
# do
# 	mkdir "CufflinksOut_"${i:10:9}

# 	# format: cufflinks -p threads -o outputdir -L label_prefix(default:CUFF) hits.sam
# 	cufflinks -p 8 -o "CufflinksOut_"${i:10:9} -L ${i:10:9} $OUTPUTDIR/TophatOut_"${i:10:9}"/accepted_hits.bam 

# 	# echo "${i:10:9} transcripts assembly completed"
# done

# stime=`date +"%Y-%m-%d %H:%M:%S"`
# echo "[$stime]  step2: transcripts assembly by cufflinks completed......."

cufflinks -p 8 -o CufflinksOut_HL60-4h-2 -L HL60-4h-2 $OUTPUTDIR/TophatOut_HL60-4h-2/accepted_hits.bam 
cufflinks -p 8 -o CufflinksOut_HL60-4h-3 -L HL60-4h-3 $OUTPUTDIR/TophatOut_HL60-4h-2/accepted_hits.bam 


cd $OUTPUTDIR
####————————————————— step3: merge assemblies from multiple RNA-seq libraries into a master transcriptome by cuffmerge. ——————————————————————————————————————————————————————————————————————————————————————————————
stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step3: merge assemblies from multiple RNA-seq libraries into a master transcriptome by cuffmerge"

# ls  CufflinksOut_*/*.gtf > ./mergelist.txt  produce mergelist.txt /////////good  idea/////////////
ls  $OUTPUTDIR/CufflinksOut_*/*cripts.gtf > ./assemblies_gtf_list.txt
cuffmerge -g $GTF -o CuffmergeOut -s $GenomeSEQ -p 8 assemblies_gtf_list.txt 

stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step3: merge assemblies from multiple RNA-seq libraries into a master transcriptome by cuffmerge  completed ......."

cd $OUTPUTDIR
####—————————————————  step4: differential expression analysis by cuffdiff—————————————————————————————————————————————————————————————————————————————————————————————
stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step4: differential expression analysis by cuffdiff "
mkdir Cuffdiff_Out

cuffdiff -o Cuffdiff_Out -b $GENOMESEQ -p 8 -L HL60_ctrl,HL60_4h,HL60_2d,HL60_4d -u CuffmergeOut/merged.gtf \
	 ./TophatOut_HL60-1*/accepted_hits.bam,./TophatOut_HL60-2*/accepted_hits.bam,./TophatOut_HL60-3*/accepted_hits.bam \
	 ./TophatOut_HL60-4h-1/accepted_hits.bam,./TophatOut_HL60-4h-2/accepted_hits.bam,./TophatOut_HL60-4h-3/accepted_hits.bam \
	 ./TophatOut_HL60-2d-1/accepted_hits.bam,./TophatOut_HL60-2d-2/accepted_hits.bam,./TophatOut_HL60-2d-3/accepted_hits.bam \
	 ./TophatOut_HL60-4d-1/accepted_hits.bam,./TophatOut_HL60-4d-2/accepted_hits.bam,./TophatOut_HL60-4d-3/accepted_hits.bam


stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step4: differential expression analysis by cuffdiff  completed......."

####—————————————————  step5: RNAseq_Visualization—————————————————————————————————————————————————————————————————————————————————————————————
stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step5: RNAseq_Visualization "

Rscript /lustre1/lch3000_pkuhpc/chenq/03_code/rnaseq/05_RNAseq_Visualization.R

stime=`date +"%Y-%m-%d %H:%M:%S"`
echo "[$stime] step5: RNAseq_Visualization completed......."
