#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>, Tomas Germade <tgermade>
# @Date:   Friday, June 14, 2019 17:55
# @Email:  tobias.jakobi@med.uni-heidelberg.de, tgermade@student-net.ethz.ch
# @Project: MSc Semester Project, ETH Zuerich, Laboratory of Systems Neuroscience
# @Last modified by:   tgermade
# @Last modified time: Friday, June 14, 2019 17:55
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=250G
#SBATCH -J "circtools alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de


# check number of arguments to display help
case $# in 5|7) ;;
	*)
		echo "Usage paired: $0 [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read 1 file] [Read 2 file] [Read 1 marker, e.g. R1]"
		echo "Usage unpaired: $0 [STAR index] [target dir] [gene annotation GTF file] [thread number] [Read]"
		exit
esac


# $1 -> Genome index
# $2 -> target directory
# $3 -> gene annotation file
# $4 -> number of threads to use
# $5 -> read / read 1
# $6 -> read 2
# $7 -> read 1 marker

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
if [ $# = 5 ]; then
      format="unpaired" &&
      target=`echo $5 | sed "s%\.fastq\.gz$%%"` &&
      # create the target directory, STAR will not do that for us
      mkdir -pv $2/$target

elif [ $# = 7 ]; then
	format="paired" &&
	target=`expr ${5/$7/} : '\(.*\)\..*\.'` &&
	#target=`echo $5 | sed "s%\.fastq\.gz$%%;s%$7%%"` &&
	# create the target directory, STAR will not do that for us
	mkdir -pv $2/$target &&
	mkdir -pv $2/$target/mate1/ &&
	mkdir -pv $2/$target/mate2/
fi

OLD_PATH=`pwd`

# main mapping part

# Run on unpaired data
if [ $format = "unpaired" ]; then
	STAR  --runThreadN $4\
	      --genomeDir $1\
	      --genomeLoad NoSharedMemory\
	      --readFilesIn $5\
	      --readFilesCommand zcat\
	      --outFileNamePrefix $2/$target/\
	      --outReadsUnmapped Fastx\
	      --outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
		--outSJfilterOverhangMin 15   15   15   15\
        	--outFilterMultimapNmax 20\
        	--outFilterScoreMin 1\
        	--outFilterMatchNminOverLread 0.7\
        	--outFilterMismatchNmax 999\
        	--outFilterMismatchNoverLmax 0.05\
        	--alignIntronMin 20\
        	--alignIntronMax 1000000\
        	--alignMatesGapMax 1000000\
        	--alignSJoverhangMin 15\
     		--alignSJDBoverhangMin 10\
        	--alignSoftClipAtReferenceEnds No\
        	--chimSegmentMin 15\
        	--chimScoreMin 15\
        	--chimScoreSeparation 10\
        	--chimJunctionOverhangMin 15\
        	--sjdbGTFfile $3\
        	--quantMode GeneCounts\
        	--twopassMode Basic\
        	--chimOutType Junctions SeparateSAMold &&

	cd $2/$target &&

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam &&

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam &&

	grep "^@" Aligned.out.sam > header.txt &&

	rm -f Aligned.out.sam &&
	rm -f Chimeric.out.sam &&

	rm -f -r _STARgenome &&
	rm -f -r _STARpass1 &&

	samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin &&
	samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp ;
	mv Aligned.noS.tmp Aligned.noS.bam &&
	samtools index Aligned.noS.bam &&

	samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin &&
	samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp &&
	mv Chimeric.noS.tmp Chimeric.noS.bam &&
	samtools index Chimeric.noS.bam &&

	rm -f Aligned.noS.sam &&
	rm -f Chimeric.noS.sam &&

	cd $OLD_PATH

# Run on paired data
elif [ $format = "paired" ]; then
	STAR	--runThreadN $4\
		--genomeDir $1\
		--genomeLoad NoSharedMemory\
		--readFilesIn $5 $6\
		--readFilesCommand zcat\
		--outFileNamePrefix $2/$target/\
		--outReadsUnmapped Fastx\
		--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
		--outSJfilterOverhangMin 15   15   15   15\
		--outFilterMultimapNmax 20\
		--outFilterScoreMin 1\
		--outFilterMatchNminOverLread 0.7\
		--outFilterMismatchNmax 999\
		--outFilterMismatchNoverLmax 0.05\
		--alignIntronMin 20\
		--alignIntronMax 1000000\
		--alignMatesGapMax 1000000\
		--alignSJoverhangMin 15\
		--alignSJDBoverhangMin 10\
		--alignSoftClipAtReferenceEnds No\
		--chimSegmentMin 15\
		--chimScoreMin 15\
		--chimScoreSeparation 10\
		--chimJunctionOverhangMin 15\
		--sjdbGTFfile $3\
		--quantMode GeneCounts\
		--twopassMode Basic\
		--chimOutType Junctions SeparateSAMold

	cd $2/$target

	gzip Unmapped.out.mate1
	gzip Unmapped.out.mate2

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

	grep "^@" Aligned.out.sam > header.txt

	rm -f Aligned.out.sam
	rm -f Chimeric.out.sam

	rm -f -r _STARgenome
	rm -f -r _STARpass1

	samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
	samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
	mv Aligned.noS.tmp Aligned.noS.bam
	samtools index Aligned.noS.bam

	samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
	samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
	mv Chimeric.noS.tmp Chimeric.noS.bam
	samtools index Chimeric.noS.bam

	rm -f Aligned.noS.sam
	rm -f Chimeric.noS.sam

	cd $OLD_PATH

	## done with main mapping

	## mapping mate1 now


	STAR	--runThreadN $4\
		--genomeDir $1\
		--genomeLoad NoSharedMemory\
		--readFilesIn $2/$target/Unmapped.out.mate1.gz\
		--readFilesCommand zcat\
		--outFileNamePrefix $2/$target/mate1/ \
		--outReadsUnmapped Fastx\
		--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
		--outSJfilterOverhangMin 15   15   15   15\
		--outFilterMultimapNmax 20\
		--outFilterScoreMin 1\
		--outFilterMatchNminOverLread 0.7\
		--outFilterMismatchNmax 999\
		--outFilterMismatchNoverLmax 0.05\
		--alignIntronMin 20\
		--alignIntronMax 1000000\
		--alignMatesGapMax 1000000\
		--alignSJoverhangMin 15\
		--alignSJDBoverhangMin 10\
		--alignSoftClipAtReferenceEnds No\
		--chimSegmentMin 15\
		--chimScoreMin 15\
		--chimScoreSeparation 10\
		--chimJunctionOverhangMin 15\
		--sjdbGTFfile $3\
		--quantMode GeneCounts\
		--twopassMode Basic\
		--chimOutType Junctions SeparateSAMold

	cd $2/$target/mate1/

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

	grep "^@" Aligned.out.sam > header.txt

	rm -f Aligned.out.sam
	rm -f Chimeric.out.sam

	rm -f -r _STARgenome
	rm -f -r _STARpass1

	samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
	samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
	mv Aligned.noS.tmp Aligned.noS.bam
	samtools index Aligned.noS.bam

	samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
	samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
	mv Chimeric.noS.tmp Chimeric.noS.bam
	samtools index Chimeric.noS.bam

	rm -f Aligned.noS.sam
	rm -f Chimeric.noS.sam

	cd $OLD_PATH

	## done with mate1 mapping

	## mapping mate2 now

	STAR	--runThreadN $4\
		--genomeDir $1\
		--genomeLoad NoSharedMemory\
		--readFilesIn $2/$target/Unmapped.out.mate2.gz\
		--readFilesCommand zcat\
		--outFileNamePrefix $2/$target/mate2/ \
		--outReadsUnmapped Fastx\
		--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
		--outSJfilterOverhangMin 15   15   15   15\
		--outFilterMultimapNmax 20\
		--outFilterScoreMin 1\
		--outFilterMatchNminOverLread 0.7\
		--outFilterMismatchNmax 999\
		--outFilterMismatchNoverLmax 0.05\
		--alignIntronMin 20\
		--alignIntronMax 1000000\
		--alignMatesGapMax 1000000\
		--alignSJoverhangMin 15\
		--alignSJDBoverhangMin 10\
		--alignSoftClipAtReferenceEnds No\
		--chimSegmentMin 15\
		--chimScoreMin 15\
		--chimScoreSeparation 10\
		--chimJunctionOverhangMin 15\
		--sjdbGTFfile $3\
		--quantMode GeneCounts\
		--twopassMode Basic\
		--chimOutType Junctions SeparateSAMold

	cd $2/$target/mate2/

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

	awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Chimeric.out.sam > Chimeric.noS.sam

	grep "^@" Aligned.out.sam > header.txt

	rm -f Aligned.out.sam
	rm -f Chimeric.out.sam

	rm -f -r _STARgenome
	rm -f -r _STARpass1

	samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Aligned.noS.bam /dev/stdin
	samtools reheader header.txt Aligned.noS.bam > Aligned.noS.tmp
	mv Aligned.noS.tmp Aligned.noS.bam
	samtools index Aligned.noS.bam

	samtools view -bS Chimeric.noS.sam | samtools sort -@ 10 -m 2G -T tempo -o Chimeric.noS.bam /dev/stdin
	samtools reheader header.txt Chimeric.noS.bam > Chimeric.noS.tmp
	mv Chimeric.noS.tmp Chimeric.noS.bam
	samtools index Chimeric.noS.bam

	rm -f Aligned.noS.sam
	rm -f Chimeric.noS.sam
fi
