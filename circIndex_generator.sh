#!/bin/bash


# $1 -> ".exon_counts.bed" file
# $2 -> output directory
# $3 -> raw read length (used for determining range of backsplice junction overlap)

############################################################################################
# Usage
#######

usage() {
    cat <<EOF

Usage:
        $0 ['.exon_counts.bed' file] [output dir] [raw read length]

Options:
	- ['.exon_counts.bed' file]	Exon count FUCHS output file or merged FUCHS output files (BED);
					currently needs this suffix!
	- [output dir]			Directory where circRNA reference (.circIndex.fa) output should be saved (FASTA)
	- [raw read length]		Optional input of raw data read lengths (found in .fastq.gz file)
					to determine the range of backsplice junction overlap; default is 125 bps

Examples:
	./circIndex_generator.sh ../circtools/03_reconstruct/1746_D_test/ ./mouse_RNAseR_hippocampus/
	parallel -j1 --xapply ./circIndex_generator.sh ../circtools/03_reconstruct/1746_{} ./mouse_RNAseR_hippocampus 250 ::: D H

EOF
}

if [ ! $# = 2 ] && [ ! $# = 3 ]; then
    usage
    exit 1
fi


############################################################################################
# Allocation
############

gene_reference_mouse="/reference/Mus_musculus/GENCODE/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"

in_dir=$(dirname `realpath $1`)
if [[ $1 = *exon_counts.bed ]]; then
	name=`echo $(basename $1) | sed "s%\.exon\_counts\.bed$%%"`
else 	name=`echo $(basename $1) | sed "s%\.bed$%%"`
fi
name_UID=`echo $(basename $1) | sed "s%\.bed$%\_UID\.bed%"`

out_dir=`realpath $2`

if [ -z $3 ]; then
	start_len=`echo $(( 125 / 2 ))`
	end_len=`echo $(( 125 - $start_len ))`
else
	start_len=`echo $(( $3 / 2 ))`
	end_len=`echo $(( $3 - $start_len ))`
fi


############################################################################################
# Process
#########

mkdir -pv $out_dir

cd $in_dir

# remove empty rows & remove comments & concatenate a unique name (UID) & remove instances with empty exon coordinates
grep -v "^$" $(basename $1) | grep -v "^#" | bioawk -c bed '{$name=$name":"$11":"$12; print}' | grep -v '::' > $out_dir/$name_UID.tmp

# add chr prefix if not present
if [[ `head $(basename $1) | grep "^chr" | wc -l` = 0 ]]; then
	sed -i 's/^/chr/' $out_dir/$name_UID.tmp
fi
# remove bars from dataset, they generate problems later on
sed -i 's/|/;/g' $out_dir/$name_UID.tmp

# the transcript coordinates of our bed files are out of sync
# shift the transcript ranges 'up' by 1 nucleotide & adjust the blocksize of the transcripts
#bioawk -c bed '{$start=$start-1; $thickstart=$thickstart-1; print}' $out_dir/$name_UID.tmp | \
#	awk '{
#		split($11,a,",")
#		a[length(a)]=a[length(a)]+2
#		a[
#		a2=""
#		for(i in a){
#		  if(a2!=""){ a2=a2 ","}
#		  a2=a2 a[i]
#		}
#		$11=a2
#		print
#		}' | sed 's/ /\t/g' > $out_dir/$name_UID

bioawk -c bed '{$start=$start-1; $thickstart=$thickstart-1; print}' $out_dir/$name_UID.tmp | \
	awk '{
            split($11,a,",")
		split($12,b,",")
#            a[2]=a[2]-1
            a2=""; b2=""
		if(a[1]==1){
		  a[1]=a[2]
		  delete a[2]
		  delete b[2]
		  $10=$10-1
		}
            for(i in a){
              if(i==1){
			a2=a[i]+2
		  } else { a2=a2 ","; a2=a2 a[i]+2
		  }
            }
		for(i in b){
		  if(i==1){
			b2=b[i]
		  } else { b2=b2 ","; b2=b2 b[i]-1
		  }
		}
            $11=a2; $12=b2
            print
            }' | sed 's/ /\t/g' > $out_dir/$name_UID

#rm $out_dir/$name_UID.tmp

# create fasta reference file from the .exon_counts_UID.bed file

## The fasta reference created with bedtools getfasta is modified with bioawk to make the start and end of each sequence overlap
## This way we'll be able to map any read over the backsplice junctions of the reference sequences.

cd $out_dir

echo "start of bedtools getfasta"

bedtools getfasta -fi $gene_reference_mouse \
		-bed $name_UID \
		-name \
		-split \
	-fo 1746_getfasta.tmp \
		-s #| \
#		bioawk -c fastx -v s_len=$start_len -v e_len=$end_len \
#				'{start=substr($seq, 1, s_len)
#                		end=substr($seq, length($seq)-e_len+1, length($seq))
#                		$seq=end $seq start
#               		print(">" $name "\n" $seq)}' > $name.circIndex.fa

echo "bedtools getfasta & sequence manipulation finished"

#rm $out_dir/$name_UID

echo "$name.circIndex.fa created and saved in $out_dir"
