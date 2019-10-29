#!/bin/bash
############################################################################################
# Script to adjust transcript coordinates created by FUCHS. Do not use this for Salmon processing.
############################################################################################

# $1 -> ".exon_counts.bed" file
# $2 -> output directory

############################################################################################
# Usage
#######

usage() {
    cat <<EOF

Usage:
        $0 [FUCHS bed file] [output dir]

Options:
      - [FUCHS bed file]            Exon count FUCHS output file (BED)
      - [output dir]                Directory where FUCHS bed file with corrected coordinates should be saved (BED)

Important:
	NEVER use this script on dataset meant for Salmon processes. ONLY use it if the FUCHS output is what you're interested in.
	All coordinate adjustments are built into the Salmon processing scripts. Otherwise the coordinates get corrected twice.

Example:
	parallel -j1 --xapply ./fuchs_adjust_coordinates.sh ~/mouse_RNAseR_hippocampus/circtools/03_reconstruct/1746_{} ./ ::: D H

EOF
}

if [ ! $# = 2 ]; then
    usage
    exit 1
fi

############################################################################################
# Allocation
############

in_dir=$(dirname `realpath $1`)
if [[ $1 = *exon_counts.bed ]]; then
      name=`echo $(basename $1) | sed "s%\.exon\_counts\.bed$%%"`
else  name=`echo $(basename $1) | sed "s%\.bed$%%"`
fi
name_out=`echo $(basename $1) | sed "s%\.bed$%\.corrected\.bed%"`

out_dir=`realpath $2`


############################################################################################
# Process
#########

mkdir -pv $out_dir

cd $in_dir

# remove empty rows & remove comments
grep -v "^$" $(basename $1) | grep -v "^#" > $out_dir/$name_out.tmp

# add chr prefix if not present
if [[ `head $(basename $1) | grep "^chr" | wc -l` = 0 ]]; then
      sed -i 's/^/chr/' $out_dir/$name_out.tmp
fi
# remove bars from dataset, they generate problems later on
sed -i 's/|/;/g' $out_dir/$name_out.tmp

# the transcript coordinates of our bed files are out of sync
# shift the transcript ranges 'up' by 1 nucleotide
if [[ $1 = *exon_counts.bed ]]; then
  bioawk -c bed '{$start=$start-1; $thickstart=$thickstart-1; print}' $out_dir/$name_out.tmp | \
      awk '{
            # create arrays from blocksizes and blockstarts
            split($11,a,",")
            split($12,b,",")
            # get rid of 1 nt exons at end
            if(a[$10]==1){
                  delete a[$10]
                  delete b[$10]
                  ## adjust exon counts
                  $10=$10-1
                  ## adjust transcript end coordinates
                  $3=$2+b[$10]+a[$10]+1
                  $8=$3
            }
            # get rid of 1 nt exons at start & update exon counts if 1 nt exon adjacent to next exon
            if(a[1]==1 && b[2]==1){
                  a[1]=a[2]
                  delete a[2]
                  delete b[2]
                  $10=$10-1
            # get rid of 1 nt exons at start & update exon counts & shift transcript if 1 nt exon not adjacent to next exon
            } else if(a[1]==1){
                  a[1]=a[2]
                  delete a[2]
                  $10=$10-1
                  shift=b[2]
                  $2=$2+shift-1
                  $7=$7+shift-1
                  for(i in b){
                        b[i]=b[i]-shift+1
                  }
                  b[1]=0
                  delete b[2]
            }
            # adjust exon sizes and start coordinates
            a2=""; b2=""
            for(i in a){
              if(i==1){
                  a2=a[i]+2
              } else { a2=a2 ","; a2=a2 a[i]+2
              }
            }
            for(i in b){
              if(i==1){
                  b2=b[i]
              ## skip exon nr 2 if it was removed
              } else if(i==2 && length(a[2])==0){
                  continue
              } else {
                  b2=b2 ","; b2=b2 b[i]-1
              }
            }
            $11=a2; $12=b2
            print
            # tab separate list & concatenate a unique name (UID) & remove instances with empty exon coordinates
            }' | sed 's/ /\t/g' | bioawk -c bed '{$name=$name":"$11":"$12; print}' | grep -v '::' > $out_dir/$name_out

else
      awk '{
            $2=$2-1
            $7=$7-1
            # create arrays from blocksizes and blockstarts
            split($11,a,",")
            split($12,b,",")
            # adjust exon sizes and start coordinates
            a2=""; b2=""
            for(i in a){
              if(i==1){
                  a2=a[i]+1
              } else {
                  a2=a2 ","; a2=a2 a[i]
              }
            }
            for(i in b){
              if(i==1){
                  b2=b[i]
              } else {
                  b2=b2 ","; b2=b2 b[i]+1
              }
            }
            $11=a2; $12=b2
            print
            # tab separate list & concatenate a unique name (UID) & remove instances with empty exon coordinates
            }' $out_dir/$name_out.tmp | sed 's/ /\t/g' | bioawk -c bed '{$name=$name":"$11":"$12; print}' | grep -v '::' > $out_dir/$name_out

fi

rm $out_dir/$name_out.tmp

echo "$name_out created and saved in $out_dir"

