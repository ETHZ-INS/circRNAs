#!/bin/bash


# @Author: Tomas Germade <tgermade>
# @Date:   Friday, November 01, 2019 17:40
# @Email:  tgermade@student-net.ethz.ch
# @Project: ETH Zuerich semester project Laboratory of Systems Neuroscience
# @Last modified by:   tgermade
# @Last modified time: Friday, November 01, 2019 17:40

# Improved version of 'circtools_master.sh'
# One wrapper to rule them all and in a bash script bind them.
############################################################################################
# Usage
#######

usage() {
    cat <<EOF

Usage:
      $0 [options]

Options:
	-m	which modules to run		can assume values: ['all' 'star' 'dcc' 'fuchs']
	-a	organism			can assume values: ['mouse' 'rat' 'rat_assembled' 'human']
	-p	form of read data 		can assume values: ['paired' 'unpaired']
	-i	input dir			specify location of input data (only relevant if not run with 'all' modules)
	-o	output dir        		specify location where output should be saved
	-r	read dir			specify location of raw reads (only relevant if run with 'all' modules)
	-t	thread number [optional]	if value is 'default', will use max. number of cores - 8 (min 2),
							otherwise input any number of cores you want to use (integer)
	-n	batch name [optional]		input unique name for the overall read collection;
							if no name is given, it will be automatically generated
							based on the supplied [read files dir]
							(only affects DCC module output)

Examples:
	$0 -m all \\
		-a mouse \\
		-p unpaired \\
		-r /path/to/reads/ \\
		-o /path/to/output/ \\
		-t 12 \\
		-n example_name

	$0 -m dcc \\
		-a human \\
		-p paired \\
		-i /path/to/input/ \\
		-o /path/to/output/ \\
		-t default \\
		-n example_name

Note:
	Option -r is specific to -m ['all']
	Option -d is specific to -m ['star' 'dcc' 'fuchs']

EOF
}

function error_exit {
	echo "ERROR: Please select a valid $1 for the $2 option"
	usage
	exit 1
}



############################################################################################
# Allocation
############

# create lists of possible options
declare -a module_opt=("all" "star" "dcc" "fuchs")
declare -a organism_opt=("mouse" "rat" "rat_assembled" "human")
declare -a mode_opt=("paired" "unpaired")

# allocation based on received options
while [ -n "$1" ]; do
	case "$1" in
		-m) 	module=$2
			if [[ ! ${module_opt[@]} =~ $module ]]; then
				error_exit "module" "-m"
			fi
			shift ;;
		-a) 	organism=$2
			if [[ ! ${organism_opt[@]} =~ $organism ]]; then
				error_exit "organism" "-a"
			fi
			shift ;;
		-p) 	mode=$2
			if [[ ! ${mode_opt[@]} =~ $mode ]]; then
				error_exit "read mode" "-p"
			fi
			shift ;;
		-i) 	in_dir=`realpath $2`
			shift ;;
		-o) 	out_dir=`realpath $2`
			shift ;;
		-r) 	read_dir=`realpath $2`
			shift ;;
		-t) 	procs=$2
			shift ;;
		-n)	batch_name=$2
			shift ;;
		*) 	echo "ERROR: option $1 not recognized"
			usage
			exit 1 ;;
	esac
	shift
done

# allocate annotations
if [ $organism = "mouse" ]; then
      star_index="/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes_STARIndex"
      genome_reference="/reference/Mus_musculus/GENCODE/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"
      gene_annotation="/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes.gtf"
      exon_annotation="/mnt/schratt/tgermade_test/mouse_GRCm38_p5.GENCODE.exons.bed"
elif [[ $organism =~ "rat" ]]; then
      star_index="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/shortRNA/starINDEX"
      genome_reference="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa"
      if [ $organism = "rat" ]; then
            gene_annotation="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf"
            exon_annotation="/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons.bed"
      elif [ $organism = "rat_assembled" ]; then
            gene_annotation="/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes_assembled.gtf"
            exon_annotation="/mnt/schratt/tgermade_test/rat_Rnor_6.0.Ensembl.exons_assembled.bed"
      fi
elif [ $organism = "human" ]; then
      star_index="/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Genes/genes_STARIndex"
      genome_reference="/reference/Homo_sapiens/Ensembl/GRCh38.p10/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.toplevel.fa"
      gene_annotation="/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Genes/features.gtf"
      exon_annotation="/mnt/schratt/tgermade_test/human_GRCh38_p10.Ensembl.exons.bed"
fi


#  test if -t is assigned an integer, else the script will be run with default thread numbers
if [ -z $procs ] || [[ ! $procs =~ ^[0-9]+$ ]]; then
	if [[ `nproc` < 10 ]]; then
		procs=2
	else 	procs=$(expr `nproc` - 8)
	fi
fi


# test if output is set
if [ -z $out_dir ]; then
	error_exit "output directory" "-o"
fi

# get full script path name
script_path=$(dirname `realpath $0`)
# generate temporary directory
mkdir -pv $out_dir/.tmp

# Modules -m all
if [ $module = "all" ]; then
	## testing
	if [ -z $read_dir ]; then
		error_exit "read directory" "-r"
	fi

	if [ ! -z $in_dir ]; then
		echo "ERROR: when all modules are run the -r option serves as input directory"
		usage
		exit 1
	fi
	## create directories
	mkdir -pv $out_dir
	mkdir -pv $out_dir/star
	mkdir -pv $out_dir/circtools/
	mkdir -pv $out_dir/circtools/01_detect
	mkdir -pv $out_dir/circtools/02_reconstruct

	## create batch name if not defined (output directory name for DCC module)
	if [ -z $batch_name ]; then
      	if [ `echo $(basename $read_dir)` = "raw" ]; then
            	batch_name=`echo $(dirname $read_dir) | rev | cut -d "/" -f1 | rev`
      	else
            	batch_name=`echo $(basename $read_dir)`
		fi
      fi

	## read file allocation
      if [ $mode = "paired" ]; then
		### define read markers for naming convention
            if [[ `echo $read_dir/*` = *_R1.fastq.gz* ]]; then
                  r1_marker=_R1 &&
                  r2_marker=_R2
      	else
      		echo "INPUT DEMAND: Input read 1 marker of raw read data files. Examples: '_R1', '_read1', etc." &&
            	read r1_marker &&
            	r2_marker=`echo $r1_marker | tr 1 2`
            fi
		### specify set of read files to analyze (prevent computation of sets already computed)
	      for file in $read_dir/*.fastq.gz; do
	            if [[ $file = *$r1_marker* ]]; then
     	                  echo $(basename $file) | sed "s%\.fastq\.gz$%%;s%$r1_marker%%"
     		      fi
		done > $out_dir/.tmp/names.tmp

	elif [ $mode = "unpaired" ]; then
		### specify set of read files to analyze (prevent computation of sets already computed)
		for file in $read_dir/*.fastq.gz; do
	            echo $(basename $file) | sed "s%\.fastq\.gz$%%"
	      done > $out_dir/.tmp/names.tmp
	fi
fi

# Modules -m star dcc fuchs
if [[ ${module_opt[@]:1} =~ $module ]]; then
	## testing
	if [ ! -z $read_dir ]; then
		echo "ERROR: -r is not available for the $module module; please select an input directory with the -i option"
		usage
		exit 1
	fi
	if [ -z $in_dir ]; then
		error_exit "input directory" "-i"
	fi

	## create directory
	mkdir -pv $out_dir

	## module-specific allocation
	if [ $module = "star" ]; then
		echo "under construction"

	elif [ $module = "dcc" ]; then
		## identify sample names
		for name in $in_dir/*.bam; do
            	echo $(basename $name) | sed "s%\.bam$%%"
      	done > $out_dir/.tmp/s_names.tmp

	elif [ $module = "fuchs" ]; then
		echo "under construction"
	fi
fi


echo "module: $module, in_dir: $in_dir, read_dir: $read_dir, out_dir: $out_dir, org: $organism, mode: $mode, threads: $procs"


############################################################################################
# Tools
#######



# ALIGNMENT


if [ $module = "all" ] || [ $module = "star" ]; then

      cd $read_dir &&

      # different wrapper options depending on paired / unpaired read data
      if [ $mode = "paired" ]; then
            align_option=$(echo "{1} {2} $r1_marker ::: *$r1_marker* ::: *$r2_marker*")
      elif [ $mode = "unpaired" ]; then
            align_option=$(echo "{} ::: *")
      fi

      # run STAR
      # PAIRED: slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [threads] [Read 1 file] [Read 2 file] [Read 1 marker]
      # UNPAIRED: slurm_circtools_detect_mapping [STAR index] [target dir] [gene annotation GTF file] [threads] [Read file]
      echo "Running STAR."
      parallel -j1 --xapply $script_path/slurm_circtools_detect_mapping.sh $star_index $out_dir/star/ $gene_annotation $procs $align_option
fi



# DETECTION


if [ $module = "all" ] || [ $module = "dcc" ]; then

      # setup links to /star folder
      echo "Setting up links."

	if [ $module = "all" ]; then

		dcc_dir="$out_dir/circtools/01_detect"
		star_dir="$out_dir/star/"

		declare -a from=("$star_dir/{}/*.bam" "$star_dir/{}/*.bam.bai" "$star_dir/{}/[Cc]himeric*")
		declare -a to=("$dcc_dir/{}.bam" "$dcc_dir/{}.bam.bai" "$dcc_dir/{}.Chimeric.out.junction")
		parallel_opt="$out_dir/.tmp/names.tmp"

		if [ $mode = "paired" ]; then

			from+=("$star_dir/{}/mate1/*.bam" "$star_dir/{}/mate2/*.bam" "$star_dir/{}/mate1/[Cc]himeric*" "$star_dir/{}/mate2/[Cc]himeric*")
			to+=("$dcc_dir/{}.mate1.bam" "$dcc_dir/{}.mate2.bam" "$dcc_dir/{}.mate1.Chimeric.out.junction" "$dcc_dir/{}.mate2.Chimeric.out.junction")
		fi

	elif [ $module = "dcc" ]; then

		dcc_dir="$out_dir"
		star_dir="$in_dir"

		declare -a from=("$star_dir/{}.bam" "$star_dir/{}.bam.bai" "$star_dir/{}.[Cc]himeric*")
		declare -a to=("$dcc_dir/{}.bam" "$dcc_dir/{}.bam.bai" "$dcc_dir/{}.bam.Chimeric.out.junction")
		parallel_opt="$out_dir/.tmp/s_names.tmp"

		if [ $mode = "paired" ]; then

			from+=("$star_dir/mate1/{}.bam" "$star_dir/mate2/{}.bam" "$star_dir/mate1/{}.[Cc]himeric" "$star_dir/mate2/{}.[Cc]himeric")
			to+=("$dcc_dir/{}.mate1.bam" "$dcc_dir/{}.mate2.bam" "$dcc_dir/{}.mate1.Chimeric.out.junction" "$dcc_dir/{}.mate2.Chimeric.out.junction")
		fi
	fi

	len=$(expr ${#from[@]} - 1)
	for i in $(seq 0 $len); do
		parallel --plus ln -s ${from[$i]} ${to[$i]} :::: $parallel_opt
	done

      cd $dcc_dir

      # create txt files
      if [ $mode = "paired" ]; then
            ls *bam | grep -v mate > bam_files.txt &&
            ls *[Cc]himeric* | grep  mate1 > mate1 &&
            ls *[Cc]himeric* | grep  mate2 > mate2 &&
            ls *[Cc]himeric* | grep -v mate > samplesheet

      elif [ $mode = "unpaired" ]; then
            ls *bam | grep -v mate > bam_files.txt &&
            ls *[Cc]himeric* | grep -v mate > samplesheet
      fi


	# check if data is unstranded, first-stranded or second-stranded
	cd $star_dir

	if [ -e dataset.tsv ]; then

		echo "Checking sample strandedness:"

            all_str=`tail -n+2 dataset.tsv | cut -f 9` &&
            first_str=`echo $all_str | cut -d' ' -f1` &&
            for str in $all_str; do
                  if [[ $str != $first_str ]]; then
                        echo "strandedness differs between samples! please subset samples according to strandedness and rerun"
                        exit 1
                  fi
            done

            if [ $first_str = "sense" ]; then
                  strandcall="firststrand"
            elif [ $first_str = "antisense" ]; then
                  strandcall="secondstrand"
            else
                  strandcall="unstranded"
            fi

	else
		echo "Checking sample strandedness (.85 cutoff):"

	      decideStrandedness () {
	            if (( $(echo "$1 >= 0.85" | bc -l) )); then
	                  echo "firststrand";
	            elif (( $(echo "$1 <= 0.15" | bc -l) )); then
	                  echo "secondstrand";
	            else
	                  echo "unstranded";
	            fi
	      }

		## for each sample, calculate ratio btw. first strand alignments and total alignments
	      sarr=()
	      parr=()

	      for dir in `cat $out_dir/.tmp/names.tmp`; do
	            cd $star_dir/$dir
	            p=`awk '{if(NR >= 5) { sumF+=$3; sumS+=$4 }} END{print sumF/(sumF+sumS)}' ReadsPerGene.out.tab;`
	            strandcall=`decideStrandedness $p`
	            parr+=($p)
	            sarr+=($strandcall)
	            echo "firststrand ratio: " $f $p "-->" $strandcall
	            cd ..
	      done

	      # check if all samples of batch share the same strandedness
	      if [ `printf '%s\n' "${sarr[@]}" | sort | uniq | wc -l` -gt 1 ]; then
	            echo "WARNING: not all libraries appear to have the same strand specificity!"
	            # check if (rough) median ratio passes a threshold
	            n=${#parr[@]}
	            strandcall=decideStrandedness `printf '%s\n' "${parr[@]}" | sort -n | tail -n $(( $n / 2 + 1 )) | head -1`
			echo "Falling back on $strandcall mode, but consider splitting the libraries."
	      else
	            # all calls are identical, take the first
	            strandcall=${sarr[0]}
	      fi
	fi

      echo "Running DCC on $strandcall mode."


	# allocate DCC options
	cd $dcc_dir

      if [ $strandcall = "unstranded" ]; then
            strand_option="-N"
      elif [ $strandcall = "secondstrand" ]; then
            strand_option="-ss"
      else strand_otpion=""
      fi

      if [ $mode = "paired" ]; then
            paired_option=$(echo "-mt1 @mate1 -mt2 @mate2 -Pi")
      else paired_option=""
      fi

	# run DCC
      source activate python27

      circtools detect  @samplesheet \
                  -D \
                  -an $gene_annotation \
                  -F \
                  -Nr 1 1 \
                  -fg \
                  -G \
                  -A $genome_reference \
                  -O ./output \
                  -T $procs \
                  -B @bam_files.txt \
                  $paired_option $strand_option

      conda deactivate
fi



# RECONSTRUCTION


if [ $module = "all" ] || [ $module = "fuchs" ]; then

      cd $out_dir/circtools/01_detect

############################################################################################

      circtools_reconstruct () {

            source activate python27

            # allocation
            sample_name=$1
            main_out=$2/
            bed_file=$3
            dcc_dir=$4
            dcc_out_dir=$5
            format=$6
            procs=$7
            tmp_folder=/tmp/global_tmp/

            main_bam=$dcc_dir/${sample_name}.bam &&
            main_junction=$dcc_dir/${sample_name}.Chimeric.out.junction &&

            mkdir -p $main_out/${sample_name}

            if [ $format = "paired" ]; then
                  mate1_bam=$dcc_dir/${sample_name}.mate1.bam &&
                  mate1_junction=$dcc_dir/${sample_name}.mate1.Chimeric.out.junction &&

                  mate2_bam=$dcc_dir/${sample_name}.mate2.bam &&
                  mate2_junction=$dcc_dir/${sample_name}.mate2.Chimeric.out.junction.fixed &&

                  merged_bam=$main_out/${sample_name}/${sample_name}_merged.bam
            fi

		# call FUCHS
            if [ $format = "paired" ]; then
                  # merge both mate BAM files into one new BAM file
                  samtools merge -l 9 -@ 8 $merged_bam $main_bam $mate1_bam $mate2_bam &&
                  # re-index the newly aggregated BAM file
                  samtools index $merged_bam &&
                  # run FUCHS
                  FUCHS -N $sample_name -D $dcc_out_dir/CircRNACount -B $merged_bam -A $bed_file -O $main_out/${sample_name} -F $mate1_junction -R $mate2_junction -J $
main_junction -T $tmp_folder -p ensembl -r 2 -e 1 -q 2 -P $procs

            elif [ $format = "unpaired" ]; then
                  # run FUCHS
                  FUCHS -N $sample_name -D $dcc_out_dir/CircRNACount -B $main_bam -A $bed_file -O $main_out/${sample_name} -J $main_junction -T $tmp_folder -p ensembl
-r 2 -e 1 -q 2 -P $procs
            fi

            conda deactivate
      }

############################################################################################

      export -f circtools_reconstruct

      # Parallel reconstruction:
      echo "Running FUCHS."
      parallel -j1 --xapply circtools_reconstruct {} $out_dir/circtools/03_reconstruct $exon_annotation ./ ./output $mode $procs :::: $out_dir/.tmp/names.tmp

fi

