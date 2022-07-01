#!/bin/bash

########## BATCH Lines for Resource Request ##########
#SBATCH --time=04:00:00		# limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1		# number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1		# number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=50	# number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=80G		# Specify the real memory required per node. (20G)
#SBATCH --job-name TriSE1	# you can give your job a name for easier identification (same as -J)
#SBATCH -A egl

### Modules required

#module purge
#module load GCC/8.2.0
#module load parallel-fastq-dump/0.6.5-Python-3.7.2
#module load Trimmomatic/0.39-Java-11
#module load GCC/4.7.2
#module load hisat2/2.1.0
#module load GCC/9.3.0
#module load SAMtools/1.11

################### Script ###################

######## Choose additional processes #######

donwloading=true;             ## Enabling the processes from downloading to read mapping
geno_index=false;		## Input the reference genome to indexing
Q20=false;                     ## Filter the mapped reads by 20 quality 
Unique_mapping=false;         ## Delete the multimapping reads
Drop_Duplicates=false;        ## Drop the duplicated reads
DeDuplicates_Q20=false;         ## Drop the duplicated Q20 reads

####### Number of ${threads} #######

threads=50;

####### Sample list ########

input_list='test_list.txt';

####### Reference Genome ########
indexing_genome(){

	if ${geno_index}; then
        path_genome='path/file/reference/genome';       ## Input the reference genome path file
	hisat2-build ${path_genome} ${path_genome//.fa/};
	index_genome=Index_${path_genome//.fa/};
	
	else
        index_genome='Index_TAIR10_Hisat2';     ## Input the base name of reference genome index
	fi
}

###################################### SCRIPT #########################################

############# Single functions ################

donwload_single() {
	
	line=$1
	parallel-fastq-dump -s ${line} -t ${threads} --gzip -O /mnt/home/gomezcan/Projects/Arabidopsis_Atlas/SecondPhase/;
}

cleaning_single(){
	line=$1
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads ${threads} ${line} Clean.${line} \
	ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
}

read_mapping_single(){

	file=$1;
        file_out=${file//Clean./};
	file_out=${file_out//.fastq.gz/.sam};
        hisat2 -p ${threads} --no-unal --max-intronlen 20000 -x ${index_genome} -U ${file} -S $file_out;

	samtools view -@ $((threads-10)) -h -S -b ${file_out} | samtools sort -@ $((threads-10)) -o ${file_out//.sam/.bam};
	rm ${file_out};
}

######################################################
######################################################
############### Paired functions #####################

donwload_paired() {
	line=$1;
	parallel-fastq-dump -s ${line} -t ${threads} --gzip --split-files -O /mnt/home/gomezcan/Projects/Arabidopsis_Atlas/SecondPhase/;
}

cleaning_paired() {

	file=$1
        file1=${file};
        file2=${file//1.fastq.gz/2.fastq.gz};

        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${threads}  $file1 $file2 \
        Clean.$file1 Unpaired.$file1 Clean.$file2 Unpaired.$file2 \
        ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
}

read_mapping_paired() {
	file=$1;

        file1=${file};
        file2=${file//1.fastq.gz/2.fastq.gz};
        file_out=${file//Clean./};
        file_out=${file_out//_1.fastq.gz/.sam};
        hisat2 -p ${threads} --no-unal --max-intronlen 20000 --no-discordant -x ${index_genome} -1 ${file1} -2 ${file2} -S ${file_out};


        samtools view -@ $((threads-10)) -h -S -b ${file_out} | samtools sort -@ $((threads-10)) -o ${file_out//.sam/.bam};
	rm ${file_out};
}

################################################################
################################################################
######################### Single Module ########################

single_module() {
	sample=$1;
	module purge
	ml GCC/8.2.0 parallel-fastq-dump/0.6.5-Python-3.7.2
	if $donwloading; then
		echo 'donwloading the sample';
		donwload_single ${sample} ;
	else
		echo 'dont donwloading';
	fi
	module purge
	ml GCC/10.2.0 Trimmomatic/0.39-Java-11 hisat2/2.1.0 SAMtools/1.11
	cleaning_single ${sample}.fastq.gz;

	read_mapping_single Clean.${sample}.fastq.gz;

}
#################################################################
#################################################################
######################### Paired Module #########################

paired_module() {
	sample=$1;
        module purge
        ml GCC/8.2.0 parallel-fastq-dump/0.6.5-Python-3.7.2
	if ${donwloading}; then
		echo 'Paired downloading the sample';
		donwload_paired ${sample};
	else
		echo 'dont downloading';
	fi

        module purge
        ml GCC/10.2.0 Trimmomatic/0.39-Java-11 hisat2/2.1.0 SAMtools/1.11
	cleaning_paired ${sample}_1.fastq.gz;

	read_mapping_paired Clean.${sample}_1.fastq.gz;
}

############################################################################
############################################################################
###################### Additional processes Module #########################

Q20_filter() {
	file=$1;
	samtools view -@ $((threads-10)) -q 20 ${file} > Q20.${file};
}

DeDu_Q20() {
	file=$1;
	samtools fixmate -@ $((threads-10)) -O bam -m ${file} - | samtools sort -@ $((threads-10)) -o sort.${file} - | samtools markdup \
	-@ $((threads-10)) -r - DeDu_Q20.${file};
}

DeDup() {
	file=$1;
	samtools fixmate -@ $((threads-10)) -O bam -m ${file} - | samtools sort -@ $((threads-10)) -o sort.${file} - | samtools markdup \
	-@ $((threads-10)) -r - DeDup.${file};
}

unique_mapped(){
	file=$1;
	samtools view -@ $((threads-10)) -h  ${file} | grep -v 'ZS:' | samtools view -@ $((threads-10)) -h -o unic.${file};
}



################################ COMPLETE MODULE ##########################

indexing_genome;
while read -r -a line; do

	file=${line[0]};
	if [[ ${line[1]} == "Single" ]]  ## Ask if the sample is single
	then
		echo 'executing single'
		single_module ${file}; #Execute the single script

		###### The additional processes ###########
		if ${Q20}; then
			Q20_filter ${file}.bam;
		fi

		if ${DeDuplicates_Q20}; then
			if ${Q20}; then 
				Dedu_Q20 Q20.${file}.bam;
			else
				echo 'The Q20 filtered data must be generated first';
			fi
		fi

		if ${Drop_Duplicates}; then
		DeDup ${file}.bam;
		fi

		if ${Unique_mapping}; then
			unique_mapped ${file}.bam;
		fi

	elif [[ ${line[1]} == "Paired" ]]   ## Ask if the sample is paired
	then
		echo 'executing paired'
		paired_module ${file}; #Execute the paired script
		
		##### The additonal processes #########
             	if ${Q20}; then
                        Q20_filter ${file}.bam;
                fi

                if ${DeDuplicates_Q20}; then
                        if ${Q20}; then 
                                Dedu_Q20 Q20.${file}.bam;
                        else
                                echo 'The Q20 filtered data must be generated first';
                        fi
                fi

                if ${Drop_Duplicates}; then
                        DeDup ${file}.bam;
                fi

                if ${Unique_mapping}; then
                        unique_mapped ${file}.bam;
                fi

	else
		echo "Error when input samples name"
	fi
done < ${input_list}
