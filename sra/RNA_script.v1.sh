#!/bin/bash

## Modules required 
#module purge
#module load GCCcore/8.2.0
#module load GGCcore/4.7.2
#module load hisat2/2.1.0
#module load SAMtools/1.11
#module load parallel-fastq-dump/0.6.5-Python-3.7.2
#module load Trimmomatic/0.39-Java-11
#module load FastQC/0.11.7-Java-1.8.0_162

################### Script ###################

######## Choose additional processes #######

processing = true;		## Enabling the processes from downloading to read mapping
fastqc_report = false;           ## Ask the fastqc report for every enabled filter process
Q20 = false;			## Filter the mapped reads by 20 quality
Unique_mapping = false;		## Delete the multimapping reads
Drop_Duplicates = false;	## Drop the duplicated reads
DeDuplicates_Q20 false;		## Drop the duplicated Q20 reads | Remember to set Q20 variable to "true"

####### Reference Genome #######

*** Pending index condition
path_genome = "path/file/genome";   # Input the Reference genome path file

####### Paired or single ends #######

list_file = "file_list/path/file";	## Input the donwload list path file

*** read -r firtsline < ${list_file};	## Read the first line on the list


####### Number of threads #######

threads = 50;	# It's a recurrent parameter

####### Choose the way #######

###################### SINGLE END SCRIPT ##########################
#### Genome Indexing ####
	
# hisat2-build ${path_genome} index_genome;

if [ ${firstline[1]} == "Single" ]
	then
	#### Donwload the data####

	if ${processing}; then

	while read -r line; do 
	parallel-fastq-dump --gzip -t ${threads} -s $line -O . ;
	
	echo 'cleaning' $line
	
	java -jar trimmomatic-0.39.jar SE -t ${threads} ${line}.fastq.gz Clean.${line}.fastq.gz ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
	echo 'mapping' $line
	hisat2 -p ${threads} --no-unal --add-chrname --no-mixed --no-discordant -x index_genome -U $file -S $file_out;
	
	samtools view -@ ${threads} -h -S -b ${file_out} | samtools sort -@ ${threads} -o sorted.${file_out//.sam/.bam};
	
	# $line.fastq.gz
	# $line_1.fastq.gz
	# $line_2.fastq.gz

	done; 
	./my_counter


	#### Read Mapping ####

	for file in Clean.*; do
	file_out = ${file//Clean./};
	file_out = ${file_out//fastq.gz/.sam};
	hisat2 -p ${threads} --no-unal --add-chrname --no-mixed --no-discordant -x index_genome -U $file -S $file_out;
	samtools view -@ ${threads} -h -S -b ${file_out} | samtools sort -@ ${threads} -o ${file_out//.sam/.bam};
	done;

	fi;

	#### Q20 filtered Mapped reads ####
	if ${Q20}; then
		for file in *.bam; do
		file_out = ${file//sorted./};
		samtools view -@ ${threads} -q 20 $file > Q20.${file_out};
		done;		
	
		if ${fastqc_report}; then
			for file in Q20.*; do
			fastqc -o fastqc.$file -t ${threads} -f bam  $file;
			done;
		fi;
	fi;

	### Drop Duplicated Q20 Mapped reads ####

	if ${DeDuplicates_Q20}; then
		for file in Q20*; do
		samtools fixmate -@ ${threads} -O bam -m $file - | samtools sort -@ ${threads} -o sort.$file - | samtools markdup \
		-@ ${threads} -r - DeDu.$file;
		done;
		
		if ${fastqc_report}; then
			for file in DeDu.*; do
			fastqc -o fastqc.$file -t ${threads} -f bam $file;
			done;
		fi;
	fi;

	#### Drop Duplicated reads ####

	if ${Drop_Duplcates}; then
		for file in sorted*; do
		file_out = ${file//sorted./DeDuplicate};
		samtools fixmate -@ ${threads} -O bam -m $file - | samtools sort -@ ${threads} -o sort.$file - | samtools markdup \
		-@ ${threads} -r - ${file_out};
		done;

                if ${fastqc_report}; then
                        for file in DeDuplicate.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam $file;
                        done;
                fi;
	fi; 
	
	#### Unic Mapping Reads ####

	if ${unic_mapping}; then
		for file in sorted.*; do
		samtools view -@ ${threads} -h  $file | grep -v 'ZS:' | samtools view -@ ${threads} -h -o unic.$file;
		done;

                if ${fastqc_report}; then
                        for file in unic.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam $file;
                        done;
                fi;

	fi;


############################## PAIRED END SCRIPT ################################

else

	if ${processing}; then

	##### Donwload the data ######

	while red -r line; do
	parallel-fastq-dump --gzip -t 30 -s $line --split-files -O .;
	done < ${list_file};

	##### Clean the files #######

	for file in *1.fastq.gz; do
	file1 = $file;
	file2 = ${file//1.fastq.gz/2.fastq.gz};

	java -jar trimmomatic-0.39.jar PE -threads ${threads}  $file1 $file2 \
	Clean.$file1 Unpaired.$file1 Clean.$file2 Unpaired.$file2 \
	ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
	done;

	##### Genome Indexing #######

	hisat2-build -p ${threads} ${path_genome} Index_AGPv4;

	##### Reads Mapping ########

	for file in Clean.*; do
	file1 = $file
	file2 = ${file//1.fastq.gz/2.fastq.gz};
	file_out = ${file//Clean./};
	file_out = ${file_out//1.fastq.gz/.sam};
	hisat2 -p ${threads} --no-unal --add-chrname --no-mixed --no-discordant -x Index_AGPv4 -1 ${file1} -2 ${file2} -S ${file_out};
	
	samtools view -@ ${threads} -h -S -b ${file_out} | samtools sort -@ ${threads} -o sorted.${file_out//.sam/.bam};
	done;

	fi;
	
       #### Q20 filtered Mapped reads ####
        if ${Q20}; then
                for file in *.bam; do
                file_out = ${file//sorted./};
                samtools view -@ ${threads} -q 20 $file > Q20.${file_out};
                done;

                if ${fastqc_report}; then
                        for file in Q20.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam  $file;
                        done;
                fi;
        fi;

        ### Drop Duplicated Q20 Mapped reads ####

        if ${DeDuplicates_Q20}; then
                for file in Q20*; do
                samtools fixmate -@ ${threads} -O bam -m $file - | samtools sort -@ ${threads} -o sort.$file - | samtools markdup \
		-@ ${threads} -r - DeDu.$file;
                done;

                if ${fastqc_report}; then
                        for file in DeDu.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam $file;
                        done;
                fi;
        fi;

        #### Drop Duplicated reads ####

        if ${Drop_Duplcates}; then
                for file in sorted*; do
                file_out = ${file//sorted./DeDuplicate};
                samtools fixmate -@ ${threads} -O bam -m $file - | samtools sort -@ ${threads} -o sort.$file - | samtools markdup \ 
		-@ ${threads} -r - ${file_out};
                done;

                if ${fastqc_report}; then
                        for file in DeDuplicate.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam $file;
                        done;
                fi;
        fi;

        #### Unic Mapping Reads ####

        if ${unic_mapping}; then
                for file in sorted.*; do
                samtools view -@ ${threads} -h  $file | grep -v 'ZS:' | samtools view -@ ${threads} -h -o unic.$file;
                done;

                if ${fastqc_report}; then
                        for file in unic.*; do
                        fastqc -o fastqc.$file -t ${threads} -f bam $file;
                        done;
                fi;

        fi;
fi

# Note
#for i in SRR_*; do n=$(wc -l $i | cut -d' ' -f1); echo $n | awk -v f=$i '{print f, $1/4, $1/4}' ; done;
