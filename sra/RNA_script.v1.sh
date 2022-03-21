
### Modules required

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

donwloading= true;             ## Enabling the processes from downloading to read mapping
geno_index= true;		## Input the reference genome to indexing
Q20= true;                     ## Filter the mapped reads by 20 quality 
Unique_mapping= false;         ## Delete the multimapping reads
Drop_Duplicates= false;        ## Drop the duplicated reads
DeDuplicates_Q20= false;         ## Drop the duplicated Q20 reads

####### Reference Genome ########

if ${geno_index}; then
	path_genome='path/file/reference/genome';	## Input the reference genome path file
else
	index_genome='base_name/reference/genome';	## Input the base name of reference genome index
fi

#if ${geno_index}; then
#		indexing_genome ${path_genome};
#fi;

###################################### SCRIPT #########################################

############# Single functions ################

donwload_single () {
	line= $1;
	parallel-fastq-dump --gzip -t 50 -s ${line} -O .
}

cleaning_single() {
	line= $1;
	java -jar trimmomatic-0.39.jar SE -t 50 ${line} Clean.${line} ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
}

indexing_genome(){
	genome= $1;
	hisat2-build ${genome} index_genome;
	index_genome= 'index_genome';
}

read_mapping_single(){
	file=$1;
        file_out=${file//Clean./};
	file_out=${file_out//fastq.gz/.sam};
        hisat2 -p 50 --max-intronlen 20000 --no-unal -x ${index_genome} -U ${file} -S $file_out;

        samtools view -@ 50 -h -S -b ${file_out} | samtools sort -@ 50 -o ${file_out//.sam/.bam};
}

######################################################
######################################################
############### Paired functions #####################

donwloading_paired() {
	line= $1;
	parallel-fastq-dump --gzip -t 50 -s ${line} --split-files -O .;
}

cleaning_paired() {

	file=$1
        file1=${file};
        file2=${file//1.fastq.gz/2.fastq.gz};

        java -jar trimmomatic-0.39.jar PE -threads 50  $file1 $file2 \
        Clean.$file1 Unpaired.$file1 Clean.$file2 Unpaired.$file2 \
        ILLUMINACLIP:Adapter.fastq:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30;
}

read_mapping_paired() {
	
	file= $1;

        file1= ${file};
        file2= ${file//1.fastq.gz/2.fastq.gz};
        file_out= ${file//Clean./};
        file_out= ${file_out//_1.fastq.gz/.sam};
        hisat2 -p 50 --max-intronlen 20000 --no-unal --no-discordant -x ${index_genome} -1 ${file1} -2 ${file2} -S ${file_out};

        samtools view -@ 50 -h -S -b ${file_out} | samtools sort -@ 50 -o ${file_out//.sam/.bam};

}

################################################################
################################################################
######################### Single Module ########################

single_module() {
	sample= $1;
	if ${downloading}; then
		download_single ${sample} ;
	fi;
	
	cleaning_single ${sample}.fastq.gz;

	read_mapping_single Clean.${sample}.fastq.gz;

}
#################################################################
#################################################################
######################### Paired Module #########################

paired_module() {
	sample= $1;
	if ${downloading}; then
		donwload_paired ${sample};
	fi;

	cleaning_paired ${sample}_1.fastq.gz;

	read_mapping_pared Clean.${sample}_1.fastq.gz;
}

############################################################################
############################################################################
###################### Additional processes Module #########################

Q20_filter() {
	file= $1;
	samtools view -@ 50 -q 20 ${file} > Q20.${file};
}

DeDu_Q20() {
	file= $1;
	samtools fixmate -@ 50 -O bam -m ${file} - | samtools sort -@ 50 -o sort.${file} - | samtools markdup -@ 50 -r - DeDu_Q20.${file};
}

DeDup() {
	file= $1;
	samtools fixmate -@ 50 -O bam -m ${file} - | samtools sort -@ 50 -o sort.${file} - | samtools markdup -@ 50 -r - DeDup.${file};
}

unique_mapped(){
	file= $1;
	samtools view -@ 50 -h  ${file} | grep -v 'ZS:' | samtools view -@ 50 -h -o unic.${file};
}



################################ COMPLETE MODULE ##########################

while read -r line; do

	file= ${line[0]};
	if [${line[1]} == "Single"]  ## Ask if the sample is single
	then
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

	elif [${line[1]}== "Paired"]   ## Ask if the sample is paired
	then

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
done;



