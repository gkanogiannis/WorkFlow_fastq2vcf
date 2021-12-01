#
# VariantDiscovery_WF
# Cromwell workflow for calling variants with GATK
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

import "HaplotypeCaller.wdl" as subHC
import "GenotypeGVCFs.wdl" as subGT

workflow VariantDiscovery_WF {
	String genome_ref

	String listChromosomeFile
	String listChromosomeLengthFile
	Array[String] listChromosome = read_lines(listChromosomeFile)
	String? listScaffoldFile
	String? listScaffoldLengthFile
	Boolean withScaffolds
	Boolean? onlyGVCF=false

	String inputSamplesFile
	Array[Array[String]] inputSamples= read_tsv(inputSamplesFile)

	String basename

	String workDir
	String log
	String tmp
	
	String threads
	String java
	String picard
	String gatk
	String tabix

	String gvcf_out = "${workDir}/gvcf"
	String gDB_out = "${workDir}/gDB"
	String vcf_out = "${workDir}/vcf"

	call checkExistGvcf {
		input:
			samples = inputSamplesFile,
			workDir = gvcf_out,
			basename=basename,
			log = log
	}

	call checkExistGvcfIndex {
		input:
			samples = inputSamplesFile,
			workDir = gvcf_out,
			basename=basename,
			log = log
	}

	call checkNotExist {
		input:
			samples = inputSamplesFile,
			workDir = gvcf_out,
			basename=basename,
			log = log
	}

	scatter (sample in checkNotExist.notExist) {
		if(length(sample)>1){
			String name = sample[0]
			String bam = sample[1]
			String bamIndex = sample[2]

			call subHC.HaplotypeCaller_WF as HaplotypeCaller {
				input:
					genome_ref = genome_ref,
					name = name,
					bam = bam,
					bamIndex = bamIndex,
					listChromosomeFile = listChromosomeFile,
					listScaffoldFile = listScaffoldFile,
					withScaffolds = withScaffolds,
					workDir = gvcf_out,
					basename=basename,
					log = log,
					tmp = tmp,
					threads = threads,
					gatk = gatk
			}
		}
	}

	# This will run when there are HC unpocessed samples with some gvcf files nonexist
	if(defined(HaplotypeCaller.gvcf[0])){
		call concat as concatGvcf {
			input:
				a = checkExistGvcf.exist,
				b = select_all(HaplotypeCaller.gvcf),
				basename=basename,
				log = log
		}

		call concat as concatGvcfIndex {
			input:
				a = checkExistGvcfIndex.exist,
				b = select_all(HaplotypeCaller.gvcfIndex),
				basename=basename,
				log = log
		}

		if(!onlyGVCF){
			scatter (chr in listChromosome) {
				call GenomicsDB_Chr as gDB_Chr_Exist_nonExist {
					input:
						outDir = gDB_out,
						log = log,
						tmp = tmp,
						chr = chr,
						gvcfs = concatGvcf.out,
						gvcfsIndex = concatGvcfIndex.out,
						basename=basename,
						threads = 8,
						gatk = gatk
				}
			}

			if(withScaffolds){
				call GenomicsDB_Sca as gDB_Sca_Exist_nonExist {
					input:
						outDir = gDB_out,
						log = log,
						tmp = tmp,
						listScaffoldFile = listScaffoldFile,
						gvcfs = concatGvcf.out,
						gvcfsIndex = concatGvcfIndex.out,
						basename=basename,
						threads = 8,
						gatk = gatk
				}
			}

			call subGT.GenotypeGVCFs_WF as GT_Exist_nonExist {
				input:
					genome_ref = genome_ref,
					listChromosomeFile = listChromosomeFile,
					listChromosomeLengthFile = listChromosomeLengthFile,
					listScaffoldFile = listScaffoldFile,
					listScaffoldLengthFile = listScaffoldLengthFile,
					withScaffolds = withScaffolds,
					basename=basename,
					workDir = vcf_out,
					gDB_base=gDB_out,
					gDB_chr=gDB_Chr_Exist_nonExist.gDB,
					gDB_Sca = gDB_Sca_Exist_nonExist.gDB,
					log = log,
					tmp = tmp,
					threads = threads,
					gatk = gatk,
					tabix = tabix
			}
		}
	}

	# This will run when there aren't HC unpocessed samples and all gvcf files exist
	if(!defined(HaplotypeCaller.gvcf[0]) && !onlyGVCF){
		scatter (chr in listChromosome) {
			call GenomicsDB_Chr as gDB_Chr_Exist {
				input:
					outDir = gDB_out,
					log = log,
					tmp = tmp,
					chr = chr,
					gvcfs = checkExistGvcf.exist,
					gvcfsIndex = checkExistGvcfIndex.exist,
					basename=basename,
					threads = 8,
					gatk = gatk
			}
		}
		if(withScaffolds){
			call GenomicsDB_Sca as gDB_Sca_Exist {
				input:
					outDir = gDB_out,
					log = log,
					tmp = tmp,
					listScaffoldFile = listScaffoldFile,
					gvcfs = checkExistGvcf.exist,
					gvcfsIndex = checkExistGvcfIndex.exist,
					basename=basename,
					threads = 8,
					gatk = gatk
			}
		}

		call subGT.GenotypeGVCFs_WF as GT_Exist {
			input:
				genome_ref = genome_ref,
				listChromosomeFile = listChromosomeFile,
				listChromosomeLengthFile = listChromosomeLengthFile,
				listScaffoldFile = listScaffoldFile,
				listScaffoldLengthFile = listScaffoldLengthFile,
				withScaffolds = withScaffolds,
				basename=basename,
				workDir = vcf_out,
				gDB_base=gDB_out,
				gDB_chr=gDB_Chr_Exist.gDB,
				gDB_Sca = gDB_Sca_Exist.gDB,
				log = log,
				tmp = tmp,
				threads = threads,
				gatk = gatk,
				tabix = tabix
		}
	}
}

############################################################################
task checkExistGvcf {
	File samples
	String workDir
	String basename
	String log

	command <<<
		mkdir -p ${log} ;
		while IFS=$'\t' read -r i j k ; do 
			GVCF="${workDir}/$i.g.vcf.gz" ;
			GVCFINDEX="${workDir}/$i.g.vcf.gz.tbi" ;
			if [[ -f $GVCF && -f $GVCFINDEX ]]; then
				echo -e "$GVCF"
			fi ;
		done < ${samples} ;
		exit 0 ;
	>>>
	runtime { 
		cpu: 2
		memory: "512 MB"
		my_job_name: "VD_CheckExistGvcf_${basename}"
		my_log: "${log}/VD_CheckExistGvcf_${basename}"
	}

	output {
		Array[File] exist = read_lines(stdout())
	}
}

task checkExistGvcfIndex {
	File samples
	String workDir
	String basename
	String log

	command <<<
		mkdir -p ${log} ;
		while IFS=$'\t' read -r i j k ; do 
			GVCF="${workDir}/$i.g.vcf.gz" ;
			GVCFINDEX="${workDir}/$i.g.vcf.gz.tbi" ;
			if [[ -f $GVCF && -f $GVCFINDEX ]]; then
				echo -e "$GVCFINDEX"
			fi ;
		done < ${samples} ;
		exit 0 ;
	>>>
	runtime { 
		cpu: 2
		memory: "512 MB"
		my_job_name: "VD_CheckExistGvcfIndex_${basename}"
		my_log: "${log}/VD_CheckExistGvcfIndex_${basename}"
	}

	output {
		Array[File] exist = read_lines(stdout())
	}
}

task checkNotExist {
	File samples
	String workDir
	String basename
	String log

	command <<<
		mkdir -p ${log} ;
		while IFS=$'\t' read -r i j k ; do 
			GVCF="${workDir}/$i.g.vcf.gz" ;
			GVCFINDEX="${workDir}/$i.g.vcf.gz.tbi" ;
			if [[ ! -f $GVCF || ! -f $GVCFINDEX ]]; then
				echo -e "$i\t$j\t$k"
			fi ;
		done < ${samples} ;
		exit 0 ;
	>>>
	runtime { 
		cpu: 2
		memory: "512 MB"
		my_job_name: "VD_CheckNotExist_${basename}"
		my_log: "${log}/VD_CheckNotExist_${basename}"
	}

	output {
		Array[Array[String]] notExist = read_tsv(stdout())
	}
}

task concat {
	Array[File] a
	Array[File] b
	String basename
	String log

	command <<<
		cat ${write_lines(a)} ;
		cat ${write_lines(b)} ;
		exit 0 ;
	>>>
	runtime { 
		cpu: 2
		memory: "512 MB"
		my_job_name: "VD_Concat_${basename}"
		my_log: "${log}/VD_Concat_${basename}"
	}
	
	output {
		Array[File] out = read_lines(stdout())
	}
}

task GenomicsDB_Chr {
	String outDir
	String chr
	Array[File] gvcfs
	Array[File] gvcfsIndex
	String basename
	String log
	String tmp
	String gatk

	Int threads
	Int mem=35000
	Int mems=34000

	command <<<
		mkdir -p ${tmp} ${outDir};
		FILE="${outDir}/${chr}/callset.json"
		if [[ ! -f $FILE ]]; then 
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
				--genomicsdb-workspace-path ${outDir}/${chr} \
				--tmp-dir ${tmp} \
				--reader-threads 2 \
				-L ${chr} \
				--variant ${sep=' --variant ' gvcfs} ;
		else
			echo -e "${outDir}/${chr} exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VD_GenomicsDB_${basename}_${chr}"
		my_log: "${log}/VD_GenomicsDB_${basename}_${chr}"
		maxRetries: 1
	}

	output {
		File gDB = "${outDir}/${chr}"
	}
}

task GenomicsDB_Sca {
	String outDir
	File listScaffoldFile
	Array[File] gvcfs
	Array[File] gvcfsIndex
	String basename
	String log
	String tmp
	String gatk

	Int threads
	Int mem=35000
	Int mems=34000

	command <<<
		mkdir -p ${tmp} ${outDir};
		FILE="${outDir}/Scaffolds/callset.json"
		if [[ ! -f $FILE ]]; then 
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" GenomicsDBImport \
				--genomicsdb-workspace-path ${outDir}/Scaffolds \
				--tmp-dir ${tmp} \
				--reader-threads 2 \
				-L ${listScaffoldFile} \
				--variant ${sep=' --variant ' gvcfs} ;
		else
			echo -e "${outDir}/Scaffolds exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VD_GenomicsDB_${basename}_Scaffolds"
		my_log: "${log}/VD_GenomicsDB_${basename}_Scaffolds"
		maxRetries: 1
	}

	output {
		File gDB = "${outDir}/Scaffolds"
	}
}
