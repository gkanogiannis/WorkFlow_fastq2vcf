#
# VariantFiltering_WF
# Cromwell workflow for filtering variants
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

workflow VariantFiltering_WF {
	String genome_ref
	String genome_ref_index
	String genome_ref_dict
	String? repeats

	String inputVcfFile
	String inputVcfFileIndex
	String basename

	Boolean? onlyBiallelic=true
	Int? filterDPValue=3
	Float? filterMAFValue=0.05
	Float? filterMissingValue=0.90
	Float? filterLDValue=0.50

	String filterDensityScript
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	String workDir
	String log
	String tmp
	
	String threads
	String java
	String picard
	String gatk
	String tabix
	String vcftools
	String bcftools
	String groovy

	String vcf_out = "${workDir}/vcf"

	call getSNPs {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = inputVcfFile,
				inputVcfFileIndex = inputVcfFileIndex,
				basename = basename,
				onlyBiallelic = onlyBiallelic,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				bcftools = bcftools,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
	}

	if(defined(repeats)) {
		call getRepeats {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				repeats = repeats,
				inputVcfFile = inputVcfFile,
				inputVcfFileIndex = inputVcfFileIndex,
				basename = basename,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools
		}

		call filterINFO {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = getSNPs.vcf,
				inputVcfFileIndex = getSNPs.index,
				repMaskVcfFile = getRepeats.vcf,
				repMaskVcfFileIndex = getRepeats.index,
				basename = basename(getSNPs.vcf,".vcf.gz"),
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterDP {
			input:
				inputVcfFile = filterINFO.vcf,
				inputVcfFileIndex = filterINFO.index,
				basename = basename(filterINFO.vcf,".vcf.gz"),
				filterDPValue = filterDPValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterMissing {
			input:
				inputVcfFile = filterDP.vcf,
				inputVcfFileIndex = filterDP.index,
				basename = basename(filterDP.vcf,".vcf.gz"),
				filterMissingValue = filterMissingValue,
				filterDensityScript = filterDensityScript,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				groovy = groovy,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterMAF {
			input:
				inputVcfFile = filterMissing.vcf,
				inputVcfFileIndex = filterMissing.index,
				basename = basename(filterMissing.vcf,".vcf.gz"),
				filterMAFValue = filterMAFValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}
	}

	if(!defined(repeats)) {
		call filterINFO_noRM {
			input:
				ref = genome_ref,
				ref_index = genome_ref_index,
				ref_dict = genome_ref_dict,
				inputVcfFile = getSNPs.vcf,
				inputVcfFileIndex = getSNPs.index,
				basename = basename(getSNPs.vcf,".vcf.gz"),
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				gatk = gatk,
				vcftools = vcftools,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterDP as filterDP_noRM {
			input:
				inputVcfFile = filterINFO_noRM.vcf,
				inputVcfFileIndex = filterINFO_noRM.index,
				basename = basename(filterINFO_noRM.vcf,".vcf.gz"),
				filterDPValue = filterDPValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterMissing as filterMissing_noRM {
			input:
				inputVcfFile = filterDP_noRM.vcf,
				inputVcfFileIndex = filterDP_noRM.index,
				basename = basename(filterDP_noRM.vcf,".vcf.gz"),
				filterMissingValue = filterMissingValue,
				filterDensityScript = filterDensityScript,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				groovy = groovy,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}

		call filterMAF as filterMAF_noRM {
			input:
				inputVcfFile = filterMissing_noRM.vcf,
				inputVcfFileIndex = filterMissing_noRM.index,
				basename = basename(filterMissing_noRM.vcf,".vcf.gz"),
				filterMAFValue = filterMAFValue,
				outDir = vcf_out,
				log = log,
				tmp = tmp,
				threads = 2,
				vcftools = vcftools,
				tabix = tabix,
				calculateStatsScript = calculateStatsScript,
				calculateIStatsScript = calculateIStatsScript,
				calculateHistogramsScript = calculateHistogramsScript
		}
	}

}

############################################################################
task getRepeats {
	File ref
	File ref_index
	File ref_dict
	File repeats
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.reps.vcf.gz" ;
		INDEX="${outDir}/${basename}.reps.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				-L ${repeats} \
 				-V ${inputVcfFile} \
 				-O ${outDir}/${basename}.reps.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.reps.vcf.gz > ${outDir}/${basename}.reps.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.reps.vcf.gz exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_GetRepeats"
		my_log: "${log}/VF_GetRepeats"
	}

	output {
		File vcf = "${outDir}/${basename}.reps.vcf.gz"
		File index = "${outDir}/${basename}.reps.vcf.gz.tbi"
	}
}

task getSNPs {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	Boolean? onlyBiallelic=true
	String bcftools
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=4
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.snps.vcf.gz" ;
		INDEX="${outDir}/${basename}.snps.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			if ${onlyBiallelic} ; then
				ALLELES="-m2 -M2 --types snps" ;
			else
				ALLELES="--types snps" ;
			fi ;
			${bcftools} view \
				$ALLELES \
				-Oz -o ${outDir}/${basename}.snps.vcf.gz \
				${inputVcfFile} ;
			${tabix} -p vcf ${outDir}/${basename}.snps.vcf.gz ;
			${vcftools} --gzvcf ${outDir}/${basename}.snps.vcf.gz > ${outDir}/${basename}.snps.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.snps.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}.snps.vcf.gz && ! -f ${outDir}/${basename}.snps.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}.snps.vcf.gz ${outDir}/${basename}.snps.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}.snps.vcf.stats && ! -f ${outDir}/${basename}.snps.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}.snps.vcf.stats ${outDir}/${basename}.snps.vcf ;
		fi ;
		#if [[ -f ${outDir}/${basename}.snps.vcf.gz && ! -f ${outDir}/${basename}.snps.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.snps.vcf.gz ${outDir}/${basename}.snps.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_GetSNPs"
		my_log: "${log}/VF_GetSNPs"
	}

	output {
		File vcf = "${outDir}/${basename}.snps.vcf.gz"
		File index = "${outDir}/${basename}.snps.vcf.gz.tbi"
	}
}

task filterINFO {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	File repMaskVcfFile
	File repMaskVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.filter_info.vcf.gz" ;
		INDEX="${outDir}/${basename}.filter_info.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" VariantFiltration \
				-R ${ref} \
 				-V ${inputVcfFile} \
 				--filter-name "LowQD" --filter-expression "QD<2.0" \
 				--filter-name "LowMQ" --filter-expression "MQ<40.0" \
 				--filter-name "HiFS" --filter-expression "FS>60.0" \
 				--filter-name "HiSOR" --filter-expression "SOR>3.0" \
 				--filter-name "LowMQRankSum" --filter-expression "MQRankSum<-12.5" \
 				--filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum<-8.0" \
 				--mask-name "RepMask" --mask ${repMaskVcfFile} \
 				-O ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				--exclude-filtered \
 				-V ${outDir}/${basename}.filter_info_marked.vcf.gz \
 				-O ${outDir}/${basename}.filter_info.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.filter_info.vcf.gz > ${outDir}/${basename}.filter_info.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.filter_info.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}.filter_info.vcf.stats && ! -f ${outDir}/${basename}.filter_info.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}.filter_info.vcf.stats ${outDir}/${basename}.filter_info.vcf ;
		fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterINFO"
		my_log: "${log}/VF_FilterINFO"
	}

	output {
		File vcf = "${outDir}/${basename}.filter_info.vcf.gz"
		File index = "${outDir}/${basename}.filter_info.vcf.gz.tbi"
	}
}

task filterINFO_noRM {
	File ref
	File ref_index
	File ref_dict
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	String outDir
	String log
	String tmp
	String gatk
	String vcftools
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}.filter_info.vcf.gz" ;
		INDEX="${outDir}/${basename}.filter_info.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" VariantFiltration \
				-R ${ref} \
 				-V ${inputVcfFile} \
 				--filter-name "LowQD" --filter-expression "QD<2.0" \
 				--filter-name "LowMQ" --filter-expression "MQ<40.0" \
 				--filter-name "HiFS" --filter-expression "FS>60.0" \
 				--filter-name "HiSOR" --filter-expression "SOR>3.0" \
 				--filter-name "LowMQRankSum" --filter-expression "MQRankSum<-12.5" \
 				--filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum<-8.0" \
 				-O ${outDir}/${basename}.filter_info_marked.vcf.gz ;
 			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" SelectVariants \
				-R ${ref} \
				--exclude-filtered \
 				-V ${outDir}/${basename}.filter_info_marked.vcf.gz \
 				-O ${outDir}/${basename}.filter_info.vcf.gz ;
 			${vcftools} --gzvcf ${outDir}/${basename}.filter_info.vcf.gz > ${outDir}/${basename}.filter_info.vcf.log 2>&1 ;
		else
			echo -e "${outDir}/${basename}.filter_info.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}.filter_info.vcf.stats && ! -f ${outDir}/${basename}.filter_info.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}.filter_info.vcf.stats ${outDir}/${basename}.filter_info.vcf ;
		fi ;
		#if [[ -f ${outDir}/${basename}.filter_info.vcf.gz && ! -f ${outDir}/${basename}.filter_info.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}.filter_info.vcf.gz ${outDir}/${basename}.filter_info.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterINFO"
		my_log: "${log}/VF_FilterINFO"
	}

	output {
		File vcf = "${outDir}/${basename}.filter_info.vcf.gz"
		File index = "${outDir}/${basename}.filter_info.vcf.gz.tbi"
	}
}

task filterDP {
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	Int filterDPValue
	String outDir
	String log
	String tmp
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}_DP${filterDPValue}.vcf.gz" ;
		INDEX="${outDir}/${basename}_DP${filterDPValue}.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${vcftools} --gzvcf ${inputVcfFile} --minDP ${filterDPValue} --recode --recode-INFO-all --stdout \
				| bgzip -c > ${outDir}/${basename}_DP${filterDPValue}.vcf.gz ;

 			${vcftools} --gzvcf ${outDir}/${basename}_DP${filterDPValue}.vcf.gz > ${outDir}/${basename}_DP${filterDPValue}.vcf.log 2>&1 ;
 			${tabix} -p vcf ${outDir}/${basename}_DP${filterDPValue}.vcf.gz ;
		else
			echo -e "${outDir}/${basename}_DP${filterDPValue}.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}_DP${filterDPValue}.vcf.gz && ! -f ${outDir}/${basename}_DP${filterDPValue}.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}_DP${filterDPValue}.vcf.gz ${outDir}/${basename}_DP${filterDPValue}.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}_DP${filterDPValue}.vcf.stats && ! -f ${outDir}/${basename}_DP${filterDPValue}.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}_DP${filterDPValue}.vcf.stats ${outDir}/${basename}_DP${filterDPValue}.vcf ;
		fi ;
		#if [[ -f ${outDir}/${basename}_DP${filterDPValue}.vcf.gz && ! -f ${outDir}/${basename}_DP${filterDPValue}.vcf.istats ]]; then
		#	bash ${calculateIStatsScript} ${outDir}/${basename}_DP${filterDPValue}.vcf.gz ${outDir}/${basename}_DP${filterDPValue}.vcf.istats 32 ;
		#fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterDP"
		my_log: "${log}/VF_FilterDP"
	}

	output {
		File vcf = "${outDir}/${basename}_DP${filterDPValue}.vcf.gz"
		File index = "${outDir}/${basename}_DP${filterDPValue}.vcf.gz.tbi"
	}
}

task filterMissing {
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	Float filterMissingValue
	String filterDensityScript
	String outDir
	String log
	String tmp
	String vcftools
	String groovy
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}_den${filterMissingValue}.vcf.gz" ;
		INDEX="${outDir}/${basename}_den${filterMissingValue}.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			export JAVA_OPTS="$JAVA_OPTS -Xms${mems}M -Xmx${mem}M" ;
			${groovy} ${filterDensityScript} ${inputVcfFile} ${outDir}/${basename}_den${filterMissingValue}.vcf ${filterMissingValue} ;
			bgzip ${outDir}/${basename}_den${filterMissingValue}.vcf ;

 			${vcftools} --gzvcf ${outDir}/${basename}_den${filterMissingValue}.vcf.gz > ${outDir}/${basename}_den${filterMissingValue}.vcf.log 2>&1 ;
 			${tabix} -p vcf ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ;
		else
			echo -e "${outDir}/${basename}_den${filterMissingValue}.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.gz && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.stats && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.stats ${outDir}/${basename}_den${filterMissingValue}.vcf ;
		fi ;
		if [[ -f ${outDir}/${basename}_den${filterMissingValue}.vcf.gz && ! -f ${outDir}/${basename}_den${filterMissingValue}.vcf.istats ]]; then
			bash ${calculateIStatsScript} ${outDir}/${basename}_den${filterMissingValue}.vcf.gz ${outDir}/${basename}_den${filterMissingValue}.vcf.istats 32 ;
		fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterMissing"
		my_log: "${log}/VF_FilterMissing"
	}

	output {
		File vcf = "${outDir}/${basename}_den${filterMissingValue}.vcf.gz"
		File index = "${outDir}/${basename}_den${filterMissingValue}.vcf.gz.tbi"
	}
}

task filterMAF {
	File inputVcfFile
	File inputVcfFileIndex
	String basename
	Float filterMAFValue
	String outDir
	String log
	String tmp
	String vcftools
	String tabix
	String calculateStatsScript
	String calculateIStatsScript
	String calculateHistogramsScript

	Int? threads=2
	Int mem=30000
	Int mems=25000

	command <<<
		mkdir -p ${log} ${tmp} ${outDir} ;
		FILE="${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz" ;
		INDEX="${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz.tbi" ;
		if [[ ! -f $FILE || ! -f $INDEX ]]; then
			${vcftools} --gzvcf ${inputVcfFile} --maf ${filterMAFValue} --recode --recode-INFO-all --stdout \
				| bgzip -c > ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ;

 			${vcftools} --gzvcf ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz > ${outDir}/${basename}_MAF${filterMAFValue}.vcf.log 2>&1 ;
 			${tabix} -p vcf ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ;
		else
			echo -e "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz exists" ;
		fi ;

		if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ]]; then
			bash ${calculateStatsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ;
		fi ;
		if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.Ho_He_PIC_MAF_Mis.hist.pdf ]]; then
			bash ${calculateHistogramsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.stats ${outDir}/${basename}_MAF${filterMAFValue}.vcf ;
		fi ;
		if [[ -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz && ! -f ${outDir}/${basename}_MAF${filterMAFValue}.vcf.istats ]]; then
			bash ${calculateIStatsScript} ${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz ${outDir}/${basename}_MAF${filterMAFValue}.vcf.istats 32 ;
		fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "VF_FilterMAF"
		my_log: "${log}/VF_FilterMAF"
	}

	output {
		File vcf = "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz"
		File index = "${outDir}/${basename}_MAF${filterMAFValue}.vcf.gz.tbi"
	}
}
