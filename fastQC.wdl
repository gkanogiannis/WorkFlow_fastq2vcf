#
# fastQC_WF
# Cromwell workflow for performing QC on samples fastq
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

workflow fastQC_WF {
	String workDir
	String log
	String tmp

	String inputSamplesFile
	String inputFastqsFile
	Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)
	Array[Array[String]] inputFastqs = read_tsv(inputFastqsFile)

	String fastqc

	scatter (i in range(length(inputSamples))) {
		call fastqcTask {
			input:
				name = inputSamples[i][0],
				fastq_file = inputFastqs[i][0],
				outDir = "${workDir}/fastQC",
				log = log,
				fastqc = fastqc
		}
	}
}

#######################################################################
task fastqcTask {
	String name 
	String fastq_file
	String outDir
	String log
	String fastqc

	command <<<
		mkdir -p ${outDir} ;
		if ls ${outDir}/${name}_*fastqc.zip 1> /dev/null 2>&1; then
			echo -e "fastqc exists" ;
			unzip -o ${outDir}/${name}_*fastqc.zip -d ${outDir} ;
		else
			${fastqc} -o ${outDir} ${fastq_file} ;
			unzip -o ${outDir}/${name}_*fastqc.zip -d ${outDir} ;
		fi ;
		enc=$(grep -i encoding `ls ${outDir}/${name}_*fastqc/fastqc_data.txt`) ;
		echo -e "$enc" > ${outDir}/${name}_encoding.txt ;
		exit 0 ;
	>>>

	runtime {
		cpu: "1"
		memory: "4000 MB"
		my_job_name: "fastqc_${name}"
		my_log: "${log}/fastqc_${name}"
	}

	output {
		File dir = "${outDir}"
	}
}
