//
// filterDensity.groovy
// Groovy script to fiter (g)vcf files according to average density(missingness) threshold
//
// Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//

// groovy filterDensity.groovy input.vcf(gz) output.vcf density_threshold

import groovy.transform.Field

format_index = 0
samples_size = 0

def inputvcf = args[0]
def outputvcfFile = new File(args[1])
def outputvcf = outputvcfFile.newWriter()
double density_thresh = args[2] as Double
def signature = "##source=Output of filterDensity.groovy script by Anestis Gkanogiannis to select the variants with the given average density ${density_thresh}."

def variants = new ArrayList<VariantLight>()
def variantsIndexKeep = new ArrayList<Integer>()
def comments = new ArrayList<String>()
def header

def reader	
if(inputvcf.endsWith('.vcf')){
	reader = new BufferedReader(new FileReader(inputvcf))
}
else if(inputvcf.endsWith('.vcf.gz')){
	reader = new BufferedReader(new InputStreamReader(new java.util.zip.GZIPInputStream(new FileInputStream(inputvcf)), "UTF-8"));
}	

println "Start reading variants from ${inputvcf}"
def variant_counter = 0
reader.eachLine { String line ->
	if(line.startsWith('##'))
		comments.add(line)
	else if(line.startsWith('#')) {
		header = line
		def counter = 0
		line.split('\\t').each {
			if(it.equals('FORMAT')) {
				format_index = counter
				counter = -1
			}
			else if(counter<0)
				samples_size++
			else
				counter++
		}
	}	
	else {
		def variant = new VariantLight(variantIndex: variant_counter++, formatIndex: format_index)
		variant.initialize(line)
		variants.add(variant)
		//println variant.density
		if(variant_counter%10000==0)
			println "Read ${variant_counter} variants."
	}	
}

println "Sorting..."
variants.sort{a,b->b.density<=>a.density}

println "Writing filtered variants to ${outputvcfFile}"
outputvcf << comments.join('\n')
outputvcf << '\n'
outputvcf << signature
outputvcf << '\n'
outputvcf << header
outputvcf << '\n'
double current_selected = 0.0
double current_density = 0.0
variants.find{it ->
	current_density += it.density
	current_selected++
	if(current_density/current_selected>=density_thresh){
		variantsIndexKeep.add(it.variantIndex)
		print ''
		//println "${current_selected} ${current_density/current_selected} ${it.variantIndex} ${it.density}"
	}
}
println variantsIndexKeep.size()
variant_counter = 0
if(inputvcf.endsWith('.vcf')){
	reader = new BufferedReader(new FileReader(inputvcf))
}
else if(inputvcf.endsWith('.vcf.gz')){
	reader = new BufferedReader(new InputStreamReader(new java.util.zip.GZIPInputStream(new FileInputStream(inputvcf)), "UTF-8"));
}	
reader.eachLine { String line ->
	if(!line.startsWith('#')) {
		if(variantsIndexKeep.contains(variant_counter++)){
			outputvcf << line
			outputvcf << '\n'
		}
	}	
}
outputvcf.flush()


class VariantLight {
	def variantIndex
	def formatIndex
	def density

	def initialize(String line) {
		def format
		def genotypes = new ArrayList<String>()
		def counter = 0
		line.split('\\t').each {	
			if(counter==formatIndex) {
				format=it.split(':').toList()
			}
			else if(counter>formatIndex) {
				genotypes.add(it)
			}
			counter++
		}

		density=0.0
		def GT_index = format.findIndexOf(){ it == 'GT' }
		def DP_index = format.findIndexOf(){ it == 'DP' }
		genotypes.eachWithIndex { genotype, index ->
			String GT = genotype.split(':')[GT_index]
			if(!GT.equals('./.') && DP_index>=0){
				def DP = genotype.split(':')[DP_index]
				if(DP.isInteger()){
					int DP_value = DP as Integer
					if(DP_value>=2)
						density++
				}
			}
			else if(!GT.equals('./.')){
				density++
			}
		}	
		density /= genotypes.size()
		genotypes.clear()
	}

}
