#
# calculate_stats.sh
# Bash script for calculating per variant statistics of vcf files
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

#$1 input
#$2 output

cat $1 |\
awk -F"\t" '{line=$0} BEGIN { \
        print "CHR\tPOS\tID\tA\ta\tAA\tAa\taa\taa_e\tAa_e\taa_e\tMissing\tMissPercentage\tAF_A\tAF_a\tMAF\tHo\tHe\tPIC\tF\tchi_sq_2df" \
    } !/^#/ { AA=gsub(/\t(0\/0|0\|0)/,"") ;
              Aa=gsub(/\t(0\/1|1\/0|0\|1|1\|0)/,"") ;
              aa=gsub(/\t(1\/1|1\|1)/,"") ;
              mis=gsub(/\t(.\/.|.\|.)/,"") ;
              if (AA+Aa+aa==0) {
                AA_e=0 ;
                Aa_e=0 ;
                aa_e=0 ;
                hetperc="NA" ;
                misperc=1 ;
                AF_A="NA" ;
                AF_a="NA" ;
                MAF="NA" ;
                Ho="NA" ;
                He="NA" ;
                PIC="NA" ;
                Fst="NA" ;
                chi_sq_1df="NA" ;
                chi_sq_2df="NA" ;
              }
              else {
                misperc=mis/(AA+Aa+aa+mis) ;
                AF_A=(AA+0.5*Aa)/(AA+Aa+aa) ;
                AF_a=(aa+0.5*Aa)/(AA+Aa+aa) ;
                AA_e=(AA+Aa+aa)*(AF_A^2) ;
                aa_e=(AA+Aa+aa)*(AF_a^2) ;
                Aa_e=(AA+Aa+aa)*2.0*AF_A*AF_a ;
                MAF=AF_a ;
                if (AF_a>AF_A) MAF=AF_A ;
                Ho=Aa/(AA+Aa+aa) ;
                He=1.0-MAF^2-(1.0-MAF)^2 ;
                PIC=He-2.0*(MAF^2)*((1.0-MAF)^2) ;
                if (He==0) {Fst="NA" ;}
                else {Fst=(He-Ho)/He ;}
                if (AA_e!=0 && Aa_e!=0 && aa_e!=0){
                  chi_sq_2df=(((AA_e-AA)^2)/AA_e)+(((Aa_e-Aa)^2)/Aa_e)+(((aa_e-aa)^2)/aa_e);
                }
                else{
                  chi_sq_2df="NA";
                }
              }
              print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t" \
              AA "\t" \
              Aa "\t" \
              aa "\t" \
              AA_e "\t" \
              Aa_e "\t" \
              aa_e "\t" \
              mis "\t" \
              misperc "\t" \
              AF_A "\t" \
              AF_a "\t" \
              MAF "\t" \
              Ho "\t" \
              He "\t" \
              PIC "\t" \
              Fst "\t" \
              chi_sq_2df \
    }' > $1.stats.tmp ;

unset R_HOME ;
Rscript --vanilla <(echo -e "\
  data<-read.table(file='$1.stats.tmp',\
    colClasses=c('NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',\
                 'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',NA),header=T,fill=T)
  res.p2<-pchisq(data\$chi_sq_2df, df=2, lower.tail=FALSE)
  write.table(res.p2,file='$1.stats.tmp2',sep='\t',col.names=F,row.names=F)") ;

sed -i '1iP2_value' $1.stats.tmp2 ;

paste $1.stats.tmp $1.stats.tmp2 > $2 ;

rm -f $1.stats.tmp $1.stats.tmp2 ;
