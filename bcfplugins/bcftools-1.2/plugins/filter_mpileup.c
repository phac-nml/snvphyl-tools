#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <htslib/vcf.h>
#include <htslib/regidx.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"

/* 
    This short description is used to generate the output of `bcftools annotate -l`.
*/
const char *about(void)
{
    return 
        "A plugin to soft filter position if they have a SNP or have low coverage.\n";
}

const char *usage(void)
{
  return "";

}


int min_dp;
bcf_hdr_t *in_hdr;
bcf_hdr_t *out_hdr;
int flag_coverage;
int flag_mpileup;
int verbose;


/* 
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)

{
 
    min_dp=-1;
    flag_coverage=-1;
    flag_mpileup=-1;
    verbose=-1;

    static struct option loptions[] =
    {
        {"verbose",1,0,'v'},
        {"dp",1,0,'d'},
        {0,0,0,0}
    };
    char c, *tmp;
    while ((c = getopt_long(argc, argv, "d:v",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': verbose = 1; break; 
            case 'd': 
                min_dp = strtol(optarg,&tmp,10); 
                if (*tmp) error("Unexpected argument to --dp: %s\n", optarg); break; 
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    in_hdr  = in;
    out_hdr  = out;


    //add this bogus flag otherwise for some reason the 'DP4' flag is the same id as first filtered 
    bcf_hdr_append(out_hdr,"##FILTER=<ID=DP44,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");
    int dp4 = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "DP44");
    dp4=-1;
    //needs to golf down into a simple plugin and reported as a bug. Phil Feb 16,2015 


    char info[100];
    sprintf(info, "##FILTER=<ID=filtered-coverage,Description=\"Set true if DP is less then %d\">", min_dp);
    bcf_hdr_append(out_hdr, info);

    //get flag index id so we can mark as such
    flag_coverage = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"filtered-coverage");

  
          
    bcf_hdr_append(out_hdr, "##FILTER=<ID=filtered-mpileup,Description=\"Set if not true: SNP not found in Freebayes\">");
    
    //get flag index id so we can mark as such
    flag_mpileup = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "filtered-mpileup");

    return 0;
}



int pl_type = BCF_HT_INT;
int *buf = NULL;
int nbuf = 0;   // NB: number of elements, not bytes


bcf1_t *process(bcf1_t *rec)
{

      
  if (bcf_get_info_values(in_hdr,rec,"DP",(void**)&buf,&nbuf,pl_type)) {
    int DP = *buf;

    //check to see if we have or over the minimum depth of coverage
    if (DP >= min_dp) {
      //since we are over the coverage, mark as PASS
      bcf_add_filter(out_hdr,rec,0);
    }
    else {
      //since we are under the coverage, mark as 'filtered-coverage'
      bcf_add_filter(out_hdr,rec,flag_coverage);
    }
  }

  
  int type = bcf_get_variant_types(rec);
  //if we have SNP at this point, implies not found in freebayes so should mark as filtered-mileup
  if ( type&VCF_SNP ){
    bcf_add_filter(out_hdr,rec,flag_mpileup); 
  }
  

  return rec;
}





/*
    Clean up.
*/
void destroy(void)
{

}


