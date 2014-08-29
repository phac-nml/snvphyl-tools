#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include "config.h"
#include <wordexp.h>

/* 
    This short description is used to generate the output of `bcftools annotate -l`.
*/
const char *about(void)
{
    return 
        "A minimal plugin which counts number of SNPs, Indels, and\n"
        "total number of sites.\n";
}

int min_dp =-1;
bcf_hdr_t *in_hdr = NULL;
bcf_hdr_t *out_hdr = NULL;
int flag_coverage=-1;
int flag_mpileup=-1;

/* 
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{
  char *temp = config_get_string(opts,"dp");

  wordexp_t wexp;
  wordexp(temp, &wexp, 0);

  free(temp);
  min_dp = atoi(wexp.we_wordv[0]);
  printf("Received minimum coverage of '%d'\n",min_dp);
  
  in_hdr  = in;
  out_hdr  = out;

  char info[80];
  sprintf(info, "##FILTER=<ID=filtered-coverage,Description=\"Set true if DP < %d\">", min_dp);
  bcf_hdr_append(out_hdr, info);

  //get flag index id so we can mark as such
  flag_coverage = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"filtered-coverage");


  bcf_hdr_append(out_hdr, "##FILTER=<ID=filtered-mpileup,Description=\"Set if not true: SNP not found in Freebayes\">");

  //get flag index id so we can mark as such
  flag_mpileup = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"filtered-mpileup");
  
  return 0;
}



int pl_type = BCF_HT_INT;
int *buf = NULL;
int nbuf = 0;   // NB: number of elements, not bytes

/*
    Called for each VCF record after all standard annotation things are finished.
    Return 0 on success, 1 to suppress the line from printing, -1 on critical errors.
*/
int process(bcf1_t *rec)
{

  if (bcf_get_info_values(in_hdr,rec,"DP",(void**)&buf,&nbuf,pl_type)) {
    int DP = *buf;

    //check to see if we have or over the minimum depth of coverage
    if (DP >= min_dp) {
      //printf("DP is '%d'\n",DP);
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
  
  
  return 0;
}





/*
    Clean up.
*/
void destroy(void)
{

}


