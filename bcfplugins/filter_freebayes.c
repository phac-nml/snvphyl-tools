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
      "A plugin which filters on freebayes based on depth of coverage (INFO/DP), mean mapping quality (MQM) and \n"
        "Percentage of alternative alleles from total DP of coverage (INFO/AO / INFO/DP ) * 100 \n";
}

int dp_type = BCF_HT_INT;
int mqm_type = BCF_HT_REAL;
int type_type = BCF_HT_STR;
int min_dp =-1;
double min_ao =-1;
float min_mqm = -1;
bcf_hdr_t *in_hdr = NULL;
bcf_hdr_t *out_hdr = NULL;
int flag_coverage=-1;

/* 
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{

  wordexp_t wexp;
  
  char *temp = config_get_string(opts,"dp");


  wordexp(temp, &wexp, 0);
  min_dp = atoi(wexp.we_wordv[0]);

  temp = config_get_string(opts,"mqm");
  wordexp(temp, &wexp, 0);


  min_mqm = atoi(wexp.we_wordv[0]);
  temp = config_get_string(opts,"ao");


  wordexp(temp, &wexp, 0);
  min_ao = atoi(wexp.we_wordv[0]);
  

  free(temp);
  
  in_hdr  = in;
  out_hdr  = out;

  char info[80];
  sprintf(info, "##FILTER=<ID=filtered-coverage,Description=\"Set true if DP < %d\">", min_dp);
  bcf_hdr_append(out_hdr, info);
    
  //get flag index id so we can mark as such
  flag_coverage = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"filtered-coverage");

  
  return 0;
}




/*
    Called for each VCF record after all standard annotation things are finished.
    Return 0 on success, 1 to suppress the line from printing, -1 on critical errors.
*/
int process(bcf1_t *rec)
{

  int *buf = NULL;
  int nbuf = 0;   // NB: number of elements, not bytes
  int DP =-1;
  if (bcf_get_info_values(in_hdr,rec,"DP",(void**)&buf,&nbuf,dp_type)) {
    DP = *buf;

    //check to see if we have or over the minimum depth of coverage
    //if not we are softclipping out of file
    if (DP < min_dp) {
      bcf_add_filter(out_hdr,rec,flag_coverage);
    }
  }
  else {
    printf ("A record does not contain DP flag in INFO column\n");
    return -1;
  }


  float *buf2 = NULL;
  int nbuf2 = 0;   // NB: number of elements, not bytes
  if (bcf_get_info_values(in_hdr,rec,"MQM",(void**)&buf2,&nbuf2,mqm_type)) {
    float MQM = *buf2;

    //check to see if we have good mean mapping quality of observed for alternate alleles
    //if not we are hardclipping out of file
    if (MQM < min_mqm) {
      return 1;
    }
  }
  else {
    printf ("A record does not contain MQM flag in INFO column\n");
    return -1;
  }


  //get INFO/AO and see if over % of total reads

  if (bcf_get_info_values(in_hdr,rec,"AO",(void**)&buf,&nbuf,dp_type)) {
    int AO = *buf;

    //determine percentage of reads that match the alternative allele
    double test = (double)AO/DP * 100;
    if (test < min_ao ) {
      return 1; 
    }
  }
  else {
    printf ("A record does not contain AO flag in INFO column\n");
    return -1;
  }



/* turn off complex/mnp region filtering completely for now. 
 * Want to get the same results as Aaron version
  char *buf3 = NULL;
  int nbuf3 = 0;   // NB: number of elements, not bytes
  if (bcf_get_info_values(in_hdr,rec,"TYPE",(void**)&buf3,&nbuf3,type_type)) {
    char *type = buf3;

    //we are throwing out complex or mnp regions
    if (strcmp(type,"complex") == 0 || strcmp(type,"mnp") == 0) {
      return 1;
    }
  }
  else {
    printf ("A record does not contain TYPE flag in INFO column\n");
    return -1;
  }
*/
  //add that we passed all filtering!
  bcf_add_filter(out_hdr,rec,0);
  return 0;
}





/*
    Clean up.
*/
void destroy(void)
{

}

