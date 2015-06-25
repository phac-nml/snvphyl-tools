#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>
#include <htslib/vcf.h>
#include <htslib/regidx.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"

const char *usage(void);
int mark_density_snp(bcf1_t *bcf_current, bcf_hdr_t *header, int flag_density);
int check_density(bcf1_t *current, bcf1_t *next, int density_threshold);
bcf1_t *process(bcf1_t *rec);
void destroy(void);

//define the script variables
char *filename = NULL;
int density_threshold = 10;
int flag_density = -1;
int verbose = -1;
bcf_hdr_t *in_hdr = NULL;
bcf_hdr_t *out_hdr = NULL;

//for look ahead functionality
htsFile *htsAhead=NULL; 
bcf_hdr_t *htsAheadHdr = NULL;
bcf1_t *bcf_ahead = NULL;

const char *about(void)
{
    return 
      "A plugin which filters on freebayes for SNP's deemed to be within high density regions of the genome.\n"; 
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
   static struct option loptions[] =
   {
	   {"filename",1,0,'f'},		   
       {"threshold",1,0,'t'},
       {0,0,0,0}
   };
   char c, *tmp;
   while((c = getopt_long(argc, argv, "f:t:", loptions, NULL)) >= 0)
   {
       switch (c)
      {  
         case 'f':
	       filename = optarg; 
	       break;
         case 't':	 
           density_threshold = strtol(optarg, &tmp, 10); break;
           if (*tmp) error("Unexpected argument to -t: %s\n", optarg); break;
         case '?':
         default: error("%s", usage()); break;
      }
   }
         
   in_hdr = in;
   out_hdr = out;

   //input the bcf file in order to "look ahead" when iterating over the records
   htsAhead = hts_open(filename, "ru");
   htsAheadHdr = bcf_hdr_read(htsAhead);
   bcf_ahead = bcf_init();
   //skip the first record so that we are always looking ahead when the 
   //process function is called
   if(bcf_read(htsAhead, htsAheadHdr, bcf_ahead) == -1){
       printf("Unable to read look ahead bcf record for analysis.\n");    	   
   }
   
   //generate the new filter value
   char info[100];
   sprintf(info, "##FILTER=<ID=filtered-density,Description=\"Set true if spacing is < %d bp\">", density_threshold);
   bcf_hdr_append(out_hdr, info);
   flag_density = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"filtered-density");
   
   return 0;
}

bcf1_t *process(bcf1_t *hts){
   if(bcf_read(htsAhead, htsAheadHdr, bcf_ahead) == -1){
      printf("Unable to read look ahead bcf record for analysis.\n");    	   
   }
   
   if(check_density(hts, bcf_ahead, density_threshold)){
       mark_density_snp(hts, out_hdr, flag_density);
   }
   //return the bcf record
   return hts;
}

/*
 * Method to mark the input bcf records as filtered-density
 * @return int 1 for success, 0 for error
 */
int mark_density_snp(bcf1_t *bcf_current, bcf_hdr_t *header, int flag_density){
	if(bcf_add_filter(header, bcf_current, flag_density)){
		return 1;
	}
	else{
		return 0;
	}
}

/*
 * Method to determine if a particular position is within a high density
 * region, based on the density threshold value.
 * @return int 1 if it is a high density position, 0 if not
 */
int check_density(bcf1_t *current, bcf1_t *next, int density_threshold){
	int current_position = current->pos;
	int next_position = next->pos;

	if(current->rid != next->rid){
	    if((next_position - current_position) <= density_threshold){
		    return 1;
	    }
	}
	return 0;
}

const char *usage(void)
{
  return "";
}

/*
    Clean up.
*/
void destroy(void)
{

}

