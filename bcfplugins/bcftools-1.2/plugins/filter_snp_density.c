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
void print_density_region(int start, int end, const char *chromosome);
void destroy(void);

//define the script variables
char *filename = NULL;
int density_threshold = 10;
int flag_density = -1;
int verbose = -1;
int eof = 0;
bcf_hdr_t *in_hdr = NULL;
bcf_hdr_t *out_hdr = NULL;
int in_density_region = 0;
char *regionFilepath=NULL;
FILE *regionsFile;
//for look ahead functionality
htsFile *htsAhead=NULL; 
bcf_hdr_t *htsAheadHdr = NULL;
bcf1_t *bcf_ahead = NULL;

//track the density regions
int density_start_position = 0;
int density_end_position = 0;

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
       {"region_file",1,0,'r'},	   
       {"threshold",1,0,'t'},
       {0,0,0,0}
   };
   char c, *tmp;
   while((c = getopt_long(argc, argv, "f:r:t", loptions, NULL)) >= 0)
   {
       switch (c)
      {  
         case 'f':
	       filename = optarg; 
	       break;
         case 'r':
               regionFilepath = optarg;
               break;
         case 't':	 
           density_threshold = strtol(optarg, &tmp, 10);
           break;
           if (*tmp) error("Unexpected argument to -t: %s\n", optarg); break;
         case '?':
         default: error("%s", usage()); break;
      }
   }
   
   //perform som validation of input parameters
   if(density_threshold <= 0){
	   density_threshold = 10;
	   printf("Negative density_threshold value found, setting to default value of %d.", density_threshold);
   }
   
   in_hdr = in;
   out_hdr = out;

   //input the bcf file in order to "look ahead" when iterating over the records
   htsAhead = hts_open(filename, "r");
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
      printf("Unable to read look ahead bcf record for analysis OR we have reached the end of the bcf file.\n");
      eof = 1;    	   
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
            printf("Unable to add density flag to vcf record.");       
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
    
    //check the absolute difference between snv's to determine if it is within a dense region 
    if(((next_position - current_position) <= density_threshold) && ((current_position - next_position) <= density_threshold) && (current->rid == next->rid) && (!eof)){
        //track the start and end positions for the density regions appropriately:
        //if not previously in a denisty region, mark the start of density region:
        if(!in_density_region){
           density_start_position = current_position;
        }
    	in_density_region = 1;
	return 1;
    }

    if(in_density_region){
        density_end_position = current_position;
        //print the newly discovered density region to a tab seperated file:
        print_density_region(density_start_position, density_end_position, in_hdr->id[BCF_DT_CTG][current->rid].key);
        //if we have come to the end of a density region, we still must ensure the last snv
        //in the region is marked as filtered-density.
        in_density_region = 0;
        return 1;
    }
    in_density_region = 0;
    return 0;
}

/**
* Prints the high density regions to the output regions file in a tab seperated format.
* @start The start position of the high density region
* @end The end position of the high density region
* @chromosome The name of the chromosome containing the region
*/
void print_density_region(int start, int end, const char *chromosome){
   //we need to add 1 position as the output appears to start counting from 0, but this
   //will create an off by one error in the tab seperated regions file.
   start = start +1;
   end = end +1;
   regionsFile = fopen(regionFilepath, "a");
   fprintf(regionsFile, "%s\t%d\t%d\n", chromosome, start, end);
   fclose(regionsFile);
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
   hts_close(htsAhead);
}

