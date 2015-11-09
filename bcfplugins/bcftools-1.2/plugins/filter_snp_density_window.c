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
#include <sys/queue.h>
#include "bcftools.h"

const char *usage(void);
bcf1_t *process(bcf1_t *rec);
void destroy(void);
struct entry *initializeStruct(bcf1_t *record);
void update_window();
void update_density();
void refresh_density();
void output_density_region();
void print_density_region(int start, int end, const char *chromosome);



//define the script variables
char *filename = NULL;
int density_threshold = 10;
int window_size = 100;
//output regions file
char *regionFilepath=NULL;
FILE *regionsFile;

//htslib requirements:
bcf_hdr_t *in_hdr = NULL;
bcf_hdr_t *out_hdr = NULL;

//tracking variables:
int queue_size = 0;
int density_flag = 0;
int density_start_position = 0;
int density_end_position = 0;

//define the tail queue:
TAILQ_HEAD(tailhead, entry) head;
struct tailhead *headp;
struct entry{
   TAILQ_ENTRY(entry) entries;
   char *contig;
   int position;   
};

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
       {"window_size", 1, 0, 'w'},
       {0,0,0,0}
   };
   char c, *tmp;
   while((c = getopt_long(argc, argv, "f:r:w:t", loptions, NULL)) >= 0)
   {
       switch (c)
      {  
         case 'f':
	       filename = optarg; 
	       break;
         case 'r':
               regionFilepath = optarg;
               break;
         case 'w':
               window_size = strtol(optarg, &tmp, 10);
         case 't':	 
           density_threshold = strtol(optarg, &tmp, 10);
           break;
           if (*tmp) error("Unexpected argument to -t: %s\n", optarg); break;
         case '?':
         default: error("%s", usage()); break;
      }
   }
   
   //perform some validation of input parameters
   if(density_threshold <= 0){
	   density_threshold = 10;
	   printf("Negative density_threshold value found, setting to default value of %d.", density_threshold);
   }
  
   if(window_size <= 0){
      window_size = 100;
      printf("Invalid window_size value found, setting to the default value of %d.", window_size);
   }

   in_hdr = in;

   //initialize the queue:
   TAILQ_INIT(&head);

   return 0;
}

bcf1_t *process(bcf1_t *hts){

   //create the new struct for this variant:
   struct entry *tempStruct = initializeStruct(hts);

   //for the first record, add the record to the head of the queue:
   if(queue_size == 0){
      TAILQ_INSERT_HEAD(&head, tempStruct, entries);
      queue_size++;
   }
   //otherwise, add to the tail of the queue and conduct density analysis:
   else{
      TAILQ_INSERT_TAIL(&head, tempStruct, entries);
      queue_size++;
      update_density(); 
   }   
  
   //return the bcf record
   return hts;
}

/*
* Method to allocate memory for the variant struct and assign the correct
*position and contig values. 
*/
struct entry *initializeStruct(bcf1_t *record){
   struct entry *tempStruct = malloc(sizeof(struct entry));
   tempStruct->contig = in_hdr->id[BCF_DT_CTG][record->rid].key;
   tempStruct->position = record->pos;
   return tempStruct;
}

/*
* Method to update the density calculations for the new variant entry.
*/
void update_density(){
   //update queue to ensure all of the entries are within the window size:
   update_window();
   //check the density for the updated window:
   refresh_density();
}

/*
* Deletes entries from the head of the list that are no longer
* within the window of interest.
*/
void update_window(){
   //update the window when a contig change occurs and remove all but tail:
   if(strcmp(TAILQ_LAST(&head, tailhead)->contig, head.tqh_first->contig)!=0){
      //if currently in density region, output:
      output_density_region();
      while(strcmp(head.tqh_first->contig, TAILQ_LAST(&head, tailhead)->contig)!=0){
          TAILQ_REMOVE(&head, head.tqh_first, entries);
          queue_size--;
      }
   }
   else{
      //if no contig change, do a regular window update:
      while((TAILQ_LAST(&head, tailhead)->position - head.tqh_first->position) > window_size){
         output_density_region();
         TAILQ_REMOVE(&head, head.tqh_first, entries);
         queue_size--;
      }   
   }
}

/*
*
*/
void refresh_density(){
   if(queue_size > density_threshold){
      density_flag = 1;
      density_start_position = head.tqh_first->position;
      density_end_position = TAILQ_LAST(&head, tailhead)->position;
   }
   else{
      density_end_position = TAILQ_LAST(&head, tailhead)->position; 
   }
}


/*
* Outputs a new density region to the output file, if applicable. Resets the start and end positions for a new analysis. 
*/
void output_density_region(){
   //check if the density flag is set:
   if(density_flag){
      //there is something to output:
      print_density_region(density_start_position, density_end_position, head.tqh_first->contig);
      density_flag = 0;
   }
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
   
}

