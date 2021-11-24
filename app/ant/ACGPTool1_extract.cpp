#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUMSTTCOLS 20

//#define DEBUG


void doc(char prog[]);
void openFiles(char *argv[], FILE **inP, FILE **outP);
int countSubpops(FILE *in);
bool getGenData(FILE *in, int numPops, float data[][NUMSTTCOLS]); 
void processData(FILE *out, int numSubpops, float data[][NUMSTTCOLS], int colsToProcess[], int numColsToProcess, bool average); 
float avg(float data[][NUMSTTCOLS], int col, int numRows);

float(*colFs[NUMSTTCOLS])(float [][NUMSTTCOLS],int,int) =
        {0,0,avg,avg,avg,0,0,avg,avg,avg,
         avg,0,avg,avg,0,0,avg,avg,avg,avg};

char *colTitles[NUMSTTCOLS] = {
        "Gen", "Subpop", "FitMeanOfGen", "FitBestOfGen", "FitWorstOfGen",
        "SizMeanOfGen", "DepMeanOfGen", "SizBestOfGen", "DepBestOfGen", "SizWorstOfGen",
        "DepWirstOfGen", "FitMeanOfRun", "FitBestOfRun", "FitWorstOfRun", "SizMeanOfRun",
        "DepMeanOfRun", "SizBestOfRun", "DepBestOfRun", "SizWorstOfRun", "DepWorstOfRun"};

void doc(char prog[]) {
   printf("\n\n\
   Cezary Z Janikow\n\
   8/21/08\n\
   Tool for command-line extract of columns from stt file\n\n");
   printf("Invocation: > %s {$In.stt|-} {$Out.stte|-} {avg|best} {col}+\n\n",prog);
   printf("\
   Extract columns (stt columns counted from 1, can ask 3rd or higher):\n\
   - input taken from $In.stt or stdin if '-'\n\
   - output written to $Out.stte or stdout if '-'\n\n\
     - first write title for each column\n\
   - columns are extracted and order preserved\n\
   - {col}+ lists at least one column, will be written in the order listed\n\
     - gen# is always written\n\
   - if mult subpops\n\
     - only the summary '0' row is extracted and written \n\
     - if 'avg' then columns producing best-of/worst-of are computed as averages across the subpopulation\n\
       and thus using the subpopulations as independent runs being averaged\n");
}

int main(int argc, char *argv[]) {
  if (argc <5 || argc>(4+NUMSTTCOLS) || !(!strcmp(argv[3],"best") || !strcmp(argv[3],"avg"))) {
    doc(argv[0]);
    exit(1);
  }

  int numSubpops;
  FILE *in;
  FILE *out;

  openFiles(argv, &out, &in);
  numSubpops=countSubpops(in);
  bool average;
  int numColsToProcess=argc-4;
  
  if (!strcmp(argv[3],"best"))
    average=false;
  else
    average=true;
  int *colsToProcess = new int[numColsToProcess];
  int col;
  for (int i=0; i<numColsToProcess; i++) {
    col=atoi(argv[4+i]);
    if (col<3 || col>NUMSTTCOLS) {
      printf("Wrong column number [%d]\n",col);
      exit(1);
    }
    else {
      colsToProcess[i]=col;
    }
  }

#ifdef DEBUG
  printf("numSubpops=%d\n",numSubpops);
#endif  
  
  if (numSubpops>1 && average) {
    fprintf(out,"Average of %d independent settings\n",numSubpops);
  }
  fprintf(out,"%s ","Gen");
  for (int i=0; i<numColsToProcess; i++) {
    fprintf(out,"%s ",colTitles[colsToProcess[i]-1]);
  }
  fprintf(out,"\n");

  typedef float oneRowData_t[NUMSTTCOLS];
  oneRowData_t *data=new oneRowData_t[numSubpops+1]; // +1 for the "0" row if there
  bool moreData;
  moreData=getGenData(in,numSubpops,data);
  while (moreData) {
    processData(out, numSubpops, data, colsToProcess, numColsToProcess, average);
    moreData=getGenData(in,numSubpops,data);
  }
  return 0;
}

/* open files, names are in argv[1] for in and argv[2] for out */
void openFiles(char *argv[], FILE **outP, FILE **inP) {
  if (!strcmp(argv[1],"-")) {
    char line[256];
    fgets(line,256,stdin);
    FILE *fp=tmpfile();
    while (!feof(stdin)) {
      fputs(line,fp);
      fgets(line,256,stdin);
    }
    rewind(fp);
    *inP=fp;
  }
  else {
    *inP=fopen(argv[1],"r");
  }
  if (!strcmp(argv[2],"-")) {
    *outP=stdout;
  }
  else {
    *outP=fopen(argv[2],"w");
  }
  if (!(*inP) || !(*outP)) {
    printf("Could not open %s or %s\n",argv[1],argv[2]);
    exit(1);
  }
}

/* peek into in file and return the number of subpopulations, revind in file */
int countSubpops(FILE *in) {
  int numSubpops=1;
  int gen1, gen2, pop;
  fscanf(in,"%d%d%*[^\n]",&gen1,&pop);
  fscanf(in,"%d%d%*[^\n]",&gen2,&pop);
  while (!feof(in) && gen1==gen2 && pop!=0) {
    fscanf(in,"%d%d%*[^\n]",&gen2,&pop);
    numSubpops++;
  }
  rewind(in);
  return numSubpops;
}

/* process one generation in file in opened and properly positioned */
/* this file has numSubpops */
/* write into data , includiing summary "0" row if numSubpops>1 */
/* return true if new complete data was read */
bool getGenData(FILE *in, int numSubpops, float data[][NUMSTTCOLS]) {
  for (int i=0; i<numSubpops; i++) {
    for (int j=0; j<NUMSTTCOLS; j++) {
      fscanf(in,"%f",&data[i][j]); 
    }
  }
  if (numSubpops>1) {
    for (int j=0; j<NUMSTTCOLS; j++) {
      fscanf(in,"%f",&data[numSubpops][j]); 
    }
  }
  return !feof(in);
}

/* process data from data write to out */
/* write gen# followed by values of the columns in colsToProcess */
/* process according to average */
void processData(FILE *out, int numSubpops, float data[][NUMSTTCOLS], int colsToProcess[], int numColsToProcess, bool average) {
  int col;
  fprintf(out,"%-3d ",(int)data[0][0]);
  if (numSubpops==1) { // write from data[0]
    for (int i=0; i<numColsToProcess; i++) {
      col=colsToProcess[i]-1;
      fprintf(out,"%f ",data[0][col]);
    }
    fprintf(out,"\n");
  }
  else if (!average) { // write from data[numSubpops]
    for (int i=0; i<numColsToProcess; i++) {
      col=colsToProcess[i]-1;
      fprintf(out,"%.3f ",data[numSubpops][col]);
    }
    fprintf(out,"\n");
  }
  else { // write from data[*] calling colsFs()
    for (int i=0; i<numColsToProcess; i++) {
      col=colsToProcess[i]-1;
      if (colFs[col]==0) {
        fprintf(out,"%.3f ",data[numSubpops][col]);
      }
      else {
        fprintf(out,"%.3f ",colFs[col](data,col,numSubpops));
      }
    }
    fprintf(out,"\n");
  }
}

/* return avg on column col in the array data of numRows rows */
float avg(float data[][NUMSTTCOLS], int col, int numRows) {
  float sum=0;
  for (int i=0; i<numRows; i++) {
    sum+=data[i][col];
  }
  return sum/numRows;
}

