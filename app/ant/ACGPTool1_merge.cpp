#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUMSTTCOLS 20

//#define DEBUG


void doc(char prog[]);
void openFiles(int numFiles, char *argv[], FILE *fps[], FILE **outP);
int countSubpops(int numFiles, FILE *fps[], int popsPerFile[]);
bool processFile(FILE *fps, int popsPerFile, int curSubpop, float data[][NUMSTTCOLS]);
float min(float data[][NUMSTTCOLS], int col, int numRows);
float max(float data[][NUMSTTCOLS], int col, int numRows);
float avg(float data[][NUMSTTCOLS], int col, int numRows);
float best4(float data[][NUMSTTCOLS], int col, int numRows);
float worst4(float data[][NUMSTTCOLS], int col, int numRows);

float(*colFs[NUMSTTCOLS])(float [][NUMSTTCOLS],int,int) = 
	{0,0,avg,max,min,avg,avg,best4,best4,worst4,
         worst4,avg,max,min,avg,avg,best4,best4,worst4,worst4};

void doc(char prog[]) {
   printf("\n\n\
   Cezary Z Janikow\n\
   8/21/08\n\
   Tool for command-line merge of stt files\n\n");
   printf(" Invocation:\n > %s {$Out.stt|-} {List of stt files}+\n\n",prog);
   printf("\
   Merging multiple stt files:\n\
   - files are assumed to be produced on the same application and the\n\
     same parameters except for random seed\n\
   - each file can have single or multiple populations\n\
   - if a file terminates sooner than longest file, it will reproduce its\n\
     last result for the remaining generations (the assumption is runs for\n\
     the same number of generations so shorter file implies termination)\n\n\
   Output:\n\
   - $Out.stt has format identical to what lilgp would have produced\n\
     if all runs were done in multiple subpopulation\n\n\
   - Print to stdout if '-' \n\n");
}

int main(int argc, char *argv[]) {
  if (argc<3) {
    doc(argv[0]);
    exit(1);
  }

  int numFiles=argc-2;
  int totSubpops;
  FILE **fps = new FILE*[numFiles];
  FILE *out;
  int *popsPerFile=new int[numFiles];

  openFiles(numFiles, argv, fps, &out);
   
  totSubpops=countSubpops(numFiles, fps, popsPerFile);

#ifdef DEBUG
  printf("totSubpops=%d\n",totSubpops);
  for (int i=0; i<numFiles; i++)
    printf("File%d=%d\n",i,popsPerFile[i]);
#endif  

  typedef float oneRowData_t[NUMSTTCOLS];
  oneRowData_t *data=new oneRowData_t[totSubpops];
  int curGen=0;
  int curSubpop=0;
  bool moreInput=true;
  while (moreInput) {
    curSubpop=0;
    moreInput=false;
    for (int i=0; i<numFiles; i++) {
      moreInput = processFile(fps[i],popsPerFile[i],curSubpop,data) || moreInput;
      curSubpop+=popsPerFile[i];
    }
#ifdef DEBUG
    for (int i=0; i<totSubpops && moreInput; i++) {
      for (int j=0; j<NUMSTTCOLS; j++) {
        printf("%f ",data[i][j]);
      }
      printf("\n");
    }
#endif
    for (int i=0; i<totSubpops && moreInput; i++) {
      fprintf(out,"%-3d %-2d ",curGen,i+1);
      for (int j=2; j<NUMSTTCOLS; j++) {
        fprintf(out,"%6.3f ",data[i][j]);
      }
      fprintf(out,"\n");
    }
    if (totSubpops>1 && moreInput) {
      fprintf(out,"%-3d %-2d ",curGen,0);
      for (int i=2; i<NUMSTTCOLS; i++)
        fprintf(out,"%6.3f ",colFs[i](data,i,totSubpops));
      fprintf(out,"\n");
    }
    curGen++;
  }
  return 0;
}

/* open numFiles, names are in argv[2]+, put file pointers in fps */
/* also open outP for output */
void openFiles(int numFiles, char *argv[], FILE *fps[], FILE **outP) {
  for (int i=0; i<numFiles; i++) {
    fps[i]=fopen(argv[i+2],"r");
    if (!fps[i]) {
      printf("Could not open %s\n",argv[i+2]);
      exit(1);
    }
  }
  if (!strcmp(argv[1],"-")) {
    *outP=stdout;
  }
  else {
    *outP=fopen(argv[1],"w");
  }
  if (!(*outP)) {
    printf("Could not open %s\n",argv[1]);
    exit(1);
  }
}

int countSubpops(int numFiles, FILE *fps[], int popsPerFile[]) {
  int totSubpops=0;
   int gen1, gen2, pop;
   int subPop;
  for (int i=0; i<numFiles; i++) {
    fscanf(fps[i],"%d%d%*[^\n]",&gen1,&pop);
    fscanf(fps[i],"%d%d%*[^\n]",&gen2,&pop);
    if (gen1==gen2) { // multiple subpops
      subPop=1;
      while (!feof(fps[i]) && gen2==0 && pop!=0) {
        fscanf(fps[i],"%d%d%*[^\n]",&gen2,&pop);
        subPop++;
      }
    }
    else { // one population
      subPop=1;
    }
    popsPerFile[i]=subPop;
    totSubpops+=subPop;
    rewind(fps[i]);
  }
  return totSubpops;
}

/* process one generation in file fps opened and properly positioned */
/* this file has popsPerFile generations but it can end early (no more generations) */
/* and if so, use previous generation (data from data) */
/* write into data starting at row curSubpop */
/* return true if new complete date was read */
bool processFile(FILE *fps, int popsPerFile, int curSubpop, float data[][NUMSTTCOLS]) {
  if (feof(fps)) {
    return false; // previous eof
  }
  for (int i=0; i<popsPerFile; i++) {
    for (int j=0; j<NUMSTTCOLS; j++) {
       if (fscanf(fps,"%f",&data[curSubpop+i][j])<1) { // current eof
         return false;
       }
    }
  }
  if (popsPerFile>1) {
    for (int i=0; i<NUMSTTCOLS; i++) {
      fscanf(fps,"%*f");	// skip summary row for mult subpopulations
    }
  }
  return true;
}

/* return min on column col in the array data of numRows rows */
float min(float data[][NUMSTTCOLS], int col, int numRows) {
  float minData=data[0][col];
  for (int i=1; i<numRows; i++) {
    if (data[i][col]<minData) {
      minData=data[i][col];
    }
  }
  return minData;
}

/* return max on column col in the array data of numRows rows */
float max(float data[][NUMSTTCOLS], int col, int numRows) {
  float maxData=data[0][col];
  for (int i=1; i<numRows; i++) {
    if (data[i][col]>maxData) {
      maxData=data[i][col];
    }
  }
  return maxData;
}

/* return avg on column col in the array data of numRows rows */
float avg(float data[][NUMSTTCOLS], int col, int numRows) {
  float sum=0;
  for (int i=0; i<numRows; i++) {
    sum+=data[i][col];
  }
  return sum/numRows;
}

/* return column col  value from element with the best column 4 in data in numRows */
float best4(float data[][NUMSTTCOLS], int col, int numRows) {
  int max4I=0;
  float max4V=data[max4I][4-1];
  for (int i=1; i<numRows; i++) {
    if (data[i][4-1]>max4V) {
      max4V=data[i][4-1];
      max4I=i;
    }
  }
  return data[max4I][col];
}

/* return column col value from element with the worst column 4 in data in numRows */
float worst4(float data[][NUMSTTCOLS], int col, int numRows) {
  int min4I=0;
  float min4V=data[min4I][4-1];
  for (int i=1; i<numRows; i++) {
    if (data[i][4-1]<min4V) {
      min4V=data[i][4-1];
      min4I=i;
    }
  }
  return data[min4I][col];
}

