/* cgp_czj.c version 1.2 implements CGP preprocessing for lil-gp */
/* using lil-gp 1.02                                             */
/* it also allows weights for mutation set members and adapts    */
/* the weights based on distribution of counted F/T              */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lilgp.h" 
#include "types.h"

/* 
#define DEBUG_SORT
#define DEBUG_WHEELS
*/

/* acgp1.1.1 parameters:                                         */
/* Note that variables have the same name except that first .    */
/*   is replaced with _                                          */
/*   acgp.use_trees_prct (0..1]                                  */
/*     the effective rate for distribution sampling              */
/*     that is this prct of all trees (all pops) will be taken   */
/*   acgp.select_all  [0,1]                                      */
/*     1 - extract acgp.use_trees_prct best (after sort) out of  */
/*         each population then take them all for sampling       */
/*     0 - extract sqrt(acgp.use_trees_prct) best of each pop    */
/*         then resort and take again sqrt(acgp.use_trees_prct)  */
/*         resulting in acgp.use_trees_prct effective rate       */
/*   acgp.extract_quality_prct [0..1]                            */
/*     two trees with fitness diff by no more than               */
/*     (1-acgp.extract_quality_prct) is considered same fitness  */
/*     and thus compared on size                                 */
/*   acgp.gen_start_prct [0..1]                                  */
/*     start extracting at gen=acgp.gen_start_prct*MaxGen        */
/*   acgp.gen_step [1..MaxGen]                                   */
/*     after starting extracting, extract at this gen interval   */
/*   acgp.gen_slope [0,1,2]                                      */
/*     0 - use extracted heuristics to update and the rate of    */
/*         change is constant at                                 */
/*         sqrt(1/numIterations) if acgp.gen_slope_prct==0       */
/*         else at acgp.gen_slope_prct                           */
/*         note that acgp.gen_slope==0 and acgp.gen_slope_prct=1 */
/*         give complete greedy 100% change of the heuristics    */
/*     1 - use extracted heuristics to update old heuristics     */
/*         and the rate of change increases with iter#           */
/*   acgp.gen_slope_prct [0..1]                                  */
/*     see above                                                 */
/*   acgp.0_threshold_prct [0..1]                                */
/*     if a weight drops to weight such that weight/mutSetSize   */
/*     less than acgp.threshold_prct, then drop weight to 0      */
/*   acgp.what [0,1,2,3]                                         */
/*     0 - CGP run, no adjustments, no *.cnt/.wgt file           */
/*     1 - CGP, but also compute distribution *.cnt/*.wgt files  */
/*     2 - ACGP run, extract and adjust heuristics on this run   */
/*     3 - as 2, but after adjusting heuristics, regrow all pops */
/*   acgp.stop_on_term [0,1]                                     */
/*     0 - continue all generation even on solving (term met)    */
/*     1 - stop generation upon solving                          */
/*   acgp.use_expressed [0,1,2]                                  */
/*     0 - collect distribution from all lnodes in a tree        */
/*     1 - skip over subtress that are not expressed             */
/*     2 - as 0, but use weight proportional to the number of    */
/*         times the node is expressed (and thus skip w/0)       */



/* default values for ACGP parameters */
#define ACGP_use_trees_prct_DEFAULT 0.1   
#define ACGP_select_all_DEFAULT 1       
#define ACGP_extract_quality_prct_DEFAULT 0.99
#define ACGP_gen_start_prct_DEFAULT 0.0    /* from the beginning */
#define ACGP_gen_step_DEFAULT 1 
#define ACGP_gen_slope_DEFAULT 0             /* constant changes */     
#define ACGP_gen_slope_prct_DEFAULT 0.1      /* change at at 10% */
#define ACGP_0_threshold_prct_DEFAULT 0.025
#define ACGP_what_DEFAULT 0                  /* straight CGP run */
#define ACGP_stop_on_term_DEFAULT 1       
#define ACGP_use_expressed_DEFAULT 0            /* use all nodes */


/* global dependent variables for acgp */
static double Acgp_extract_quality_diff;
           /* set on every gen, the actual diff for same quality */
static int Acgp_use_expressed;
           /* copy of the acgp.use_expressed parameter           */
int Acgp_adj_regrow; 
          /* set on every generation:
             0 - no adjusting not regrowing run or simply not
	         regrow this generation
             1 - new iteration without regrowing this generation
	     2 - new iteration and regrowing this generation     */
int Acgp_stop_on_term;
                      /* copy of the acgp.stop_on_term parameter */


/* local static functions/macros */

static void isort(individual **data,int array_size);

static void acgp_print_counters_czj(int iterNum, const char *basefile, 
                                    int newIteration);
/* print the counters for root and parent-child distribution czjmmh */

static void acgp_reset_count(void);
/* resets the counters for count_czj() */

#define ACGP_EXPRESSED_VAL(expressed_node_count) \
(Acgp_use_expressed==0) ? (1) : ((Acgp_use_expressed==1) ?\
!!(expressed_node_count) : (expressed_node_count))

static void acgp_count_czj ( lnode *data, int toUse);
/* computes the
   distribution of parent-child function(terminal) contexts for
   the single tree 'data'. The individual instances are counted in
   First reset all counters
   Using recursive acgp_count_recurse() to count functions
   Here count the label of Root MS[NumF][0].counters[indexOfRoot]++
   then recurse 
   Counters are counted according to the expressed_czj field and 
   the parameter ACGP_use_expressed 
   If toUse==0 it means that walking through unexpressed subtree
   but must walk due to tree representation */

static void acgp_count_recurse ( lnode **l, int toUse);
/* the recursive utility
   for countFI() above
   The individual instances are counted in
   the MS structure as follow:
   MS[i][a].counters[j]++ if we encounter parent 'i' with
   child 'j' as argument 'a'. 
   Counters are counted according to the expressed_czj field and
   the parameter ACGP_use_expressed 
   If toUse==0 it means that walking through unexpressed subtree
   but must walk due to tree representation */ 


static void acgp_reset_expressed_recurse( lnode **l);
/* recursive utility for acgp_reset_expressed_czj */

static double acgp_new_weight(double oldW, double prctChnge, double statW,
                       int mutSetSize,double);
/* return new weight
   assume: oldW and statW are normalized to probs in their mut set */

static void acgp_print_wghts(int curGen, FILE *fp, int newIteration);
/* print weights from MS_czj into the open file fp */

static void acgp_normWghtsSetWheels(void);
/* normalize weights for each mut set to probability and set the wheel */

static int cmpQuality(const individual *p, const individual *q);
static int cmpSize(const void *p, const void *q);

/* end of local static functions */

/* locally global stuff */

#define SMALL 0.000000001                        /* delta for double compare */
#define MINWGHT 0.00100           /* min. maintanble wght for a member of MS */

static int NumF;                                    /* num functions in fset */
static int NumT;                         /* num typeII/III functions in fset */

const int RepeatsSrc_czj=5;    /* max repeats in crossover on infeasible src */
const int RepeatsBad_czj=5;                  /* same on bad (size) offspring */
static double WghtsExt;   /* sum of wghts of cross-feasible leaves of a tree */
static double WghtsInt;                       /* same for the internal nodes */

MS_czj_t MS_czj;                      /* the global table with mutation sets */

int Function_czj; 
        /* the node to be modified has this parent; if the node  is the root 
	   then this number is the number of type I functions */

int Argument_czj; 
        /* the node to be modified is this argument of its parent
           uninitialized if the node is the root */

int MinDepth_czj;
        /* if depth=m-n is given in eithe initialization or mutation
           and the corresponding depth_abs=true, this is used to grow only
           trees at least as deep as 'm' (if possible) */

typedef int *specArray_t; /* dynamic for NumF+NumT; array[i]==1 indicates that
                             function indexed i is present in a set */
typedef struct 
{ int numF;                         /* number of F set elements in specArray */
  int numT;
  specArray_t members; 
} specArrayArg_t;
typedef specArrayArg_t *specArrays_t;     /* dynamic array of specArrayArg_t */
typedef struct 
{ int arity;
  specArrays_t Tspecs;    /* dynamic array of ptrs to specArray_t for Tspecs */
  specArrays_t Fspec;                                       /* just one here */
  specArrays_t Fspecs;
} constraint_t;                     /* NOTE: Root will use Tspecs and Fspecs */
typedef constraint_t *constraints_t; /* dynamic array for functions and Root */

static constraints_t Constraints;



static int funNumber(const char* funName, int treeNumber)
/* given funName, return its index in fset[treeNumber] or -1 if not found */
{ int i;
  for (i=0; i<fset[treeNumber].size; i++)
    if (!strcmp(funName,fset[treeNumber].cset[i].string))
      return(i);
  return(-1);
}

static void displayHeader(void)
{ printf("\n\n\t\tWELCOME TO acgp/lilgp 1.1/1.02\n");
  printf("\n\t\tdeveloped by\n");
  printf("\tCezary Z. Janikow\n");
  printf("\tUniversity of Missouri - St. Louis\n");
  printf("\temailto:cjanikow@ola.cs.umsl.edu\n");
  printf("\thttp://www.cs.umsl.edu/~janikow\n");
  printf("\n\n\n\n\tThis is distributed as addition to lil-gp\n");
  printf("\n\tNo explicit/implicit warranty\n");
  printf("\n\tPlease send bug reports to above\n");
  printf("\n\n\n");
}

static void displayConstraints(int Ts, int F, int Fs)
	/* for debugging, arguments state which to display */
{ int fun, arg, entry;
  printf("\n\n\t\tCONSTRAINTS\n");
  for (fun=0; fun<NumF; fun++)
  { printf("Function \"%s\" [#%d]:\n",fset[0].cset[fun].string,fun);
    if (F)
    { printf("\tF_%s [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,
          Constraints[fun].Fspec[0].numF,Constraints[fun].Fspec[0].numT);
      for (entry=0; entry<NumF; entry++)
        if (Constraints[fun].Fspec[0].members[entry])
           printf(" %s",fset[0].cset[entry].string);
      printf(" ||");
      for (; entry<NumF+NumT; entry++)
        if (Constraints[fun].Fspec[0].members[entry])
           printf(" %s",fset[0].cset[entry].string);
      printf("\n");
    }
    if (Fs)
      for (arg=0; arg<Constraints[fun].arity; arg++)
      { printf("\tF_%s_%d [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,arg,
          Constraints[fun].Fspecs[arg].numF,Constraints[fun].Fspecs[arg].numT);
        for (entry=0; entry<NumF; entry++)
          if (Constraints[fun].Fspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf(" ||");
        for (; entry<NumF+NumT; entry++)
          if (Constraints[fun].Fspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf("\n");
      } 
    if (Ts)
      for (arg=0; arg<Constraints[fun].arity; arg++)
      { printf("\tT_%s_%d [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,arg,
          Constraints[fun].Tspecs[arg].numF,Constraints[fun].Tspecs[arg].numT);
        for (entry=0; entry<NumF; entry++)
          if (Constraints[fun].Tspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf(" ||");
        for (; entry<NumF+NumT; entry++)
          if (Constraints[fun].Tspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf("\n");
      } 
  }
  printf("Root:\n",fun);
  if (Fs)
  { printf("\tF_Root [#Fs=%d:#Ts=%d] = ",
          Constraints[NumF].Fspecs[0].numF,Constraints[NumF].Fspecs[0].numT);
    for (entry=0; entry<NumF; entry++)
      if (Constraints[NumF].Fspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf(" ||");
    for (; entry<NumF+NumT; entry++)
      if (Constraints[NumF].Fspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf("\n");
  } 
  if (Ts)
  { printf("\tT_Root [#Fs=%d:#Ts=%d] = ",
          Constraints[NumF].Tspecs[0].numF,Constraints[NumF].Tspecs[0].numT);
    for (entry=0; entry<NumF; entry++)
      if (Constraints[NumF].Tspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf(" ||");
    for (; entry<NumF+NumT; entry++)
      if (Constraints[NumF].Tspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf("\n");
  } 
}

static void displayMS(void)
	/* display the mutation sets from MS_czj */
{ int i,j,k;
  printf("\n\nThe following mutation sets were computed...\n\n");
  for (i=0; i<NumF; i++)                      /* First display mutation sets */
  { printf("Function \"%s\" [#%d]: %d mutation set pairs\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++) 
    { printf("\tArgument %d:\n",j);
      printf("\t\tF [%d members] =",MS_czj[i][j].numF);
      for (k=0; k<MS_czj[i][j].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[i][j].members[k]].string);
      printf("\n\t\tT [%d members] =",MS_czj[i][j].numT);
      for (k=0; k<MS_czj[i][j].numT; k++)
        printf(" %s",
          fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string);
      printf("\n");
    }
    printf("\n");
  }
  printf("Root:\n");
  printf("\t\tF [%d members] =",MS_czj[NumF][0].numF);
  for (k=0; k<MS_czj[NumF][0].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[NumF][0].members[k]].string);
  printf("\n\t\tT [%d members] =",MS_czj[NumF][0].numT);
  for (k=0; k<MS_czj[NumF][0].numT; k++)
    printf(" %s",
         fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string);
  printf("\n\n\n");
}

static void displayWeightsWheels(int weights,int wheels)
        /* display weights/wheels from MS_czj */
        /* if weights/wheels==0 then do not display weights/wheels */
{ int i,j,k;
  if (weights==0 && wheels==0)
    return;
  printf("\n\nThese are %s%s...\n\n",weights ? "weights/" : "",
         wheels ? "wheels" : "");
  for (i=0; i<NumF; i++)              
  { printf("Function \"%s\" [#%d]: %d arguments\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++)
    { printf("\tArgument %d: ",j);
      printf("F [%d members, %s]  and T [%d members, %s]\n",MS_czj[i][j].numF,
             MS_czj[i][j].areFs ? "used":"not used", MS_czj[i][j].numT,
             MS_czj[i][j].areTs ? "used":"not used");
      if (weights)
      { printf("\tWeights:");
        for (k=0; k<NumF+NumT; k++)
        /* if (MS_czj[i][j].weights[k]>=MINWGHT) */
            printf(" %.3f",MS_czj[i][j].weights[k]);
        printf("\n");
      }
      if (wheels)
      { printf("\tWheels:");
        for (k=0; k<MS_czj[i][j].numFT; k++)
          printf(" %.3f",MS_czj[i][j].wheel[k]);
        printf("\n");
      }
    }
    printf("\n");
  }
  printf("Root: ");
  printf("F [%d members, %s] and T [%d members, %s]\n",MS_czj[NumF][0].numF,
         MS_czj[NumF][0].areFs ? "used":"not used", MS_czj[NumF][0].numT,
         MS_czj[NumF][0].areTs ? "used":"not used");
  if (weights)
  { printf("\tWeights:");
    for (k=0; k<NumF+NumT; k++)
      /* if (MS_czj[NumF][0].weights[k]>MINWGHT) */
        printf(" %.3f",MS_czj[NumF][0].weights[k]);
    printf("\n");
  }
  if (wheels)
  { printf("\tWheels:");
    for (k=0; k<MS_czj[NumF][0].numFT; k++)
      printf(" %.3f",MS_czj[NumF][0].wheel[k]);
    printf("\n");
  }
  printf("\n");
}

static double readMinWghtto1(const char *prompt)
/* read (0,1], and if entry < MINWGHT then set it as MINWGHT */
{ double what;
  int res;
  printf("%s: ",prompt);
  res=scanf("%lf",&what);
  while (res<1 || what >1 || what <0)
  { if (res==EOF)
      exit(1);
    fflush(stdin);
    printf("\tInvalid weight: %s: ",prompt);
    res=scanf("%lf",&what);
  }
  if (what<MINWGHT)
    what=MINWGHT;                           /* smaller values become minimal */
  return(what);
}

static void readWeightsSetWheels(void)
/* read weights for mutation set entries and construct wheels */
/* assume weights/wheels are set for all members equalweighted */
/* assume weights for non-members are set to -1 */
/* areFs and areTs members of MS_czj are set to true if ate least one $/
/*   members has weight >MINWGHT */
{ int i,j,k;
  double adjWght, totWght;
  int areFs, areTs;
  printf("\n\nSetting initial weights for mutation set members...\n\n");
  printf("\nInitial weights are all equal. Do you accept [0 to change]: ");
  scanf("%d",&i);
  if (i) 
    return;                               /* leave inital weights and wheels */
  for (i=0; i<NumF; i++)                   
  { printf("\n");
    printf("Function \"%s\" [#%d]: %d mutation set pairs\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++)
    { areFs=0; areTs=0; 
      printf("Argument %d:\n",j);
      printf("\tF [%d members] =",MS_czj[i][j].numF);
      for (k=0; k<MS_czj[i][j].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[i][j].members[k]].string);
      printf("\n\tT [%d members] =",MS_czj[i][j].numT);
      for (k=0; k<MS_czj[i][j].numT; k++)
        printf(" %s",
               fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string);
      printf("\n\n\tReading the weights for type I functions...\n");
      for (k=0; k<MS_czj[i][j].numF; k++)
      { printf("\tFunction \"%s\" [%d]: ",
               fset[0].cset[MS_czj[i][j].members[k]].string,
               MS_czj[i][j].members[k]);
        MS_czj[i][j].weights[MS_czj[i][j].members[k]]=
                                            readMinWghtto1("give weight (0,1]");
      }
      printf("\n\tReading the weights for type II/III terminals...\n");
      for (k=0; k<MS_czj[i][j].numT; k++)
      { printf("\tTerminal \"%s\" [%d]: ",
               fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string,
               MS_czj[i][j].members[MS_czj[i][j].numF+k]);
        MS_czj[i][j].weights[MS_czj[i][j].members[MS_czj[i][j].numF+k]]=
                                            readMinWghtto1("give weight (0,1]");
      }
    /* now all non-memb weights are -1, all memb weights are in [MINWGHT..1] */
                   /* now set mut wheel skipping over weights <MINWGHT+SMALL */
      for (k=0,totWght=0; k<MS_czj[i][j].numFT; k++)
        totWght+=MS_czj[i][j].weights[MS_czj[i][j].members[k]];
      for (k=0; k<MS_czj[i][j].numF; k++) 
      { if (MS_czj[i][j].weights[MS_czj[i][j].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][j].weights[MS_czj[i][j].members[k]]/totWght;
          areFs=1;
        }
        MS_czj[i][j].wheel[k]= (k==0) ? adjWght:MS_czj[i][j].wheel[k-1]+adjWght;
      }
      for (k=MS_czj[i][j].numF; k<MS_czj[i][j].numFT; k++)
      { if (MS_czj[i][j].weights[MS_czj[i][j].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][j].weights[MS_czj[i][j].members[k]]/totWght;
          areTs=1;
        }
        MS_czj[i][j].wheel[k]= (k==0) ? adjWght:MS_czj[i][j].wheel[k-1]+adjWght;
      }
      MS_czj[i][j].areFs=areFs;
      MS_czj[i][j].areTs=areTs;
      if (!areFs && !areTs)
      { fprintf(stderr,
                "\tno member of f=%d arg=%d has any weight >MINWGHT\n",i,j);
        exit(1);
      }
    }
    printf("\n\n");
  }
  printf("Root:\n");
  areFs=0; areTs=0;
  printf("\t\tF [%d members] =",MS_czj[NumF][0].numF);
  for (k=0; k<MS_czj[NumF][0].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[NumF][0].members[k]].string);
  printf("\n\t\tT [%d members] =",MS_czj[NumF][0].numT);
  for (k=0; k<MS_czj[NumF][0].numT; k++)
    printf(" %s",
         fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string);
  printf("\n\tReading the weights for type I functions...\n");
  for (k=0; k<MS_czj[NumF][0].numF; k++)
  { printf("\tFunction \"%s\" [%d]: ",
           fset[0].cset[MS_czj[NumF][0].members[k]].string,
           MS_czj[NumF][0].members[k]);
    MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]=
                                            readMinWghtto1("give weight (0,1]");
  }
  printf("\n\tReading the weights for type II/III terminals...\n");
  for (k=0; k<MS_czj[NumF][0].numT; k++)
  { printf("\tTerminal \"%s\" [%d]: ",
           fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string,
           MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]);
    MS_czj[NumF][0].weights[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]]=
                                        readMinWghtto1("give weight (0,1]");
  }
  for (k=0,totWght=0; k<MS_czj[NumF][0].numFT; k++)
      totWght+=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]];
  for (k=0; k<MS_czj[NumF][0].numF; k++)
  { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
      adjWght=0;
    else
    { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]/totWght;
      areFs=1;
    }
    MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
  }
  for (k=MS_czj[NumF][0].numF; k<MS_czj[NumF][0].numFT; k++)
  { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
      adjWght=0;
    else
    { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]/totWght;
      areTs=1;
    }
    MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
  }
  MS_czj[NumF][0].areFs=areFs;
  MS_czj[NumF][0].areTs=areTs;
  if (!areFs && !areTs)
  { fprintf(stderr,"\tno member of Root sets has any weight >MINWGHT\n");
    exit(1);
  }
  printf("\n\n");
}

static void displayFT(void)
{ int i;
  for (i=0; i<25; i++)
    printf("\n");
  printf("%d ordinary functions: \n",NumF);
  for (i=0; i<NumF; i++)
    printf("%5d = %s\n",i,fset[0].cset[i].string); 
  printf("%d terminals (ordinary and ephemeral): \n",NumT);
  for (; i<NumF+NumT; i++)
    printf("%5d = %s\n",i,fset[0].cset[i].string); 
  printf("Separate entries by [ ,;]  Hit <ENTER> for empty set\n");
  printf("Use either function names or numbers, in any order\n\n");
}

static void readOneConstraintSet(const char*prompt, specArrays_t setP, int max)
	/* read one set from one line; max is the highest index allowed */
{ int entry, status;
  char buf[BUFSIZ];
  char sep[]=" ,;\n";
  char *p;
  for (entry=0; entry<NumF+NumT; entry++)
    setP->members[entry]=0;                            /* reset set to empty */
  setP->numF=setP->numT=0;                          /* reset member counters */
  printf("%s [0..%d] = ",prompt,max);
  if (fgets(buf,BUFSIZ,stdin)==NULL)
  { fprintf(stderr,"ERROR: failed reading constrained\n");
    exit(1);
  }
  p=strtok(buf,sep);
  while (p!=NULL)
  { if ((entry=funNumber(p,0))>=0 || (sscanf(p,"%d",&entry)>0))
    { if (entry<0 || entry >max)
        printf("\a\t\t%d entry out of range\n",entry);
      else                                /* entry is a valid function index */
      { setP->members[entry]=1;
        if (entry<NumF) 
          setP->numF++;
        else 
          setP->numT++;
      }
    }
    else               /* failed reading an integer or invalid function name */
      printf("\t\ainvalid entry\n");
    p=strtok(NULL,sep);
  }
}

static void readFTspecs(void)
{ int i, j;
  char prompt[BUFSIZ];
  
  Constraints=(constraints_t)calloc((size_t)(NumF+1),sizeof(constraint_t));
                                                  /* last entry for the Root */
  for (i=0; i<NumF; i++)                   /* first work on type I functions */
  { displayFT();
    printf("Function %d=%s:\n",i,fset[0].cset[i].string);
    Constraints[i].arity=fset[0].cset[i].arity;
    Constraints[i].Fspec=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
    Constraints[i].Fspec[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
    sprintf(prompt,"\tF_%s (exclusions)",fset[0].cset[i].string);
    readOneConstraintSet(prompt,&(Constraints[i].Fspec[0]),NumF-1);
                                            /* Note: type I only for callers */
    Constraints[i].Tspecs=(specArrays_t)calloc((size_t)Constraints[i].arity,
                                               sizeof(specArrayArg_t));
    Constraints[i].Fspecs=(specArrays_t)calloc((size_t)Constraints[i].arity,
                                               sizeof(specArrayArg_t));
    for (j=0; j<Constraints[i].arity; j++)
    { Constraints[i].Tspecs[j].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
      Constraints[i].Fspecs[j].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
      sprintf(prompt,"\tF_%s_%d (exclusions)",fset[0].cset[i].string,j);
      readOneConstraintSet(prompt,&(Constraints[i].Fspecs[j]),NumF+NumT-1);
      sprintf(prompt,"\tT_%s_%d (inclusions)",fset[0].cset[i].string,j);
      readOneConstraintSet(prompt,&(Constraints[i].Tspecs[j]),NumF+NumT-1);
    }
  }
  Constraints[NumF].arity=1;
  Constraints[NumF].Fspec=NULL;
  Constraints[i].Tspecs=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
  Constraints[i].Fspecs=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
  Constraints[i].Tspecs[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
  Constraints[i].Fspecs[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
  displayFT();
  printf("Root:");
  readOneConstraintSet("\tF^Root (exclusions)",&(Constraints[NumF].Fspecs[0]),
                       NumF+NumT-1);
  readOneConstraintSet("\tT^Root (inclusions)",&(Constraints[NumF].Tspecs[0]),
                       NumF+NumT-1);
}

static void generateNF(void)
	/* from specs of Constraints generate NF in Fspecs of Constraints */
	/* this involves creating T-extensive Fspecs */
{ int fun, arg, entry;
  for (fun=0; fun<NumF; fun++)           /* first create T-extensive F-specs */
    for (arg=0; arg<Constraints[fun].arity; arg++)
      for (entry=0; entry<NumF+NumT; entry++)
        if (Constraints[fun].Tspecs[arg].members[entry]==0)
          Constraints[fun].Fspecs[arg].members[entry]=1;
  for (entry=0; entry<NumF+NumT; entry++)               /* same for the Root */
    if (Constraints[NumF].Tspecs[0].members[entry]==0)
      Constraints[NumF].Fspecs[0].members[entry]=1;

  for (fun=0; fun<NumF; fun++)              /* now create F-extensive Fspecs */
    for (entry=0; entry<NumF; entry++)
      if (Constraints[fun].Fspec[0].members[entry]!=0)     /* must extend it */
        for (arg=0; arg<Constraints[entry].arity; arg++)
          Constraints[entry].Fspecs[arg].members[fun]=1;

  for (fun=0; fun<NumF+1; fun++)            /* recount set entries in Fspecs */
    for (arg=0; arg<Constraints[fun].arity; arg++)
    { Constraints[fun].Fspecs[arg].numF=0;
      Constraints[fun].Fspecs[arg].numT=0;
      for (entry=0; entry<NumF; entry++)
        if (Constraints[fun].Fspecs[arg].members[entry]!=0)
          Constraints[fun].Fspecs[arg].numF++;
      for (; entry<NumF+NumT; entry++)
        if (Constraints[fun].Fspecs[arg].members[entry]!=0)
          Constraints[fun].Fspecs[arg].numT++;
    }
}

static void generateMS(void)
	/* create MS from the Fspecs part (ie F-intensive) of Constraints */
{ int fun, arg, entry, k;
  MS_czj=(MS_czj_t)calloc((size_t)(NumF+1),sizeof(mutationSets_czj_t));        
                                            /* one (last) entry for the Root */

  for (fun=0; fun<NumF; fun++)                   /* set all type I functions */
  { MS_czj[fun]=(mutationSets_czj_t)calloc((size_t)fset[0].cset[fun].arity,
                                           sizeof(mutationSet_czj_t));
    for (arg=0; arg<fset[0].cset[fun].arity; arg++)
    { MS_czj[fun][arg].numF=NumF-Constraints[fun].Fspecs[arg].numF;
      MS_czj[fun][arg].numT=NumT-Constraints[fun].Fspecs[arg].numT;
      MS_czj[fun][arg].numFT=MS_czj[fun][arg].numF+MS_czj[fun][arg].numT;
      if (MS_czj[fun][arg].numFT==0)
        fprintf(stderr,"\n\tBoth sets empty for function %d=%s argument %d\n\a",
                fun,fset[0].cset[fun].string,arg);
      MS_czj[fun][arg].members=
         (int*)calloc((size_t)(MS_czj[fun][arg].numFT),sizeof(int));
      MS_czj[fun][arg].weights=
         (double*)calloc((size_t)(NumF+NumT),sizeof(double));
      MS_czj[fun][arg].counters=                         /* czjmmh */
         (int*)calloc((size_t)(NumF+NumT),sizeof(int));
      MS_czj[fun][arg].wheel=
         (double*)calloc((size_t)(MS_czj[fun][arg].numFT),sizeof(double));
      for (entry=0,k=0; k<NumF+NumT; k++)
        if (Constraints[fun].Fspecs[arg].members[k]==0)
        { MS_czj[fun][arg].members[entry]=k;
          MS_czj[fun][arg].weights[k]=1.0/MS_czj[fun][arg].numFT;
          MS_czj[fun][arg].wheel[entry]= (entry==0) ? 
            MS_czj[fun][arg].weights[k] : 
            MS_czj[fun][arg].wheel[entry-1]+MS_czj[fun][arg].weights[k];
          entry++;
        }
        else
          MS_czj[fun][arg].weights[k]= -1.0;
      MS_czj[fun][arg].areFs= !!MS_czj[fun][arg].numF;
      MS_czj[fun][arg].areTs= !!MS_czj[fun][arg].numT;
    }
  }
  MS_czj[NumF]=(mutationSets_czj_t)calloc((size_t)1,
                    sizeof(mutationSet_czj_t));              /* for the Root */
  MS_czj[NumF][0].numF=NumF-Constraints[NumF].Fspecs[0].numF;
  MS_czj[NumF][0].numT=NumT-Constraints[NumF].Fspecs[0].numT;
  MS_czj[NumF][0].numFT=MS_czj[NumF][0].numF+MS_czj[NumF][0].numT;
  if (MS_czj[NumF][0].numFT==0)
  { printf("\n\tBoth Root sets empty - no valid programs exist\n\a");
    exit(1);
  }
  MS_czj[NumF][0].members=
      (int*)calloc((size_t)(MS_czj[NumF][0].numFT),sizeof(int));
  MS_czj[NumF][0].weights=
      (double*)calloc((size_t)(NumF+NumT),sizeof(double));
  MS_czj[NumF][0].counters=                         /* czjmmh */
      (int*)calloc((size_t)(NumF+NumT),sizeof(int));
  MS_czj[NumF][0].wheel=
      (double*)calloc((size_t)(MS_czj[NumF][0].numFT),sizeof(double));
  for (entry=0,k=0; k<NumF+NumT; k++)
    if (Constraints[NumF].Fspecs[0].members[k]==0)
    { MS_czj[NumF][0].members[entry]=k;
      MS_czj[NumF][0].weights[k]=1.0/MS_czj[NumF][0].numFT;
      MS_czj[NumF][0].wheel[entry]=(entry==0) ?
            MS_czj[NumF][0].weights[k] :
            MS_czj[NumF][0].wheel[entry-1]+MS_czj[NumF][0].weights[k];
      entry++;
    }
    else
      MS_czj[NumF][0].weights[k]= -1.0;
  MS_czj[NumF][0].areFs= !!MS_czj[NumF][0].numF;
  MS_czj[NumF][0].areTs= !!MS_czj[NumF][0].numT;
  acgp_reset_count();
}

void print_fset_czj(const char *message)
{ int i;
  printf("%s\n",message);
  for (i=0; i<fset[0].size; i++)
  printf("Position %2d: fset[0].cset[%2d].%-6.6s index=%-2d arity=%-2d \n",i,i,
          fset[0].cset[i].string,fset[0].cset[i].index,fset[0].cset[i].arity);
}

void create_MS_czj(void)
        /* will access global fset function table, and will allocate and
           initialize MS_czj; check that no pair is completely empty */
	/* Must be called after fset is created and ordered,
           but before initializing the population */
{ int what=0;
  NumF=fset[0].function_count;                                      /* |F_I| */
  NumT=fset[0].terminal_count;                           /* |F_II| + |F_III| */
  displayHeader();
  readFTspecs();
/*  printf("\nRead the following constraints...\n");
    displayConstraints(1,1,1); 
*/
  generateNF(); 
  printf("\nThe normal constraints are...\n");
  displayConstraints(0,0,1);
  generateMS();
  displayMS();
  readWeightsSetWheels();
  printf("\nWheels are...\n"); 
  displayWeightsWheels(1,1); 
  printf("Send 1 to continue, anything else to quit cgp-lil-gp: ");
  scanf("%d",&what);
  if (what!=1)
    exit(1);
  printf("\n\n");
}

static int spinWheel(int startI, int stopI, double *wheel)
/* spin the wheel returning an index between startI and stopI inclusively,
   with probability proportional to wheel allocation (roulette) */
{ double mark,begining;
  begining=startI ? wheel[startI-1] : 0;
  mark=begining+random_double()*(wheel[stopI]-begining);
  while (mark > wheel[startI]) startI++;
  return(startI);
}

int random_F_czj(void)
/* return a random type I index from fset, but which appear in the */
/*    mutation set for Function_czj/Argument_czj */
/* if the set is empty, call random_FT_czj() */
/* NOTE: set is considered empty if numF==0 or each weight is <MINWGHT+SMALL */
{ int randIndex;
  randIndex=MS_czj[Function_czj][Argument_czj].numF;
  if (randIndex==0 || MS_czj[Function_czj][Argument_czj].areFs==0)
    return(random_FT_czj());    
  randIndex=spinWheel(0,randIndex-1,MS_czj[Function_czj][Argument_czj].wheel);
  return MS_czj[Function_czj][Argument_czj].members[randIndex];
}
 
int random_T_czj(void)
/* as random_F_czj, except that extract members of T */
{ int randIndex;
  if (MS_czj[Function_czj][Argument_czj].numT==0 || 
      MS_czj[Function_czj][Argument_czj].areTs==0)
    return(random_FT_czj());           
  randIndex=spinWheel(MS_czj[Function_czj][Argument_czj].numF,
                      MS_czj[Function_czj][Argument_czj].numFT-1,
                      MS_czj[Function_czj][Argument_czj].wheel); 
  return MS_czj[Function_czj][Argument_czj].members[randIndex]; 
}

int random_FT_czj(void)
        /* return a random type I/II/III index from fset, but which appear in
           the mutation set for Function_czj/Argument_czj */
        /* assume both sets (F and T) are not empty */
{ int randIndex;
  if (MS_czj[Function_czj][Argument_czj].numFT==0)
  { fprintf(stderr,"\nERROR: both sets should not be empty\n");
    exit(1);
  }
  if (MS_czj[Function_czj][Argument_czj].areFs==0 &&
      MS_czj[Function_czj][Argument_czj].areTs==0)
  { fprintf(stderr,"\nERROR: both sets should not have empty weights\n");
    exit(1);
  }
  randIndex=spinWheel(0,MS_czj[Function_czj][Argument_czj].numFT-1,
                      MS_czj[Function_czj][Argument_czj].wheel);
  return MS_czj[Function_czj][Argument_czj].members[randIndex];
} 

static int markXNodes_recurse_czj( lnode **t )
/* assume Function_czj and Argument_czj are set to indicate dest point */
/* mark and count all feasible source nodes for crossover in tree */
/* marking is done with the corresponding weights w/r to dest parent */
/*   and wheel values are accumulated */
/* clear all other marking weights and wheel values to 0 */
/* sum the weights of feasible internal/external nodes in WghtsInt/WghtsWxt */
/* return the number of marked nodes */
/* NOTE: wheel entry may be the same as that of the last node if this node */
/*   is infeasible -> this will ensure that this node is not selected later */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  double *wheelInt_czj=&((**t).wheelInt_czj);
  int j;
  double wght=0;
  int total;

  ++*t;                              /* step the pointer over the function. */

  if ( f->arity == 0 )                                    /* it is external */
  { if (f->ephem_gen)
               ++*t;                                  /* etra value to skip */
    wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
    if (wght<(MINWGHT+SMALL))     /* not in this mutation set or do not use */
      total=0;
    else
    { WghtsExt+=wght;
      total=1;
    }
    *wheelInt_czj=WghtsInt;
    *wheelExt_czj=WghtsExt;
    return total;
  }
  switch (f->type)                           /* here only for internal nodes */
  { case FUNC_DATA: case EVAL_DATA:
      wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
      if (wght<(MINWGHT+SMALL))    /* not in this mutation set or do not use */
        total=0;
      else
      { WghtsInt+=wght;
        total=1;
      }
      *wheelInt_czj=WghtsInt;
      *wheelExt_czj=WghtsExt;
      for (j=0; j<f->arity; ++j)
        total+=markXNodes_recurse_czj(t);    /* t has already been advanced */
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
      if (wght<(MINWGHT+SMALL))    /* not in this mutation set or do not use */
        total=0;
      else
      { WghtsInt+=wght;
        total=1;
      }
      *wheelInt_czj=WghtsInt;
      *wheelExt_czj=WghtsExt;
      for (j=0; j<f->arity; ++j)
      { ++*t;                   /* skip the pointer over the skipsize node. */
        total+=markXNodes_recurse_czj(t);
      }
      break;
  } /* switch */
  return total;
}

int markXNodes_czj( lnode *data )
/* assume Function_czj and Argument_czj are set to indicate dest point */
/* mark all nodes in tree which are feasible sources with their wghts */
/*   while contructing the wheels for internal and external nodes */
/* accummulate total wght marked in WghtsInt and WghtsExt */
/*   for the internal and the external nodes, respectively */
/* return the number of marked nodes */
{ lnode *t=data;
  WghtsInt=0;
  WghtsExt=0;
  return markXNodes_recurse_czj (&t);
}

static lnode *getSubtreeMarked_recurse_czj(lnode **t, double mark)
/* assume feasible internal and external nodes are marked with wheel values */
/*   and 'mark' determines which wheel entry is used */
/* this function spins the wheel looking for any node */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  double *wheelInt_czj=&((**t).wheelInt_czj);
  lnode *r;
  int i;

  if (mark < (*wheelExt_czj + *wheelInt_czj))      
    return *t;                 /* this might be either internal or external */
  ++*t;                                              /* move t to next node */
  if (f->arity==0)
  { if (f->ephem_gen)
      ++*t;                                /* skip over the terminal nodes. */
    return NULL;                                          /* not found here */
  }
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarked_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                      /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarked_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

static lnode *getSubtreeMarkedInt_recurse_czj(lnode **t, double mark)
/* same as getSubtreeMarked_recurse_czj except look only internal nodes */
{ function *f = (**t).f;
  double *wheelInt_czj=&((**t).wheelInt_czj);
  lnode *r;
  int i;

  if (f->arity==0)                               /* it is external, skip it */
  { ++*t; 
    if (f->ephem_gen)
      ++*t; 
    return NULL;
  }
  if (mark < *wheelInt_czj)                           /* return this node */
      return *t;
  ++*t;                                              /* move t to next node */
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarkedInt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                      /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarkedInt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

static lnode *getSubtreeMarkedExt_recurse_czj(lnode **t, double mark)
/* same as getSubtreeMarked_recurse_czj except look only external nodes */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  lnode *r;
  int i;

  if (f->arity==0)                           /* it is external, check it out */
  { if (mark<*wheelExt_czj)                            /* return this node */
        return *t;
    ++*t; 
    if (f->ephem_gen)
      ++*t;
    return NULL;
  }
  ++*t;                        /* if here than it is an internal node - skip */
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarkedExt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                         /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarkedExt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

lnode *getSubtreeMarked_czj(lnode *data, int intExt)
/* assume tree is filled with both internal and external wheels */
/*   accummulated in WghtsInt and WghtsExt and at least one node is feasible */
/* return a node with selection prob. proportional to its weight */
/* if no nodes found return NULL */
/* if intExt==0 then looking for both internal and external */
/* if intExt==1 then looking for internal, and switch to external only if */
/*   no internal marked */
/* the opposite for intExt==2 */
{ lnode *el = data;
  if (intExt==0 || intExt==1 && WghtsInt<SMALL ||
      intExt==2 && WghtsExt<SMALL)
    return              /* override 'intExt' parameter and look for any node */
      getSubtreeMarked_recurse_czj(&el,(WghtsInt+WghtsExt)*random_double());
  if (intExt==1)
    return
      getSubtreeMarkedInt_recurse_czj(&el,WghtsInt*random_double());
  return 
    getSubtreeMarkedExt_recurse_czj(&el,WghtsExt*random_double());
}

static int verify_tree_czj_recurse ( lnode **t )
/* return #times the tree pointed by t violates MS_czj */
/*   note: *t always points at a function node here: save the function */
{ function *f = (**t).f;
  int i;
  int total=0;

  if (MS_czj[Function_czj][Argument_czj].weights[f->index]<0)
    total++;                                                     /* invalid */
  ++*t;                              /* step the pointer over the function. */
  if (f->arity==0)
  { if (f->ephem_gen)
      ++*t;    /* skip the pointer over the ERC value if this node has one. */
    return total; 
  }
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; ++i)
      { Function_czj = f->index;
        Argument_czj = i;
        total+=verify_tree_czj_recurse (t);
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; ++i)
      { ++*t;                 /* skip the pointer over the skipsize node. */
        Function_czj = f->index;
        Argument_czj = i;
        total+=verify_tree_czj_recurse (t);
      }
      break;
  } /* switch */
  return total;
}

int verify_tree_czj ( lnode *tree )
/* return #times the tree pointed by tree violates MS_czj */
{ lnode *t = tree;
  int f_save=Function_czj;
  int a_save=Argument_czj;
  int total;

  Function_czj = fset[0].function_count;
  Argument_czj = 0;                                  /* start from the Root */
  total=verify_tree_czj_recurse (&t);
  Function_czj=f_save;
  Argument_czj=a_save;
  return total;
}

static void acgp_print_counters_czj(int curGen, const char *basefile, 
                                    int newIteration)
/* print the counters for root and parent-child distribution czj */
/* NOTE: relying on implicit file close */
{ int i, a, j;
  char fname[BUFSIZ];
  static FILE *fpczj=0;

  strcpy(fname,basefile);
  strcat(fname,".cnt");

  if (!fpczj)
  { if ((fpczj=fopen(fname,"w"))==0)
    { fprintf(stderr,"Couldnt open Counter.txt\n");
      exit(1);
    }
    fprintf(fpczj,"F %d T %d ",NumF,NumT);
    for (i=0; i<NumF+NumT; i++)
      fprintf(fpczj,"%s %d ",fset[0].cset[i].string,fset[0].cset[i].arity);
    fprintf(fpczj,"\n");
  }
  fprintf(fpczj,"Gen %3d : %6d",curGen,MS_czj[NumF][0].counters_tot);
  for (i=0; i<NumF+NumT; i++)
    fprintf(fpczj,"%6.5s",fset[0].cset[i].string);
  fprintf(fpczj,"\n");
  for (i=0; i<NumF; i++)        /* function counters */
    for (a=0; a<fset[0].cset[i].arity; a++)
    { fprintf(fpczj,"%5.5s %1d : %6d ",fset[0].cset[i].string,a,
              MS_czj[i][a].counters_tot); 
      for (j=0; j<NumF+NumT; j++)
        fprintf(fpczj,"%5d ",MS_czj[i][a].counters[j]);
      fprintf(fpczj,"\n");
    }
  fprintf(fpczj,"%5s   : %6d ","Root",MS_czj[NumF][0].counters_tot);
  for (i=0; i<NumF+NumT; i++)             /* Root counters */
    fprintf(fpczj,"%5d ",MS_czj[NumF][0].counters[i]);
  fprintf(fpczj,"\n");
  if (newIteration)
    fprintf(fpczj,"*New iteration\n");
}

static void acgp_reset_count(void)
/* Added 5/30/03 by czjmmh. 
   This function resets the counters for count_czj() */
{ int i, a, j;
  for (i=0; i<NumF+NumT; i++)   /* reset Root counters */
    MS_czj[NumF][0].counters[i]=0;
  MS_czj[NumF][0].counters_tot=0;
  for (i=0; i<NumF; i++)        /* reset function counters */
    for (a=0; a<fset[0].cset[i].arity; a++)
    { for (j=0; j<NumF+NumT; j++)
        MS_czj[i][a].counters[j]=0;
      MS_czj[i][a].counters_tot=0;
    }
}

static void acgp_count_czj ( lnode *data, int toUse)
/* This function computes the
   distribution of parent-child function(terminal) contexts for
   the single tree 'data'. The individual instances are counted in
   First reset all counters
   Using recursive acgp_count_recurse() to count functions 
   Here count the label of Root MS[NumF][0].counters[indexOfRoot]++ 
   then recurse 
   Counters are counted according to the expressed_czj field and 
   the parameter ACGP_use_expressed 
   If toUse==0 it means that walking through unexpressed subtree
   but must walk due to tree representation */ 
{
     lnode *c = data;
     int i, a, j;
     int parIndex=(*c).f->index;
     int expressed_count=ACGP_EXPRESSED_VAL((*c).expressed_czj);

     toUse = toUse && expressed_count;
     if (toUse)
     { MS_czj[NumF][0].counters[parIndex]+=expressed_count;
       MS_czj[NumF][0].counters_tot+=expressed_count;
     }
     acgp_count_recurse (&c, toUse);
}

static void acgp_count_recurse ( lnode **l, int toUse)
/* This is the recursive utility
   for countFI() above 
   The individual instances are counted in
   the MS structure as follow:
   MS[i][a].counters[j]++ if we encounter parent 'i' with
   child 'j' as argument 'a'
   Counters are counted according to the expressed_czj field and
   the parameter ACGP_use_expressed 
   If toUse==0 it means that walking through unexpressed subtree
   but must walk due to tree representation */ 
   /* NOTE: redo for proper # of trees */
   /* NOTE: a child cannot be expressed more than its parent  
      because each evaluation starts with resetting the expressed
      field */
{    function *parent;
     int i;
     int parIndex, chIndex;
     int expressed_count=ACGP_EXPRESSED_VAL((**l).expressed_czj);

     toUse = toUse && expressed_count;
     parent = (**l).f;
     parIndex=parent->index;
     ++*l;                      /* l is now first child */
     /* added to fix ERC issue - jja **************************************/
      if (parent->arity==0)
      { if (parent->ephem_gen)
          ++*l;  /* skip the pointer over the ERC value if this node has one. */
      } /* end of ERC fix code - jja **************************************/
     switch ( parent->type )
     {  case FUNC_DATA:
        case EVAL_DATA:
          for ( i = 0; i < parent->arity; ++i )
          {    chIndex=(**l).f->index;
               expressed_count=ACGP_EXPRESSED_VAL((**l).expressed_czj);
               if (toUse && expressed_count)
                                          /* expressed subtree && child */
               { MS_czj[parIndex][i].counters[chIndex]+=expressed_count;
                 MS_czj[parIndex][i].counters_tot+=expressed_count;
               }
               acgp_count_recurse(l,toUse);  
          }
          break;
        case FUNC_EXPR:
        case EVAL_EXPR:
          for ( i = 0; i < parent->arity; ++i )
          {    ++*l;           /* l was skipnode, now its the child */
               chIndex=(**l).f->index;
               expressed_count=ACGP_EXPRESSED_VAL((**l).expressed_czj);
               if (toUse && expressed_count)
               { MS_czj[parIndex][i].counters[chIndex]+=expressed_count;
                 MS_czj[parIndex][i].counters_tot+=expressed_count;
               }
               acgp_count_recurse(l,toUse);
          }
          break;
     }
}

void acgp_reset_expressed_czj(lnode *data)
/* resets the field expressed_czj in every node
   even in unexpressed subtrees (possible doe to crossover) 
   NOTE: when resetting a node can be expressed even if parent wasn't
   due to crossover */
{
     lnode *c = data;
     acgp_reset_expressed_recurse (&c);
}

static void acgp_reset_expressed_recurse(lnode **l)
/* recursive utility for above */
{
     function *parent;
     int i;
     int parIndex, chIndex;
     (**l).expressed_czj=0;

     parent = (**l).f;
     parIndex=parent->index;
     ++*l;                      /* l is now first child */
     /* added to fix ERC issue - jja ****************************************/
      if (parent->arity==0)
      { if (parent->ephem_gen)
          ++*l;  /* skip the pointer over the ERC value if this node has one. */
      } /* end of ERC fix code - jja ****************************************/
     switch ( parent->type )
     {
        case FUNC_DATA:
        case EVAL_DATA:
          for ( i = 0; i < parent->arity; ++i )
             acgp_reset_expressed_recurse(l);
          break;
        case FUNC_EXPR:
             
        case EVAL_EXPR:
          for ( i = 0; i < parent->arity; ++i )
          {    ++*l;           /* l was skipnode, now its the child */
               acgp_reset_expressed_recurse(l);
          }
          break;
     }
}

static double acgp_new_weight(double oldW, double prctChnge, double statW,
                       int mutSetSize, double thresholdPrct)
/* return new weight
   assume: oldW and statW are normalized to probs in their mut set */
{ double newW;
  newW=oldW*(1-prctChnge)+statW*prctChnge;
  if (newW<thresholdPrct/mutSetSize || newW<MINWGHT)
    return MINWGHT;
  else
    return newW;
}

void acgp_adjustWeights_czj(popstats *popStat, population **pops, int numPops,
       int curGen, int maxGen, const char *basefile)
/* called after every gp generation
   if acgp.what==0 then do nothing, plain CGP run
   first count statistics and print to *.cnt file for distribution information
   - compute counters statistics for all trees in all pops and print
   - extract (int)ceil(popSize*acgp_internal_use_trees_prct) from each pop
     onto the common extracted list
   - sort 
   - if acgp.select_all==1 select all extracted
   - else select (int)ceil(popSize*acgp_internal_use_trees_prct)
   - count statistics and print to *.cnt file
   adjust weights if appropriate
   - adjusting is done every acgp_gen_step generations after 
     acgp_gen_start=ceil(maxGen*acgp.gen_start_prct)
     if acgp.what > 1 
   on adjusting generations
   - use the counters to adjust weights, then normalize then 
     recompute the wheels and print weights into *.wgt
     (initial weights always printed up front)
   numPops is the number of populations
   popStat[0] is accumulative best from all populations */
{ int i,a,j;
  static FILE *fp;
      /* local vars for remaining acgp parameters with their default values */
  static double acgp_use_trees_prct=ACGP_use_trees_prct_DEFAULT;   
  static int acgp_select_all=ACGP_select_all_DEFAULT;
  static double acgp_gen_start_prct=ACGP_gen_start_prct_DEFAULT;
  static int acgp_gen_step=ACGP_gen_step_DEFAULT; 
  static int acgp_gen_slope=ACGP_gen_slope_DEFAULT; 
  static double acgp_gen_slope_prct=ACGP_gen_slope_prct_DEFAULT;
  static double acgp_0_threshold_prct=ACGP_0_threshold_prct_DEFAULT; 
        /* falling below this threshold/mutSetSize will be reset to MINWGHT */
  static double acgp_extract_quality_prct=ACGP_extract_quality_prct_DEFAULT;
               /* 0..1, 1-this number times the diff between best and worst
          fitness are considered same fitness for comparison (sort on size) */
  static int acgp_what=ACGP_what_DEFAULT;          /* def is no adjustments */
  
                                             /* now computed acgp parameter */
  static int acgp_gen_start;                /* generation to start adapting */
  static int popSize;
  static individual **acgp_pop_list;        /* this for sorting a whole pop */
  static int acgp_extract_from_pop;
                                           /* size to extract from each pop */
  static int acgp_extract_list_len;
                                           /* acgp_extract_from_pop*numPops */
  static individual **acgp_extract_list;
                   /* constructd, list of pointers to extracted individuals */
  static double acgp_internal_use_trees_prct;
                 /* acgp_use_trees_prct if numPops==1 || acgp_select_all==1 */
                 /* sqrt(acgp_use_trees_prct) otherwise                     */
  static double acgp_prctChnge;             /* use this prct new of weights */
  double lowFit, highFit;
 
  Acgp_adj_regrow=0;
  if (curGen==0)
  { char fname[BUFSIZ];
    char *p;     
           /* first get acgp params and set dependent params and structures */
    Acgp_stop_on_term=ACGP_stop_on_term_DEFAULT;
    Acgp_use_expressed=ACGP_use_expressed_DEFAULT;
    if (p=get_parameter( "acgp.what" ))
        acgp_what = atoi(p);
    if (acgp_what<0 || acgp_what>3)
    { error ( E_WARNING, "acgp.what reset to default" );
      acgp_what=ACGP_what_DEFAULT;
    }
    if (!acgp_what)
      return;

    popSize=atoi(get_parameter("pop_size"));
    if (p=get_parameter( "acgp.use_trees_prct" ))
      acgp_use_trees_prct = atof(p);
    if (acgp_use_trees_prct<=0 || acgp_use_trees_prct>1)
      error ( E_FATAL_ERROR, "acgp_use_trees_prct must be (0..1>" );
    binary_parameter( "acgp.select_all",  ACGP_select_all_DEFAULT );
    if (p=get_parameter( "acgp.select_all" ))
      acgp_select_all = atoi(p);
    if (numPops==1 || acgp_select_all)
      acgp_internal_use_trees_prct=acgp_use_trees_prct;
    else
      acgp_internal_use_trees_prct=sqrt(acgp_use_trees_prct);
    acgp_extract_from_pop=(int)ceil(popSize*acgp_internal_use_trees_prct);
    acgp_extract_list_len=acgp_extract_from_pop*numPops;
    acgp_pop_list=
           (individual**)MALLOC (popSize*sizeof(individual*));
    acgp_extract_list=
           (individual**)MALLOC (acgp_extract_list_len*sizeof(individual*));
    if (p=get_parameter( "acgp.extract_quality_prct" ))
        acgp_extract_quality_prct = atof(p);
    if (acgp_extract_quality_prct<0 || acgp_extract_quality_prct>1)
    { error ( E_WARNING, "acgp.extract_quality_prct reset to default" );
      acgp_extract_quality_prct=ACGP_extract_quality_prct_DEFAULT;
    }
    if (acgp_what>1)
    { if (p=get_parameter( "acgp.gen_start_prct" ))
        acgp_gen_start_prct = atof(p);
      if (acgp_gen_start_prct<0 || acgp_gen_start_prct>1)
        { error ( E_WARNING, "acgp.gen_base_prct reset to default" );
        acgp_gen_start_prct=ACGP_gen_start_prct_DEFAULT;
      }
      acgp_gen_start=acgp_gen_start_prct*1000000;
      acgp_gen_start=maxGen*acgp_gen_start/1000000;
      if (p=get_parameter( "acgp.gen_step" ))
        acgp_gen_step = atoi(p);
      if (acgp_gen_step<1)
      { error ( E_WARNING, "acgp.gen_step reset to default" );
        acgp_gen_step=ACGP_gen_step_DEFAULT;
      }
      binary_parameter( "acgp.gen_slope",  ACGP_gen_slope_DEFAULT );
      if (p=get_parameter( "acgp.gen_slope" ))
        acgp_gen_slope = atoi(p);
      if (p=get_parameter( "acgp.gen_slope_prct" ))
        acgp_gen_slope_prct = atof(p);
      if (acgp_gen_slope_prct<0 || acgp_gen_slope_prct>1)
      { error ( E_WARNING, "acgp.gen_slope_prct reset to default" );
        acgp_gen_slope_prct=ACGP_gen_slope_prct_DEFAULT;
      }
      if (p=get_parameter( "acgp.0_threshold_prct" ))
        acgp_0_threshold_prct = atof(p);
      if (acgp_0_threshold_prct<0 || acgp_0_threshold_prct>1)
      { error ( E_WARNING, "acgp.0_threshold_prct reset to default" );
        acgp_0_threshold_prct=ACGP_0_threshold_prct_DEFAULT;
      }
      binary_parameter( "acgp.stop_on_term",  ACGP_stop_on_term_DEFAULT );
      if (p=get_parameter( "acgp.stop_on_term" ))
        Acgp_stop_on_term = atoi(p);
      if (p=get_parameter( "acgp.use_expressed" ))
        Acgp_use_expressed = atoi(p);
      if (Acgp_use_expressed<0 || Acgp_use_expressed>2)
      { error ( E_WARNING, "acgp.use_expressed reset to default" );
        Acgp_use_expressed=ACGP_use_expressed_DEFAULT;
      }
      if (acgp_gen_slope==0)             /* constant heuristics change rate */
        if (acgp_gen_slope_prct<SMALL)
          acgp_prctChnge=sqrt((double)acgp_gen_step/(maxGen-acgp_gen_start+1));
        else
          acgp_prctChnge=acgp_gen_slope_prct;
    }
               /* open weight file and print parameters and initial weights */
    strcpy(fname,basefile);
    strcat(fname,".wgt");
    if ((fp=fopen(fname,"w"))==0)
    { fprintf(stderr,"Could not open %s to write weights\n",fname);
      exit(1);
    }
    fprintf(fp,"Parameters:\n");
    fprintf(fp,"\trandom_seed=%d\n",
		!(get_parameter("random_seed")) ? 1 : (atoi(get_parameter("random_seed"))));
    fprintf(fp,"\tmax_generations=%d\n",maxGen);
    fprintf(fp,"\tpop_size=%d\n",popSize);
    fprintf(fp,"\tmultiple_subpops=%d\n",numPops);
    fprintf(fp,"\tacgp.use_trees_prct=%.5f\n",acgp_use_trees_prct);
    fprintf(fp,"\tacgp_internal_use_trees_prct=%.5f\n",
            acgp_internal_use_trees_prct);
    fprintf(fp,"\tacgp.extract_quality_prct=%.5f\n",acgp_extract_quality_prct);
    fprintf(fp,"\tacgp.select_all=%d\n",acgp_select_all);
    fprintf(fp,"\tacgp.what=%d\n",acgp_what);
    fprintf(fp,"\tacgp.stop_on_term=%d\n",Acgp_stop_on_term);
    fprintf(fp,"\tacgp.use_expressed=%d\n",Acgp_use_expressed);
    if (acgp_what>1)
    { fprintf(fp,"\tacgp.gen_start_prct=%.3f\n",acgp_gen_start_prct);
      fprintf(fp,"\tacgp_gen_start=%d\n",acgp_gen_start);
      fprintf(fp,"\tacgp.gen_step=%d\n",acgp_gen_step);
      fprintf(fp,"\tacgp.gen_slope=%d\n",acgp_gen_slope);
      if (acgp_gen_slope==0)
      { fprintf(fp,"\tacgp.gen_slope_prct=%.5f\n",acgp_gen_slope_prct);
        fprintf(fp,"\tconstant change rate at =%.5f\n",acgp_prctChnge);
      }
      fprintf(fp,"\tacgp.0_threshold_prct=%.3f\n",acgp_0_threshold_prct);
    }
    fprintf(fp,"\n");
    fprintf(fp,"%-4s%5.5s%4s%6s","It#","Func","Arg","TotCt");
    for (i=0; i<NumF+NumT; i++)
      fprintf(fp,"%6.5s",fset[0].cset[i].string);
    fprintf(fp,"\n");
    acgp_print_wghts(-1,fp,Acgp_adj_regrow>0); 
                                            /* dumping the original weights */
  }                                   /* end of extra work on genereation 0 */
                                         /* here now on on every generation */
  if (!acgp_what)
    return;
  if (acgp_what>1 && curGen>=acgp_gen_start)
  { if (!((curGen-acgp_gen_start)%acgp_gen_step))
    { Acgp_adj_regrow=1;
      if (acgp_gen_slope==1)
        acgp_prctChnge= (double)(curGen-acgp_gen_start+1)/
                                (maxGen-acgp_gen_start+1);
      if (acgp_what==3)
        Acgp_adj_regrow=2;
    }
  }

/* remove to avoid printing the whole population if want to ****************/
  acgp_reset_count();
  for (i=0; i<numPops; i++)    
    for (j=0; j<popSize; j++)   
      acgp_count_czj(pops[i]->ind[j].tr->data,1);
  acgp_print_counters_czj(curGen,basefile,0);
/***************************************** end remove **********************/

  lowFit=1;      
  highFit=0;      
  for (i=0; i<numPops; i++)                    
  { for (j=0; j<popSize; j++)     /* first extract from each pop after sort */
      acgp_pop_list[j] = pops[i]->ind + j;
    if (popStat[i+1].worstfit<lowFit)
      lowFit=popStat[i+1].worstfit;
    if (popStat[i+1].bestfit>highFit)
      highFit=popStat[i+1].bestfit;
#ifdef DEBUG_SORT
    printf("Pop %d list before sort. Bestfit=%.5f worstfit=%.5f\n",i,
           popStat[i+1].bestfit,popStat[i+1].worstfit);
    for (a=0; a<popSize; a++)
      printf("%3d%8.5f%4d\n",a,acgp_pop_list[a]->a_fitness,
             acgp_pop_list[a]->tr->nodes);
#endif
    Acgp_extract_quality_diff=(popStat[i+1].bestfit-popStat[i+1].worstfit)*
                              (1-acgp_extract_quality_prct);
    isort(acgp_pop_list,popSize);
#ifdef DEBUG_SORT
    printf("Pop %d list after sort based on q_diff=%.5f\n",i,
           Acgp_extract_quality_diff);
    for (a=0; a<popSize; a++)
      printf("%3d%8.5f%4d\n",a,acgp_pop_list[a]->a_fitness,
             acgp_pop_list[a]->tr->nodes);
#endif
    for (j=0; j<acgp_extract_from_pop; j++)  /* and put best on extract list */
      acgp_extract_list[i*acgp_extract_from_pop+j]=acgp_pop_list[j];
  }

  acgp_reset_count();                                 /* now count and print */
#ifdef DEBUG_SORT
  printf("Extract list as constructed. Bestfit=%.5f=%.5f worstfit=%.5f=%.5f\n",
         highFit,popStat[0].bestfit,lowFit,popStat[0].worstfit);
  for (i=0; i<acgp_extract_list_len; i++)
    printf("%3d%8.5f%4d\n",i,acgp_extract_list[i]->a_fitness,
           acgp_extract_list[i]->tr->nodes);
#endif

                 /* if numPops==1 || acgp_select_all then take all extracted */
  if (numPops==1 || acgp_select_all)   
    for (i=0; i<acgp_extract_list_len; i++)
      acgp_count_czj(acgp_extract_list[i]->tr->data,1);
  else                           /* resort and select acgp_internal_use_prct */
  { Acgp_extract_quality_diff=(highFit-lowFit)*(1-acgp_extract_quality_prct);
    isort(acgp_extract_list,acgp_extract_list_len);
#ifdef DEBUG_SORT
    printf("Extract list after additional sort\n");
    for (i=0; i<acgp_extract_list_len; i++)
      printf("%3d%8.5f%4d\n",i,acgp_extract_list[i]->a_fitness,
             acgp_extract_list[i]->tr->nodes);
#endif
    for (i=0; i<acgp_internal_use_trees_prct*acgp_extract_list_len; i++)
      acgp_count_czj(acgp_extract_list[i]->tr->data,1);
  }                                                 

  acgp_print_counters_czj(curGen,basefile, Acgp_adj_regrow>0);  
                                                      /* on every generation */

  if (Acgp_adj_regrow<1)                   /* no adjustments this generation */
  { acgp_print_wghts(curGen,fp,Acgp_adj_regrow>0);
    return;                                   /* just print weights and done */
  }
                           /* Note: we get here only on adjusting generation */
  for (i=0; i<NumF; i++)                      /* work on functions in MS_czj */
  { for (a=0; a<fset[0].cset[i].arity; a++)
      if (MS_czj[i][a].counters_tot)
         for (j=0; j<NumF+NumT; j++)
           if (MS_czj[i][a].weights[j]>0)                  /* if weight used */
             MS_czj[i][a].weights[j]=
                 acgp_new_weight(MS_czj[i][a].weights[j],acgp_prctChnge, 
                 (double)MS_czj[i][a].counters[j]/MS_czj[i][a].counters_tot, 
                         MS_czj[i][a].numFT,acgp_0_threshold_prct);
  }
  if (MS_czj[NumF][0].counters_tot)                   /* change Root weights */
    for (j=0; j<NumF+NumT; j++ )
      if (MS_czj[NumF][0].weights[j]>0)
        MS_czj[NumF][0].weights[j]=
          acgp_new_weight(MS_czj[NumF][0].weights[j],acgp_prctChnge, 
           (double)MS_czj[NumF][0].counters[j]/MS_czj[NumF][0].counters_tot,
                   MS_czj[NumF][0].numFT,acgp_0_threshold_prct);

  acgp_normWghtsSetWheels(); 
#ifdef DEBUG_WHEELS
    printf("Generation %d\n",curGen);
    displayWeightsWheels(1,1); 
#endif
  acgp_print_wghts(curGen,fp,Acgp_adj_regrow>0);
}

static void acgp_print_wghts(int curGen, FILE *fp, int newIteration)
/* print weights from MS_czj into the open file fp */
{ int i, j, a, k;
  for (i=0; i<NumF; i++)
    for (a=0; a<fset[0].cset[i].arity; a++)
    { if (curGen<0) /* initial weights */
        fprintf(fp,"%4s","");
      else
        fprintf(fp,"%-4d",curGen);
      fprintf(fp,"%5.5s%4d %6d",fset[0].cset[i].string,a,
              MS_czj[i][a].counters_tot);
      for (k=0; k<NumF+NumT; k++)
        fprintf(fp,"%6.3f",MS_czj[i][a].weights[k]);
      fprintf(fp, "\n");
    }
  if (curGen<0) /* initial weights */
    fprintf(fp,"%4s","");
  else
    fprintf(fp,"%-4d",curGen);
  fprintf(fp,"%5.5s%4d %6d","Root",0,MS_czj[NumF][0].counters_tot);
  for (k=0; k<NumF+NumT; k++)
    fprintf(fp,"%6.3f",MS_czj[NumF][0].weights[k]);
  fprintf(fp,"\n");
  if (newIteration)
    fprintf(fp,"*New iteration\n");
}

static void acgp_normWghtsSetWheels(void)
/* normalize weights for each mut setthen set the wheel */
/* skip over those with counters_tot=0 (no change) */
/* set areTs/areFs same as in readWeightsSetWheels() */
{ int i,a,k;
  double adjWght, totWght;
  int areFs, areTs;
  for (i=0; i<NumF; i++)                   
  { for (a=0; a<fset[0].cset[i].arity; a++)
    { areFs=0; areTs=0; 
      totWght=0;
      if (!MS_czj[i][a].counters_tot)
        break;
      for (k=0; k<MS_czj[i][a].numFT; k++)
        totWght+=MS_czj[i][a].weights[MS_czj[i][a].members[k]];
                   /* now set mut wheel skipping over weights <MINWGHT+SMALL */
      for (k=0; k<MS_czj[i][a].numF; k++) 
      { if (MS_czj[i][a].weights[MS_czj[i][a].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][a].weights[MS_czj[i][a].members[k]]/totWght;
          areFs=1;
        }
        MS_czj[i][a].wheel[k]= (k==0) ? adjWght:MS_czj[i][a].wheel[k-1]+adjWght;
      }
      for (k=MS_czj[i][a].numF; k<MS_czj[i][a].numFT; k++)
      { if (MS_czj[i][a].weights[MS_czj[i][a].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][a].weights[MS_czj[i][a].members[k]]/totWght;
          areTs=1;
        }
        MS_czj[i][a].wheel[k]= (k==0) ? adjWght:MS_czj[i][a].wheel[k-1]+adjWght;
      }
      MS_czj[i][a].areFs=areFs;
      MS_czj[i][a].areTs=areTs;
      if (!areFs && !areTs)
      { fprintf(stderr,
                "\tno member of f=%d arg=%d has any weight >MINWGHT\n",i,a);
        exit(1);
      }
    }
  }
  areFs=0; areTs=0;
  totWght=0;
  if (MS_czj[NumF][0].counters_tot)
  { for (k=0; k<MS_czj[NumF][0].numFT; k++)
      totWght+=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]];
    for (k=0; k<MS_czj[NumF][0].numF; k++)
    { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
        adjWght=0;
      else
      { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]/totWght;
        areFs=1;
      }
      MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
    }
    for (k=MS_czj[NumF][0].numF; k<MS_czj[NumF][0].numFT; k++)
    { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
        adjWght=0;
      else
      { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]/totWght;
        areTs=1;
      }
      MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
    }
    MS_czj[NumF][0].areFs=areFs;
    MS_czj[NumF][0].areTs=areTs;
    if (!areFs && !areTs)
    { fprintf(stderr,"\tno member of Root sets has any weight >MINWGHT\n");
      exit(1);
    }
  }
}

static int cmpQuality(const individual *p, const individual *q)
/* Compare quality for decreasing order, but if difference less than 
   Acgp_extract_quality_diff then compare sizes
   Could extend the last 'else' to check if size diff is more
   significant than quality diff or not */
{ 
  if ((p->a_fitness - q->a_fitness)>Acgp_extract_quality_diff)
    return -1;
  else if ((q->a_fitness - p->a_fitness)>Acgp_extract_quality_diff)
    return 1;
  else          /* considered same quality, take smaller size as secondary key */
    return (p->tr->nodes - q->tr->nodes);
}

void isort(individual ** data, int array_size)
/* insertion sort redone for a list of individual*  */
{ int i, j;
  individual *indexp;
  for (i=1; i < array_size; i++)
  { indexp = data[i];
    j = i;
    while ((j > 0) && (cmpQuality(data[j-1],indexp)>0))
    { data[j] = data[j-1];
      j--;
    }
    data[j] = indexp;
  }
}

