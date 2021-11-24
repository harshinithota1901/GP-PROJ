/*  lil-gp Genetic Programming System, version 1.0, 11 July 1995
 *  Copyright (C) 1995  Michigan State University
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of version 2 of the GNU General Public License as
 *  published by the Free Software Foundation.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *  
 *  Douglas Zongker       (zongker@isl.cps.msu.edu)
 *  Dr. Bill Punch        (punch@isl.cps.msu.edu)
 *
 *  Computer Science Department
 *  A-714 Wells Hall
 *  Michigan State University
 *  East Lansing, Michigan  48824
 *  USA
 *  
 *  Modified from mutate.c by Terry Van Belle, May 25, 2001
 */

/********* czj changes
    debugging
    new member in mutate_data structure
      absDepth_czj
    changed operator_mutate_init function
      to read a new lilgp mutation parameter depth_abs which
        when false (which is default) then no changes
        when true then use the minimal depth ramp set by depth parameter
          as the aboslute minimum for grown subtrees so that smaller
          tress are treated as bad (same as those too big)
          however, MinDepth_czj will try to force generate_random_tree_grow
          to grow at least this deep
    changed operator_mutate function
      to facilitate the above
    NOTE: this operator does not check for size violations
***********/

#include <lilgp.h>

typedef struct
{
     int keep_trying;
     double internal;
     double external;
     double *tree;
     double treetotal;
     char *sname;
     sel_context *sc;
     int method;
     int mindepth, maxdepth;
     int absDepth_czj;
     double mutate_probability;  /* new */
     int guaranteed;             /* new */
} mutate_data;

/* operator_uniform_mutate_init()
 *
 * initialize a breedphase record for uniform mutation.  Modified from
 * operator_mutate_init().
 */

int operator_uniform_mutate_init ( char *options, breedphase *bp )
{
     int errors = 0;
     mutate_data *md;
     int i, j, k, m;
     double r;
     char **argv, **targv;
     int internalset = 0, externalset = 0;
     char *cp;

     md = (mutate_data *)MALLOC ( sizeof ( mutate_data ) );

     /* fill in the breedphase record. */
     bp->operator = OPERATOR_UNIFORM_MUTATE;
     bp->data = (void *)md;
     bp->operator_free = operator_uniform_mutate_free;
     bp->operator_start = operator_uniform_mutate_start;
     bp->operator_end = operator_uniform_mutate_end;
     bp->operator_operate = operator_uniform_mutate;

     /* default values for the mutation-specific data structure. */
     md->keep_trying = 0;
     md->internal = 0.9;
     md->external = 0.1;
     md->tree = (double *)MALLOC ( tree_count * sizeof ( double ) );
     for ( j = 0; j < tree_count; ++j )
          md->tree[j] = 0.0;
     md->treetotal = 0.0;
     md->sname = NULL;
     md->method = GENERATE_HALF_AND_HALF;
     md->mindepth = -1;
     md->maxdepth = -1;
	 md->absDepth_czj = 0;
     md->mutate_probability = 0.01;  /* new */
     md->guaranteed = 0;             /* new */
     
     j = parse_o_rama ( options, &argv );

     for ( i = 0; i < j; ++i )
     {
          if ( strcmp ( "keep_trying", argv[i] ) == 0 )
          {
               md->keep_trying = translate_binary ( argv[++i] );
               if ( md->keep_trying == -1 )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"%s\" is not a valid setting for \"keep_trying\".",
                           argv[i] );
               }
          }
          else if ( strcmp ( "internal", argv[i] ) == 0 )
          {
               internalset = 1;
               md->internal = strtod ( argv[++i], NULL );
               if ( md->internal < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"internal\" must be nonnegative." );
               }
          }
          else if ( strcmp ( "external", argv[i] ) == 0 )
          {
               externalset = 1;
               md->external = strtod ( argv[++i], NULL );
               if ( md->external < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"external\" must be nonnegative." );
               }
          }
          else if ( strcmp ( "select", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               FREE ( md->sname );
               md->sname = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( md->sname, argv[i] );
          }
          else if ( strcmp ( "method", argv[i] ) == 0 )
          {
               ++i;
               if ( strcmp ( argv[i], "half_and_half" ) == 0 )
                    md->method = GENERATE_HALF_AND_HALF;
               else if ( strcmp ( argv[i], "grow" ) == 0 )
                    md->method = GENERATE_GROW;
               else if ( strcmp ( argv[i], "full" ) == 0 )
                    md->method = GENERATE_FULL;
               else
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"%s\" is not a known generation method.",
                           argv[i] );
               }
          }
          else if ( strcmp ( "depth", argv[i] ) == 0 )
          {
               md->mindepth = strtol ( argv[++i], &cp, 10 );
               if ( *cp == 0 )
                    md->maxdepth = md->mindepth;
               else if ( *cp == '-' )
               {
                    md->maxdepth = strtol ( cp+1, &cp, 10 );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "uniform mutation: malformed depth string \"%s\".",
                                argv[i] );
                    }
               }
               else
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: malformed depth string \"%s\".",
                           argv[i] );
               }
          }
          else if ( strcmp ( "tree", argv[i] ) == 0 )
          {
               k = parse_o_rama ( argv[++i], &targv );
               if ( k != tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: wrong number of tree fields: \"%s\".",
                           argv[i] );
               }
               else
               {
                    for ( m = 0; m < k; ++m )
                    {
                         md->tree[m] = strtod ( targv[m], &cp );
                         if ( *cp )
                         {
                              ++errors;
                              error ( E_ERROR, "uniform mutation: \"%s\" is not a number.",
                                     targv[m] );
                         }
                    }
               }
               
               free_o_rama ( k, &targv );
          }
          else if ( strncmp ( "tree", argv[i], 4 ) == 0 )
          {
               k = strtol ( argv[i]+4, &cp, 10 );
               if ( *cp )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: unknown option \"%s\".",
                           argv[i] );
               }
               if ( k < 0 || k >= tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"%s\" is out of range.",
                           argv[i] );
               }
               else
               {
                    md->tree[k] = strtod ( argv[++i], &cp );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "uniform mutation: \"%s\" is not a number.",
                                argv[i] );
                    }
               }
          }
          else if ( strcmp ( "mutate_probability", argv[i] ) == 0 ) {  /* new */
               md->mutate_probability = strtod ( argv[++i], NULL );
               if ( md->mutate_probability < 0.0 ||
                    md->mutate_probability > 1.0 )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"mutate_probability\" must be between 0 and 1." );
               }
          }
          else if ( strcmp ( "guaranteed", argv[i] ) == 0 )
          {
               md->guaranteed = translate_binary ( argv[++i] );
               if ( md->guaranteed == -1 )
               {
                    ++errors;
                    error ( E_ERROR, "uniform mutation: \"%s\" is not a valid setting for \"guaranteed\".",
                           argv[i] );
               }
          }
		  /************* begin czj changes */
          else if ( strcmp ( "depth_abs", argv[i] ) == 0 )
          {
               md->absDepth_czj = translate_binary ( argv[++i] );
               if ( md->absDepth_czj == -1 )
               { ++errors;
                 error(E_ERROR,
         "mutation: \"%s\" is not a valid setting for \"depth_abs\".",argv[i]);
               }
          }
/************* end czj changes */
          else
          {
               ++errors;
               error ( E_ERROR, "uniform mutation: unknown option \"%s\".",
                      argv[i] );
          }
     }
     
     free_o_rama ( j, &argv );
     
     if ( internalset && !externalset )
          md->external = 0.0;
     else if ( !internalset && externalset )
          md->internal = 0.0;
     
     if ( md->sname == NULL )
     {
          ++errors;
          error ( E_ERROR, "uniform mutation: no selection method specified." );
     }

     if ( md->mindepth == -1 && md->maxdepth == -1 )
     {
          md->mindepth = 0;
          md->maxdepth = 4;
     }
     
     if ( md->mindepth < 0 || md->maxdepth < 0 ||
         md->maxdepth < md->mindepth )
     {
          ++errors;
          error ( E_ERROR, "uniform mutation: bad depth range.\n" );
     }
     
     for ( j = 0; j < tree_count; ++j )
          md->treetotal += md->tree[j];
     if ( md->treetotal == 0.0 )
     {
          for ( j = 0; j < tree_count; ++j )
               md->tree[j] = 1.0;
          md->treetotal = tree_count;
     }
          
     r = 0.0;
     for ( j = 0; j < tree_count; ++j )
          r = (md->tree[j] += r);

/*#ifdef DEBUG */
     if ( !errors )
     {
          printf ( "uniform mutation options:\n" );
          printf ( "   probability of mutation: %lf\n", md->mutate_probability ); /* new */
          printf ( "   internal: %lf  external: %lf\n", md->internal, md->external );
          printf ( "   keep_trying: %d\n", md->keep_trying );
          printf ( "   selection: %s\n", md->sname==NULL?"NULL":md->sname );
          printf ( "   method: %d    mindepth: %d   maxdepth: %d\n",
                  md->method, md->mindepth, md->maxdepth );
          printf ( "   tree total: %lf\n", md->treetotal );
          for ( j = 0; j < tree_count; ++j )
               printf ( "   tree %d: %lf\n", j, md->tree[j] );
     }
/*#endif*/
     
     return errors;
}

/* operator_uniform_mutate_free()
 *
 * copied from operator_mutate_free()
 */

void operator_uniform_mutate_free ( void *data )
{
     mutate_data * md;

     md = (mutate_data *)data;
     FREE ( md->sname );
     FREE ( md->tree );
     FREE ( md );
}

/* operator_uniform_mutate_start()
 *
 * copied from operator_mutate_start()
 */

void operator_uniform_mutate_start ( population *oldpop, void *data )
{
     mutate_data * md;
     select_context_func_ptr select_con;

     md = (mutate_data *)data;
     select_con = get_select_context ( md->sname );
     md->sc = select_con ( SELECT_INIT, NULL, oldpop, md->sname );
}

/* operator_uniform_mutate_end()
 *
 * copied from operator_mutate_end()
 */

void operator_uniform_mutate_end ( void *data )
{
     mutate_data * md;

     md = (mutate_data *)data;
     md->sc->context_method ( SELECT_CLEAN, md->sc, NULL, NULL );
}

/*
 *  Utility function that doesn't seem to be in the bag of tricks.
 *  I'm probably missing something here.
 */
void copy_to_space(lnode* tree, int space)
{
     int repcount;

     gensp_reset ( space );
     copy_tree_replace_many ( space, tree, NULL, NULL, 0, &repcount );

     if ( repcount != 0 )
     {
          error ( E_FATAL_ERROR, "botched copy_to_space" );
     }
}

/* operator_uniform_mutate()
 *
 * Works as follows:
 * 
 * int size = tree's size
 * for ( i = 0; i < size; i++ )
 * {
 *     if ( random_double() < md->mutate_probability )
 *     {
 *         do normal mutate   
 *     }
 * }
 * 
 */

void operator_uniform_mutate ( population *oldpop, population *newpop,
                           void *data )
{
     int i, j;
     int ps;
     lnode *replace[2];
     int l, ns;
     int badtree;
     int repcount;
     mutate_data * md;
     int t;
     double r;
     int depth;
     int totalnodes;
     int p;
     int forceany;
     double total;
     int mutate_count;
     int newtree_size;
 
     md = (mutate_data *)data;
     total = md->internal + md->external;

     /* choose a tree to mutate. */
     r = random_double() * md->treetotal;
     for ( t = 0; r >= md->tree[t]; ++t );

     /* select an individual to mutate. */
     p = md->sc->select_method ( md->sc ); 
     ps = tree_nodes ( oldpop->ind[p].tr[t].data );
     /* forceany = (ps==1||total==0.0); */

#ifdef DEBUG_MUTATE
     fprintf ( stdout, "the parent size is %d\n", ps );
     fprintf ( stdout, "    parent %4d: ", p );
     print_tree ( oldpop->ind[p].tr[t].data, stdout );
#endif          

     /* copy the parent tree to the offspring position. */
     duplicate_individual ( (newpop->ind)+newpop->next,
                            (oldpop->ind)+p );
    
     /* pre-calculate the number of standard mutations we'll perform */
     if (md->guaranteed) {
          mutate_count = 1;  /* at least one mutation */
     } else {
          mutate_count = 0;
     }
     for (j = 0; j < ps; j++) {
          if (random_double() < md->mutate_probability)
               mutate_count++;
     }
          
     for (j = 0; j < mutate_count; j++)
     {
          newtree_size = tree_nodes(newpop->ind[newpop->next].tr[t].data);

          forceany = (newtree_size == 1 || total==0.0);

	  if ( forceany )
	  {
	       /* choose any point. */
	       l = random_int ( newtree_size );
	       replace[0] = get_subtree ( newpop->ind[newpop->next].tr[t].data, l );
	  }
	  else if ( total*random_double() < md->internal )
	  {
	       /* choose an internal point. */
	       l = random_int ( tree_nodes_internal ( newpop->ind[newpop->next].tr[t].data ) );
	       replace[0] = get_subtree_internal ( newpop->ind[newpop->next].tr[t].data, l );
	  }
	  else
	  {
	       /* choose an external point. */
	       l = random_int ( tree_nodes_external ( newpop->ind[newpop->next].tr[t].data ) );
	       replace[0] = get_subtree_external ( newpop->ind[newpop->next].tr[t].data, l );
	  }
	  
#ifdef DEBUG_MUTATE
          fprintf ( stdout, "newsize = %d, forceany = %d\n", 
		newtree_size, forceany );
          fprintf ( std, "selected for replacement: " );
          print_tree ( replace[0], stdout );
#endif

          gensp_reset ( 1 );
	  /* pick a value from the depth ramp. */
          depth = md->mindepth + random_int ( md->maxdepth - md->mindepth + 1 );
		  if (md->absDepth_czj)
            MinDepth_czj=depth-md->mindepth;
	  /* grow the tree. */
          switch ( md->method )
          {
             case GENERATE_GROW:
               generate_random_grow_tree ( 1, depth, fset+tree_map[t].fset );
               break;
             case GENERATE_FULL:
               generate_random_full_tree ( 1, depth, fset+tree_map[t].fset );
               break;
             case GENERATE_HALF_AND_HALF:
               if ( random_double() < 0.5 )
                    generate_random_grow_tree ( 1, depth, fset+tree_map[t].fset );
               else
                    generate_random_full_tree ( 1, depth, fset+tree_map[t].fset );
               break;
          }
          MinDepth_czj=0;

#ifdef DEBUG_MUTATE
          fprintf ( stdout, "the new subtree is: " );
          print_tree ( gensp[1].data, stdout );

#endif

	  /* count the nodes in the new tree. */
          ns = newtree_size - tree_nodes ( replace[0] ) + tree_nodes ( gensp[1].data );
          totalnodes = ns;

          /* copy the selected tree, replacing the subtree at the
		  mutation point with the randomly generated tree. */
          replace[1] = gensp[1].data;
          copy_tree_replace_many ( 0, newpop->ind[newpop->next].tr[t].data,
                                       replace, replace+1, 1, &repcount );
          if ( repcount != 1 )
          {
               error ( E_FATAL_ERROR,
                      "botched mutation:  this can't happen." );
          }

	  /* free the tree selected for mutation. */
          free_tree ( newpop->ind[newpop->next].tr+t );
	  /* copy the tree to the new individual. */
          gensp_dup_tree ( 0, newpop->ind[newpop->next].tr+t );

          newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
          newpop->ind[newpop->next].flags = FLAG_NONE;
               
#ifdef DEBUG_MUTATE
          fprintf ( stdout, "    the mutated individual is: " );
          print_individual ( newpop->ind+newpop->next, stdout );
#endif
     }

     ++newpop->next;
#ifdef DEBUG_MUTATE
     printf ( "MUTATION COMPLETE.\n\n\n" );
#endif
}


