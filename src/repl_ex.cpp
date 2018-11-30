/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

//#include "gmxpre.h"

#include "repl_ex.h"
#include <mpi.h>

#include <cmath>
//#include <random>
#include <string>


#define KILO (1e3)
#define AVOGADRO         (6.02214129e23)
#define BOLTZMANN (1.3806488e-23)  /* (J/K, NIST 2010 CODATA */
#define RGAS        (BOLTZMANN*AVOGADRO)   /* (J/(mol K))  */
#define BOLTZ   (RGAS/KILO)    /* (kJ/(mol K)) */
#define PROBABILITYCUTOFF 100
#define BAR_MDUNITS      (1e5*NANO*PICO*PICO/AMU)
#define PRESFAC          (1.0/BAR_MDUNITS)
#define NANO             (1e-9)                            /* A Number  */
#define PICO             (1e-12)                           /* A Number  */
#define AMU              (1.660538921e-27)                 /* kg, NIST 2010 */

#define PROBABILITYCUTOFF 100
/* we don't bother evaluating if events are more rare than exp(-100) = 3.7x10^-44 */

// Rank in the multisimulation
#define MSRANK(ms, nodeid)  (nodeid)

enum {
    ereTEMP, ereLAMBDA, ereENDSINGLE, ereTL, ereNR
};
const char *erename[ereNR] = { "temperature", "lambda", "end_single_marker", "temperature and lambda"};
/* end_single_marker merely notes the end of single variable replica exchange. All types higher than
   it are multiple replica exchange methods */
/* Eventually, should add 'pressure', 'temperature and pressure', 'lambda_and_pressure', 'temperature_lambda_pressure'?;
   Let's wait until we feel better about the pressure control methods giving exact ensembles.  Right now, we assume constant pressure  */

typedef struct gmx_repl_ex
{
    int       repl;        /* replica ID */
    int       nrepl;       /* total number of replica */
    float      temp;       /* temperature */
    int       type;        /* replica exchange type from ere enum */
    float    **q;          /* quantity, e.g. temperature or lambda; first index is ere, second index is replica ID */
    int     bNPT;          /* use constant pressure and temperature */
    double     *pres;       /* replica pressures */
    int      *ind;         /* replica indices */
    int      *allswaps;    /* used for keeping track of all the replica swaps */
    int       nst;         /* replica exchange interval (number of steps) */
    int       nex;         /* number of exchanges per interval */
    int       seed;        /* random seed */
    int       nattempt[2]; /* number of even and odd replica change attempts */
    float     *prob_sum;   /* sum of probabilities */
    int     **nmoves;      /* number of moves between replicas i and j */
    int     *nexchange;    /* i-th element of the array is the number of exchanges between replica i-1 and i */
    int      natoms;       /* GJS - REDS - 2/14/18 */

    /* these are helper arrays for replica exchange; allocated here so they
       don't have to be allocated each time */
    int      *destinations;
    int     **cyclic;
    int     **order;
    int      *tmpswap;
    int *incycle;
    int *bEx;

    /* helper arrays to hold the quantities that are exchanged */
    float  *prob;
    float  *Epot;
    float  *beta;
    double  *Vol;
    float **de;

} t_gmx_repl_ex;


int repl_quantity(struct gmx_repl_ex *re, int ere, float q, ReplicaExchangeParameters* replExParams){
   
    
    int         bDiff;
    int         s;

    int nnodes;
    int nodeid;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &nodeid);  
  

    float       *qall = (float*)malloc(nnodes*sizeof(float));
    for(int i = 0; i < nnodes; i++){
        qall[i] = 0.0;
    }
   
    qall[nodeid] = q;
    
    MPI_Allreduce(MPI_IN_PLACE, qall, nnodes, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        

//#endif
    
bDiff = false;
    for (s = 1; s < nnodes; s++)
    {
        if (qall[s] != qall[0])
        {
            bDiff = true;
        }
    }

    if (bDiff)
    {
       //// Set the replica exchange type and quantities ///
        re->type = ere;

        for (s = 0; s < nnodes; s++)
        {
            re->q[ere][s] = qall[s];
            printf("replica %d qall[%d] = %f\n", nodeid, s, qall[s]);
        }
      
    }
    free(qall);
    return bDiff;
    
}

gmx_repl_ex_t
init_replica_exchange(FILE                            *fplog,
                      float                             temp,
                      ReplicaExchangeParameters * replExParams)
{
    double                pres;
    int                 i, j, k;
    gmx_repl_ex * re;
    int            bTemp;
    
    int nnodes;
    int nodeid;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &nodeid);
    
    re = (gmx_repl_ex*)malloc(sizeof(gmx_repl_ex));
    
    re->q = (float**)malloc(sizeof(float*)*ereENDSINGLE);
    re->q[ereTEMP] = (float*)malloc(sizeof(float)*nnodes);
    
    bTemp    = repl_quantity(re, ereTEMP, temp, replExParams);
    
    if (re->type == -1)  /* nothing was assigned */
    {
        printf("Error: The properties of the %d systems are all the same, there is nothing to exchange\n", re->nrepl);;
             exit(EXIT_FAILURE);
    }
    
    /* Make an index for increasing replica order */
    /* only makes sense if one or the other is varying, not both!
       if both are varying, we trust the order the person gave. */
    
    re->nrepl = nnodes;
    re->repl = nodeid;
    
    re->ind = (int*)malloc((re->nrepl)*sizeof(int));
    for (int i = 0; i < re->nrepl; i++)
    {
        re->ind[i] = i;
    }
    
    if (re->type < ereENDSINGLE)
    {

        for (i = 0; i < re->nrepl; i++)
        {
            for (j = i+1; j < re->nrepl; j++)
            {
                if (re->q[re->type][re->ind[j]] < re->q[re->type][re->ind[i]])
                {
                    /* Unordered replicas are supposed to work, but there
                     * is still an issues somewhere.
                     * Note that at this point still re->ind[i]=i.
                     */
                    printf("Replicas with indices %d < %d have %ss %g > %g, please order your replicas on increasing %s",
                              i, j,
                              erename[re->type],
                              re->q[re->type][i], re->q[re->type][j],
                              erename[re->type]);
                    exit(EXIT_FAILURE);
                }
                else if (re->q[re->type][re->ind[j]] == re->q[re->type][re->ind[i]])
                {
                    printf("Two replicas have identical %ss", erename[re->type]);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    /* keep track of all the swaps, starting with the initial placement. */
    re->allswaps = (int*)malloc((re->nrepl)*sizeof(int));
    for (i = 0; i < re->nrepl; i++)
    {
        re->allswaps[i] = re->ind[i];
    }

    switch (re->type)
    {
        case ereTEMP:
            printf("\nReplica exchange in temperature\n");
            fprintf(fplog, "\nReplica exchange in temperature\n");
            fflush(fplog);
    /*        for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5.1f", re->q[re->type][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        default:
            printf("Unknown replica exchange quantity");
            exit(EXIT_FAILURE);
            */
    }

    
/* 
    int nnodes;
    int nodeid;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &nodeid);  
    
    int         bDiff;
    int         s;

    float       *qall = (float*)malloc(nnodes*sizeof(float));
    for(int i = 0; i < nnodes; i++){
        qall[i] = 0.0;
    }
   
    qall[nodeid] = temp;
    
    MPI_Allreduce(MPI_IN_PLACE, qall, nnodes, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
*/        
    return re;
 
}