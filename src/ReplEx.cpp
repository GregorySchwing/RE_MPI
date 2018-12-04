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


#include "ReplEx.h"
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

using namespace std;

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


static int repl_quantity(struct gmx_repl_ex *re, int ere, float q, ReplicaExchangeParameters* replExParams){
   
    
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
    
    re->nrepl   =   nnodes;
    re->repl    =   nodeid;
    re->nst     =   replExParams->exchangeInterval;
    re->nex     =   replExParams->numExchanges;
    
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
            fprintf(fplog, "\nReplica exchange in temperature\n");
            fflush(fplog);
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5.1f", re->q[re->type][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        default:
            printf("Unknown replica exchange quantity");
            exit(EXIT_FAILURE);
    }

    if (replExParams->randomSeed == -1)
    {   if (re->repl == 0)
        {
            re->seed = rand();
        }
        else
        {
            re->seed = 0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &(re->seed), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    } else {
            re->seed    =   replExParams->randomSeed;
    }
    
    fprintf(fplog, "\nReplica exchange interval: %d\n", re->nst);
    fprintf(fplog, "\nReplica random seed: %d\n", re->seed);

    re->nattempt[0] = 0;
    re->nattempt[1] = 0;

    re->prob_sum = (float*)malloc((re->nrepl)*sizeof(float));
    re->nexchange = (int*)malloc((re->nrepl)*sizeof(int));
    re->nmoves = (int**)malloc((re->nrepl)*sizeof(int*));
    
    for (i = 0; i < re->nrepl; i++)
    {
        re->nmoves[i] = (int*)malloc((re->nrepl)*sizeof(int));
    }
    fprintf(fplog, "Replica exchange information below: ex and x = exchange, pr = probability\n");
    fflush(fplog);
    /* generate space for the helper functions so we don't have to snew each time */


    re->destinations = (int*)malloc((re->nrepl)*sizeof(int));
    re->incycle      = (int*)malloc((re->nrepl)*sizeof(bool));
    re->tmpswap      = (int*)malloc((re->nrepl)*sizeof(int));
    re->cyclic = (int**)malloc((re->nrepl)*sizeof(int*));
    re->order = (int**)malloc((re->nrepl)*sizeof(int*));

    for (int i = 0; i < re->nrepl; i++){
        re->cyclic[i] = (int*)malloc((re->nrepl)*sizeof(int));
        re->order[i] = (int*)malloc((re->nrepl)*sizeof(int));
    } 

    /* allocate space for the functions storing the data for the replicas */
    /* not all of these arrays needed in all cases, but they don't take
       up much space, since the max size is nrepl**2 */

    re->prob = (float*)malloc((re->nrepl)*sizeof(float));
    re->bEx = (int*)malloc((re->nrepl)*sizeof(int));
    re->beta = (float*)malloc((re->nrepl)*sizeof(float));
    re->Vol = (double*)malloc((re->nrepl)*sizeof(double));
    re->Epot = (float*)malloc((re->nrepl)*sizeof(float));
    
    re->de = (float**)malloc((re->nrepl)*sizeof(float*));
    
    for (int i = 0; i < re->nrepl; i++){
        re->de[i] = (float*)malloc((re->nrepl)*sizeof(float));

    }
    
    for (int i = 0; i < re->nrepl; i++){
        re->nexchange[i] = 0;
    }
    
    for (int i = 0; i < re->nrepl; i++){
        for (int j = 0; j < re->nrepl; j++){
            re->nmoves[i][j] = 0;
        }
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

bool replica_exchange(FILE *fplog, struct gmx_repl_ex *re,
                          ReplicaState *stateGlobal, float enerd, double vol_par,
                          int step, ReplicaExchangeParameters* replExParams)
{
    int j;
    int replica_id = 0;
    int exchange_partner;
    int maxswap = 0;
    /* Number of rounds of exchanges needed to deal with any multiple
     * exchanges. */
    /* Where each replica ends up after the exchange attempt(s). */
    /* The order in which multiple exchanges will occur. */
    int bThisReplicaExchanged = false;

//    if (MASTER(cr))
    replica_id  = re->repl;

    // GJS figure out where vol is stored
    double vol = 0.0;
    test_for_replica_exchange(fplog, re, enerd, vol_par, step, replExParams);
    
    prepare_to_do_exchange(re, replica_id, &maxswap, &bThisReplicaExchanged);
    
    if (bThisReplicaExchanged)
    {
            /* There will be only one swap cycle with standard replica
             * exchange, but there may be multiple swap cycles if we
             * allow multiple swaps. */

            for (j = 0; j < maxswap; j++)
            {
                exchange_partner = re->order[replica_id][j];

                if (exchange_partner != replica_id)
                {
                    /* Exchange the global states between the master nodes */
                    //if (debug)
                    //{
                        //fprintf(debug, "Exchanging %d with %d\n", replica_id, exchange_partner);
                    //}
                    //exchange_state(ms, exchange_partner, state);
                }
            }
    }
    
    
    return 1;
}

static void
test_for_replica_exchange(FILE                          *fplog,
                          struct gmx_repl_ex            *re,
                          float                         enerd,
                          double                         vol,
                          int                           step,
                          ReplicaExchangeParameters*    replExParams){ 
    
    int                                  m, i, j, a, b, ap, bp, i0, i1, tmp;
    float                                 delta = 0;
    int                             bPrint, bMultiEx;
    int                            *bEx      = re->bEx;
    float                                *prob     = re->prob;
    int                                 *pind     = re->destinations; /* permuted index */
    int                             bEpot    = false;
    int                             bDLambda = false;
    int                             bVol     = false;
    
    bMultiEx = (re->nex > 1);  /* multiple exchanges at each state */
    fprintf(fplog, "Replica exchange at step %d\n", step);
    
    if ((re->type == ereTEMP || re->type == ereTL))
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->Epot[i] = 0;
        }
        bEpot              = true;
       
        /* temperatures of different states*/
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->q[ereTEMP][i]*BOLTZ);
        }
    }
    else {
    
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->temp*BOLTZ);  /* we have a single temperature */
        }
    }
    
     if (bEpot)
    {
         re->Epot[re->repl] =   enerd;
         MPI_Allreduce(MPI_IN_PLACE, re->Epot, re->nrepl, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    
    if (bEpot)
    /* make a duplicate set of indices for shuffling */
    for (i = 0; i < re->nrepl; i++)
    {
        pind[i] = re->ind[i];
    }
    
    /* make a duplicate set of indices for shuffling */
    for (i = 0; i < re->nrepl; i++)
    {
        re->prob[i] = 0.0;
        re->prob_sum[i] = 0.0;
        re->bEx[i] = 0;
    }
    
    /* standard nearest neighbor replica exchange */

    m = (step / re->nst) % 2;
    
    for (i = 1; i < re->nrepl; i++)
    {
        a = re->ind[i-1];
        b = re->ind[i];

        bPrint = (re->repl == a || re->repl == b);
        
        if (i % 2 == m)
        {
        
        delta = calc_delta(fplog, bPrint, re, a, b, a, b);
            
            if (delta <= 0)
            {
                /* accepted */
                prob[i] = 1;
                bEx[i]  = 1;
            }
            else
            {
                if (delta > PROBABILITYCUTOFF)
                {
                    prob[i] = 0;
                }
                else
                {
                    prob[i] = exp(-delta);
                }
                // roll a number to determine if accepted. For now it is superfluous to
                // reset, but just in case we ever add more calls in different branches
                // it is safer to always reset the distribution.
                //uniformRealDist.reset();
                float randVar = ((float) rand() / (RAND_MAX)) + 1;
                bEx[i] = (int)(randVar < prob[i]);
            }
            
            re->prob_sum[i] += prob[i];

            if (bEx[i])
            {
                /* swap these two */
                tmp       = pind[i-1];
                pind[i-1] = pind[i];
                pind[i]   = tmp;
                re->nexchange[i]++;  /* statistics for back compatibility */
            }
        }
        else
        {
            prob[i] = -1;
            bEx[i]  = 0;
        }
    }
    
        /* print some statistics */
    print_ind(fplog, "ex", re->nrepl, re->ind, bEx);
    print_prob(fplog, "pr", re->nrepl, prob);
    fprintf(fplog, "\n");
    re->nattempt[m]++;

    /* record which moves were made and accepted */
    for (i = 0; i < re->nrepl; i++)
    {
        re->nmoves[re->ind[i]][pind[i]] += 1;
        re->nmoves[pind[i]][re->ind[i]] += 1;
    }
    fflush(fplog); /* make sure we can see what the last exchange was */
    
    return;
}

static void
prepare_to_do_exchange(struct gmx_repl_ex *re,
                       const int           replica_id,
                       int                *maxswap,
                       int           *bThisReplicaExchanged)
{
    //printf("called p4_replica_exchange\n");
    cout.flush();
    int i, j;
    /* Hold the cyclic decomposition of the (multiple) replica
     * exchange. */
    int bAnyReplicaExchanged = false;
    *bThisReplicaExchanged = false;

    for (i = 0; i < re->nrepl; i++)
    {
        if (re->destinations[i] != re->ind[i])
        {
            /* only mark as exchanged if the index has been shuffled */
            bAnyReplicaExchanged = true;
            break;
        }
    }
    if (bAnyReplicaExchanged)
    {
        /* reinitialize the placeholder arrays */
        for (i = 0; i < re->nrepl; i++)
        {
            for (j = 0; j < re->nrepl; j++)
            {
                re->cyclic[i][j] = -1;
                re->order[i][j]  = -1;
            }
        }

        /* Identify the cyclic decomposition of the permutation (very
         * fast if neighbor replica exchange). */
        cyclic_decomposition(re->destinations, re->cyclic, re->incycle, re->nrepl, maxswap);

        /* Now translate the decomposition into a replica exchange
         * order at each step. */
        compute_exchange_order(re->cyclic, re->order, re->nrepl, *maxswap);

        /* Did this replica do any exchange at any point? */
        for (j = 0; j < *maxswap; j++)
        {
            if (replica_id != re->order[replica_id][j])
            {
                *bThisReplicaExchanged = true;
                break;
            }
        }
    }
}

static void print_ind(FILE *fplog, const char *leg, int n, int *ind, int *bEx)
{
    int i;

    fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %c %2d", (bEx != NULL && bEx[i]) ? 'x' : ' ', ind[i]);
    }
    fprintf(fplog, "\n");
}

static void print_prob(FILE *fplog, const char *leg, int n, float *prob)
{
    int  i;
    char buf[8];

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf, "%4.2f", prob[i]);
            fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog, "     ");
        }
    }
    fprintf(fplog, "\n");
}

static float calc_delta(FILE *fplog, int bPrint, struct gmx_repl_ex *re, int a, int b, int ap, int bp)
{
    float   ediff, dpV, delta = 0;
    float  *Epot = re->Epot;
    double  *Vol  = re->Vol;
    float **de   = re->de;
    float  *beta = re->beta;
    
    /* Two cases; we are permuted and not.  In all cases, setting ap = a and bp = b will reduce
       to the non permuted case */

    switch (re->type)
    {
        case ereTEMP:
            /*
             * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
             */

            ediff = Epot[b] - Epot[a];
    
            delta = -(beta[bp] - beta[ap])*ediff;

            break;
        case ereLAMBDA:
            /* two cases:  when we are permuted, and not.  */
            /* non-permuted:
               ediff =  E_new - E_old
                     =  [H_b(x_a) + H_a(x_b)] - [H_b(x_b) + H_a(x_a)]
                     =  [H_b(x_a) - H_a(x_a)] + [H_a(x_b) - H_b(x_b)]
                     =  de[b][a] + de[a][b] */

            /* permuted:
               ediff =  E_new - E_old
                     =  [H_bp(x_a) + H_ap(x_b)] - [H_bp(x_b) + H_ap(x_a)]
                     =  [H_bp(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a) + H_a(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_b(x_b) + H_b(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a)] - [H_ap(x_a) - H_a(x_a)] + [H_ap(x_b) - H_b(x_b)] - H_bp(x_b) - H_b(x_b)]
                     =  (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b])    */
            /* but, in the current code implementation, we flip configurations, not indices . . .
               So let's examine that.
                     =  [H_b(x_ap) - H_a(x_a)] - [H_a(x_ap) - H_a(x_a)] + [H_a(x_bp) - H_b(x_b)] - H_b(x_bp) - H_b(x_b)]
                     =  [H_b(x_ap) - H_a(x_ap)]  + [H_a(x_bp) - H_b(x_pb)]
                     = (de[b][ap] - de[a][ap]) + (de[a][bp] - de[b][bp]
                     So, if we exchange b<=> bp and a<=> ap, we return to the same result.
                     So the simple solution is to flip the
                     position of perturbed and original indices in the tests.
             */

            ediff = (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b]);
            delta = ediff*beta[a]; /* assume all same temperature in this case */
            break;
        case ereTL:
            /* not permuted:  */
            /* delta =  reduced E_new - reduced E_old
                     =  [beta_b H_b(x_a) + beta_a H_a(x_b)] - [beta_b H_b(x_b) + beta_a H_a(x_a)]
                     =  [beta_b H_b(x_a) - beta_a H_a(x_a)] + [beta_a H_a(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + beta_b H_a(x_a) - beta_a H_a(x_a)] +
                        [beta_a dH_a(x_b) + beta_a H_b(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + [beta_a dH_a(x_b) +
                        beta_b (H_a(x_a) - H_b(x_b)]) - beta_a (H_a(x_a) - H_b(x_b))
                     =  beta_b dH_b(x_a) + beta_a dH_a(x_b) - (beta_b - beta_a)(H_b(x_b) - H_a(x_a) */
            /* delta = beta[b]*de[b][a] + beta[a]*de[a][b] - (beta[b] - beta[a])*(Epot[b] - Epot[a]; */
            /* permuted (big breath!) */
            /*   delta =  reduced E_new - reduced E_old
                     =  [beta_bp H_bp(x_a) + beta_ap H_ap(x_b)] - [beta_bp H_bp(x_b) + beta_ap H_ap(x_a)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                        - beta_pb H_a(x_a) + beta_ap H_a(x_a) + beta_pb H_a(x_a) - beta_ap H_a(x_a)
                        - beta_ap H_b(x_b) + beta_bp H_b(x_b) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [(beta_bp H_bp(x_a) - beta_bp H_a(x_a)) - (beta_ap H_ap(x_a) - beta_ap H_a(x_a))] +
                        [(beta_ap H_ap(x_b)  - beta_ap H_b(x_b)) - (beta_bp H_bp(x_b) - beta_bp H_b(x_b))]
             + beta_pb H_a(x_a) - beta_ap H_a(x_a) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [beta_bp (H_bp(x_a) - H_a(x_a)) - beta_ap (H_ap(x_a) - H_a(x_a))] +
                        [beta_ap (H_ap(x_b) - H_b(x_b)) - beta_bp (H_bp(x_b) - H_b(x_b))]
             + beta_pb (H_a(x_a) - H_b(x_b))  - beta_ap (H_a(x_a) - H_b(x_b))
                     =  ([beta_bp de[bp][a] - beta_ap de[ap][a]) + beta_ap de[ap][b]  - beta_bp de[bp][b])
             + (beta_pb-beta_ap)(H_a(x_a) - H_b(x_b))  */
            delta = beta[bp]*(de[bp][a] - de[bp][b]) + beta[ap]*(de[ap][b] - de[ap][a]) - (beta[bp]-beta[ap])*(Epot[b]-Epot[a]);
            break;
        default:
            printf("Unknown replica exchange quantity");
            exit(EXIT_FAILURE);
    }
    
    if (bPrint)
    {

        fprintf(fplog, "Repl %d <-> %d  dE_term = %10.3e (kT)\n", a, b, delta);
    
    }

    if (re->bNPT)
    {
        /* revist the calculation for 5.0.  Might be some improvements. */
        dpV = (beta[ap]*re->pres[ap]-beta[bp]*re->pres[bp])*(Vol[b]-Vol[a])/PRESFAC;
        if (bPrint)
        {
            fprintf(fplog, "  dpV = %10.3e  d = %10.3e\n", dpV, delta + dpV);
        }
        delta += dpV;
    }
    
    return delta;
}

static void
cyclic_decomposition(const int *destinations,
                     int      **cyclic,
                     int  *incycle,
                     const int  nrepl,
                     int       *nswap)
{

    int i, j, c, p;
    int maxlen = 1;
    for (i = 0; i < nrepl; i++)
    {
        incycle[i] = false;
    }
    for (i = 0; i < nrepl; i++)  /* one cycle for each replica */
    {
        if (incycle[i])
        {
            cyclic[i][0] = -1;
            continue;
        }
        cyclic[i][0] = i;
        incycle[i]   = true;
        c            = 1;
        p            = i;
        for (j = 0; j < nrepl; j++) /* potentially all cycles are part, but we will break first */
        {
            p = destinations[p];    /* start permuting */
            if (p == i)
            {
                cyclic[i][c] = -1;
                if (c > maxlen)
                {
                    maxlen = c;
                }
                break; /* we've reached the original element, the cycle is complete, and we marked the end. */
            }
            else
            {
                cyclic[i][c] = p;  /* each permutation gives a new member of the cycle */
                incycle[p]   = true;
                c++;
            }
        }
    }
    *nswap = maxlen - 1;
/*
    if (debug)
    {
            for (i = 0; i < nrepl; i++)
            {
                fprintf(debug, "Cycle %d:", i);
                for (j = 0; j < nrepl; j++)
                {
                    if (cyclic[i][j] < 0)
                    {
                        break;
                    }
                    fprintf(debug, "%2d", cyclic[i][j]);
                }
                fprintf(debug, "\n");
            }
            fflush(debug);
    }
*/
}

static void
compute_exchange_order(int     **cyclic,
                       int     **order,
                       const int nrepl,
                       const int maxswap)
{
    int i, j;

    for (j = 0; j < maxswap; j++)
    {
        for (i = 0; i < nrepl; i++)
        {
            if (cyclic[i][j+1] >= 0)
            {
                order[cyclic[i][j+1]][j] = cyclic[i][j];
                order[cyclic[i][j]][j]   = cyclic[i][j+1];
            }
        }
        for (i = 0; i < nrepl; i++)
        {
            if (order[i][j] < 0)
            {
                order[i][j] = i; /* if it's not exchanging, it should stay this round*/
            }
        }
    }
/*
       if (debug)
        {
            fprintf(debug, "Replica Exchange Order\n");
            for (i = 0; i < nrepl; i++)
            {
                fprintf(debug, "Replica %d:", i);
                for (j = 0; j < maxswap; j++)
                {
                    if (order[i][j] < 0)
                    {
                        break;
                    }
                    fprintf(debug, "%2d", order[i][j]);
                }
                fprintf(debug, "\n");
            }
            fflush(debug);
        }
*/
}

void print_replica_exchange_statistics(FILE *fplog, struct gmx_repl_ex *re)
{
    int  i;

    fprintf(fplog, "\nReplica exchange statistics\n");

    if (re->nex == 0)
    {
        fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
                re->nattempt[0]+re->nattempt[1], re->nattempt[1], re->nattempt[0]);

        fprintf(fplog, "Repl  average probabilities:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, NULL);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "Repl  number of exchanges:\n");
        print_ind(fplog, "", re->nrepl, re->ind, NULL);
        print_count(fplog, "", re->nrepl, re->nexchange);

        fprintf(fplog, "Repl  average number of exchanges:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  ((float)re->nexchange[i])/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, NULL);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "\n");
    }
    /* print the transition matrix */
    print_transition_matrix(fplog, re->nrepl, re->nmoves, re->nattempt);
}

static void print_transition_matrix(FILE *fplog, int n, int **nmoves, int *nattempt)
{
    int   i, j, ntot;
    float Tprint;

    ntot = nattempt[0] + nattempt[1];
    fprintf(fplog, "\n");
    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "    ");  /* put the title closer to the center */
    }
    fprintf(fplog, "Empirical Transition Matrix\n");

    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%8d", (i+1));
    }
    fprintf(fplog, "\n");

    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "Repl");
        for (j = 0; j < n; j++)
        {
            Tprint = 0.0;
            if (nmoves[i][j] > 0)
            {
                Tprint = nmoves[i][j]/(2.0*ntot);
            }
            fprintf(fplog, "%8.4f", Tprint);
        }
        fprintf(fplog, "%3d\n", i);
    }
}

static void print_count(FILE *fplog, const char *leg, int n, int *count)
{
    int i;

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %4d", count[i]);
    }
    fprintf(fplog, "\n");
}

static void exchange_floats(int b, float *v, int n)
{
    float *buf;
    int   i;

    if (v)
    {
        //snew(buf, n);
        /*
           MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(float), MPI_BYTE, MSRANK(ms, b), 0,
                      MPI_COMM_WORLD, &mpi_req);
            MPI_Recv(buf, n*sizeof(float), MPI_BYTE, MSRANK(ms, b), 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        //sfree(buf);
    }
}


static void exchange_doubles(int b, double *v, int n)
{
    double *buf;
    int     i;

    if (v)
    {
        //snew(buf, n);
        /*
           MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                      MPI_COMM_WORLD, &mpi_req);
            MPI_Recv(buf, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        //sfree(buf);
    }
}

/*static void exchange_rvecs(int b, rvec *v, int n)
{
    rvec *buf;
    int   i;

    if (v)
    {
        //snew(buf, n);
        //
        //   MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
        //   buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
        //   ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         //
        {
            MPI_Request mpi_req;

            MPI_Isend(v[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                      MPI_COMM_WORLD, &mpi_req);
            MPI_Recv(buf[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
        for (i = 0; i < n; i++)
        {
            //copy_rvec(buf[i], v[i]);
        }
       // sfree(buf);
    }
}*/

static void exchange_state(int b, ReplicaState *state)
{
    /* When t_state changes, this code should be updated. 
    int ngtc, nnhpres;
    ngtc    = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    exchange_rvecs(ms, b, state->box, DIM);
    exchange_rvecs(ms, b, state->box_rel, DIM);
    exchange_rvecs(ms, b, state->boxv, DIM);
    exchange_reals(ms, b, &(state->veta), 1);
    exchange_reals(ms, b, &(state->vol0), 1);
    exchange_rvecs(ms, b, state->svir_prev, DIM);
    exchange_rvecs(ms, b, state->fvir_prev, DIM);
    exchange_rvecs(ms, b, state->pres_prev, DIM);
    exchange_doubles(ms, b, state->nosehoover_xi.data(), ngtc);
    exchange_doubles(ms, b, state->nosehoover_vxi.data(), ngtc);
    exchange_doubles(ms, b, state->nhpres_xi.data(), nnhpres);
    exchange_doubles(ms, b, state->nhpres_vxi.data(), nnhpres);
    exchange_doubles(ms, b, state->therm_integral.data(), state->ngtc);
    exchange_doubles(ms, b, &state->baros_integral, 1);
    exchange_rvecs(ms, b, state->x.rvec_array(), state->natoms);
    exchange_rvecs(ms, b, state->v.rvec_array(), state->natoms);    */
}