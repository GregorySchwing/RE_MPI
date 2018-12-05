/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#ifndef _replica_exchange_helper_h
#define _replica_exchange_helper_h

#include <mpi.h>
#include "ReplicaState.h"

class ReplicaExchangeHelper {

public:    
    
ReplicaExchangeHelper(ReplicaState * re){
    int counter = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            this->interTens[counter]   =   re->potential->boxVirial->interTens[i][j];
            this->realTens[counter]    =   re->potential->boxVirial->realTens[i][j];
            this->recipTens[counter]   =   re->potential->boxVirial->recipTens[i][j];
            this->totalTens[counter]   =   re->potential->boxVirial->totalTens[i][j];
            this->corrTens[counter]    =   re->potential->boxVirial->corrTens[i][j];
            counter++;
        }
    }
    
    Build_derived_type_virial(  &re->potential->boxVirial->inter,
                                &re->potential->boxVirial->tc,
                                &re->potential->boxVirial->real,                                 
                                &re->potential->boxVirial->recip,
                                &re->potential->boxVirial->self,
                                &re->potential->boxVirial->correction,
                                &re->potential->boxVirial->totalElect,
                                &re->potential->boxVirial->total,
                                this->interTens,
                                this->realTens,
                                this->recipTens,
                                this->totalTens,
                                this->corrTens,
                                &this->MPI_VIRIAL
            );
}

void Build_derived_type_virial(
                            double*         inter,  //  in
                            double*         tc,  //  in
                            double*         real,  //  in
                            double*         recip,  //  in
                            double*         self,  //  in
                            double*         correction,  //  in
                            double*         totalElect,  //  in
                            double*         total,  //  in
                            double*        interTens,
                            double*        realTens,
                            double*        recipTens,
                            double*        totalTens,
                            double*        corrTens,        
                            MPI_Datatype*   mesg_mpi_t_ptr  //  out
                            //  pointer to new MPI type
)   {
    /* The number of elements in each "block" of the */
    /* new type. For us, 1 each. */
    int block_lengths[13];
    
    /* Displacement of each element from start of new */
    /* type. The "d_i's." */
    /* MPI_Aint ("address int") is an MPI defined C */
    /* type. Usually an int or a long int. */
    MPI_Aint displacements[13];
    
    /* MPI types of the elements. The "t_i's." */
    MPI_Datatype typelist[13];
    
    /* Use for calculating displacements    */
    MPI_Aint start_address;
    MPI_Aint address;
    
    for (int i = 0; i < 13; i++){
        if (i < 8)
            block_lengths[i]    =   1;
        if (i >= 8
            block_lengths[i]    =   9;
    }
    
            
    /* Build a derived datatype consisting of */
    /* three double arrays*/
            
    /* First element, a, is at displacement 0 */
    displacements[0] = 0;
    
    /* Calculate other displacements relative to a */
    MPI_Address(inter, &start_address);
    
    /* Find address of y_ptr and displacement from x_ptr */
    MPI_Address(tc, &address);
    displacements[1] = address - start_address;
    
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(real, &address);
    displacements[2] = address - start_address;
    
     /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(recip, &address);
    displacements[3] = address - start_address;
    
      
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(self, &address);
    displacements[4] = address - start_address;
    
     /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(correction, &address);
    displacements[5] = address - start_address;
      
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(totalElect, &address);
    displacements[6] = address - start_address;
    
     /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(total, &address);
    displacements[7] = address - start_address;
    
      
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(interTens, &address);
    displacements[8] = address - start_address;
    
     /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(realTens, &address);
    displacements[9] = address - start_address;
      
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(recipTens, &address);
    displacements[10] = address - start_address;
    
     /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(totalTens, &address);
    displacements[11] = address - start_address;
      
    /* Find address of z_ptr and displacement from x_ptr */
    MPI_Address(corrTens, &address);
    displacements[12] = address - start_address;

    for (int i = 0; i < 13; i++){
        typelist[i] = MPI_DOUBLE;
    }

    /* Build the derived datatype */
    MPI_Type_struct(12, block_lengths, displacements,
    typelist, mesg_mpi_t_ptr);
    
    MPI_Type_commit(mesg_mpi_t_ptr);
}

double        interTens[9];
double        realTens[9];
double        recipTens[9];
double        totalTens[9];
double        corrTens[9];
MPI_Datatype    MPI_VIRIAL;

};
#endif