#ifndef _replica_state_h
#define _replica_state_h


class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;
class Ewald;
#include "EnergyTypes.h"


class ReplicaState
{
    public:
        /*//! Constructor
        ReplicaState(){
          potential = NULL;
          coordinates = NULL; //ex
                com = NULL; //ex
          calcEwald = NULL; //ex
          cellList = NULL; //ex
        }*/
        
        ReplicaState(   SystemPotential* potential,
                        Coordinates* coordinates,
                        COM* com,
                        Ewald* *calcEwald,
                        CellList* cellList){
            
            this->potential     =   potential;
            this->coordinates   =   coordinates;
            this->com           =   com;
            this->calcEwald     =   calcEwald;
            this->cellList      =   cellList;
            
        }
          SystemPotential* potential; //ex
          Coordinates* coordinates; //ex
          COM* com; //ex
          Ewald* *calcEwald; //ex
          CellList* cellList; //ex
};

#endif
