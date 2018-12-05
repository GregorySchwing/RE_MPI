/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include "ReplDirSetup.h"
#include "ReplEx.cpp"



Simulation::Simulation(char const*const configFileName)
{ 
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  Setup set;
  set.Init(configFileName, &replExParams);
  totalSteps = set.config.sys.step.total;
  staticValues = new StaticVals(set);
  system = new System(*staticValues);
  staticValues->Init(set, *system);
  system->Init(set);
  //recal Init for static value for initializing ewald since ewald is
  //initialized in system
  staticValues->InitOver(set, *system);
  cpu = new CPUSide(*system, *staticValues);
  cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
            totalSteps);

  //Dump combined PSF
  PSFOutput psfOut(staticValues->mol, *system, set.mol.kindMap,
                   set.pdb.atoms.resKindNames);
  psfOut.PrintPSF(set.config.out.state.files.psf.name);
  std::cout << "Printed combined psf to file "
            << set.config.out.state.files.psf.name << '\n';

}

Simulation::~Simulation()
{
  delete cpu;
  delete system;
  delete staticValues;
}



void Simulation::RunSimulation(void)
{
  double startEnergy = system->potential.totalEnergy.total;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeid);
  
  const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
  if(nnodes>1 && useReplicaExchange){
  //if(nnodes>1){
   ReplDirSetup rd(staticValues->forcefield.T_in_K, replExParams);
  
    if (useReplicaExchange) {
      fplog = fopen(rd.path_to_replica_log_file.c_str(), "w");
      replEx = init_replica_exchange(fplog, staticValues->forcefield.T_in_K,
                                     &replExParams);
      
      stateLocal = new ReplicaState(&system->potential, 
                                    &system->coordinates,
                                    &system->com,
                                    &system->calcEwald,
                                    &system->cellList);
    }
  }
  
  for (ulong step = 0; step < totalSteps; step++) {
    system->moveSettings.AdjustMoves(step);
    system->ChooseAndRunMove(step);
    cpu->Output(step);

    if((step + 1) == cpu->equilSteps) {
      double currEnergy = system->potential.totalEnergy.total;
      if(abs(currEnergy - startEnergy) > 1.0e+10) {
        printf("Info: Performing total energy calculation to preserve the"
               " energy information.\n\n");
        system->calcEwald->Init();
        system->potential = system->calcEnergy.SystemTotal();
      }
    }
    
    if(step+1 == totalSteps)
        bLastStep = true;
    
    bDoReplEx = (useReplicaExchange && (step > cpu->equilSteps) && !bLastStep && 
        (step % replExParams.exchangeInterval == 0) && step >= cpu->equilSteps);
  
#if ENSEMBLE == NVT    
    if (bDoReplEx) {
        bExchanged = replica_exchange(fplog, replEx,
                                      stateGlobal, system->potential.totalEnergy.total,
                                      staticValues->boxDimensions->GetTotVolume(),
                                      step, &replExParams);
    }    
#endif
    
#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintTime();
  if (useReplicaExchange)
     print_replica_exchange_statistics(fplog, replEx);
     
}

#ifndef NDEBUG
void Simulation::RunningCheck(const uint step)
{
  system->calcEwald->UpdateVectorsAndRecipTerms();
  SystemPotential pot = system->calcEnergy.SystemTotal();

  std::cout
      << "================================================================="
      << std::endl << "-------------------------" << std::endl
      << " STEP: " << step + 1
      << std::endl << "-------------------------" << std::endl
      << "Energy       INTRA B |     INTRA NB |         INTER |           TC |         REAL |         SELF |   CORRECTION |        RECIP"
      << std::endl
      << "System: "
      << std::setw(12) << system->potential.totalEnergy.intraBond << " | "
      << std::setw(12) << system->potential.totalEnergy.intraNonbond << " | "
      << std::setw(12) << system->potential.totalEnergy.inter << " | "
      << std::setw(12) << system->potential.totalEnergy.tc << " | "
      << std::setw(12) << system->potential.totalEnergy.real << " | "
      << std::setw(12) << system->potential.totalEnergy.self << " | "
      << std::setw(12) << system->potential.totalEnergy.correction << " | "
      << std::setw(12) << system->potential.totalEnergy.recip << std::endl
      << "Recalc: "
      << std::setw(12) << pot.totalEnergy.intraBond << " | "
      << std::setw(12) << pot.totalEnergy.intraNonbond << " | "
      << std::setw(12) << pot.totalEnergy.inter << " | "
      << std::setw(12) << pot.totalEnergy.tc << " | "
      << std::setw(12) << pot.totalEnergy.real << " | "
      << std::setw(12) << pot.totalEnergy.self << " | "
      << std::setw(12) << pot.totalEnergy.correction << " | "
      << std::setw(12) << pot.totalEnergy.recip << std::endl
      << "================================================================"
      << std::endl << std::endl;

}
#endif
/*
void Simulation::initializeCR(t_commrec *cr){
  MPI_Comm_rank(MPI_COMM_WORLD, &cr->nodeid);
  /* Find out number of processes */
 /*  MPI_Comm_size(MPI_COMM_WORLD, &cr->nnodes);
  cr->ms->nsim = cr->nnodes;
  cr->ms->sim = cr->nodeid;
  cr->ms->mpi_comm_masters = MPI_COMM_WORLD;
}*/
