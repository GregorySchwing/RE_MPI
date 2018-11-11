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
//#include "gromacs/mdtypes/commrec.h"

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
  
  /*
  gmx_repl_ex_t     repl_ex = NULL;
  t_commrec *cr = (t_commrec*)malloc(sizeof(t_commrec));
  initializeCR(cr);
  */
  
  FILE *fplog;
  
  int nnodes;
  MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
  
  if(nnodes>1){
      
    directory_stream << "temp_" << staticValues->forcefield.T_in_K;
    directory_name = directory_stream.str();
    replica_directory = opendir(directory_name.c_str());
    
    if(replica_directory){
    
        //printf("Directory already exists : %s\n", directory_name.c_str());
        /* Do whatever gromacs does here w the backups #1, #2, ect.
            path_stream << "./" << out.statistics.settings.uniqueStr.val << ".console_out";
            path_string = path_stream.str();
        */
        path_stream << "./" << directory_name << "/replica_log.txt";
        path_string = path_stream.str();
    } else {
        printf("Creating directory : %s\n", directory_name.c_str());
        const int dir_err = mkdir(directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory! Writing in pwd\n");
            path_stream << "./" << directory_name << "_replica_log.txt";
            path_string = path_stream.str();
        } else {
            path_stream << "./" << directory_name << "/" << "replica_log.txt";
            path_string = path_stream.str();
        }
    }
  }
  const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
  //if (useReplicaExchange && MASTER(cr)) {
  //   repl_ex = init_replica_exchange(fplog, cr->ms, state_global, ir,
    //                                    replExParams);
  //}
  
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

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintTime();
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
