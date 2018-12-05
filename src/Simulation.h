/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SIMULATION_H
#define SIMULATION_H

//Member vars
#include "CPUSide.h"
#include "System.h"
#include "StaticVals.h"
#include "BasicTypes.h"
//#include "ReplicaState.h"
#include "ReplEx.h"


class Simulation
{
public:
  explicit Simulation(char const*const configFileName);
  ~Simulation();

  void RunSimulation(void);

#ifndef NDEBUG
  void RunningCheck(const uint step);
#endif

private:
  StaticVals * staticValues;
  System * system;
  CPUSide * cpu;
  ulong totalSteps;
  ReplicaExchangeParameters replExParams;
  FILE * fplog;
  int nnodes;
  int nodeid;
  bool bDoReplEx, bLastStep, bExchanged;
  ReplicaState * stateLocal;
  ReplicaState * stateGlobal;
  gmx_repl_ex * replEx;
};

#endif /*SIMULATION_H*/
