#ifndef _GETINTEGRATORSOLVER_H_
#define _GETINTEGRATORSOLVER_H_

// returns the string corresponding to the selected integrator solver
// "solver" must be pre-allocated
// result: PARDISO, SPOOLES or PCG
void GetIntegratorSolver(char * solver);
 
#endif

