

#ifndef CPLEXNETWORKSOLVER_H
#define CPLEXNETWORKSOLVER_H

#include "LPSolver.h"

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <iostream>

#include <vector>
#include <math.h>
#include <limits>

class CPLEXNetworkSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;



    CPXENVptr env;

    std::vector<int> sInd;
    std::vector<int> tInd;

    std::vector<double> coeff;
    std::vector<double> mass;
    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<double> colLB;
    std::vector<double> colUB;

    std::vector<int> colStatus;
    std::vector<int> rowStatus;

    int solstat;
    long iCount;
    double objValue;
    long ns;
    long nt;


  public:

    CPLEXNetworkSolver() {

        int status=-1;
        env = CPXopenCPLEX( &status );

        CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 0 );

#ifdef VERBOSE
        CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
        CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 1 );
#endif
        CPXsetdblparam( env, CPX_PARAM_NETEPOPT, 1e-5 );
        CPXsetdblparam( env, CPX_PARAM_NETEPRHS, 1e-7 );
        CPXsetintparam (env, CPX_PARAM_ADVIND, 2);
        CPXsetdblparam( env, CPX_PARAM_EPPER, 1e-7 );

        ns = 0;
        solstat = CPX_STAT_ABORT_USER;
    };


    ~CPLEXNetworkSolver(){
      deleteLP();
      CPXcloseCPLEX (&env);
    };



   virtual void solveLP(){

     int status;
     CPXNETptr prob = CPXNETcreateprob (env, &status, "netex1");

  
     //Add nodes
     CPXNETaddnodes( env, prob, mass.size(), mass.data(), NULL);

     CPXNETaddarcs( env, prob, sInd.size(), sInd.data(), tInd.data(),
                      colLB.data(), colUB.data(), coeff.data(), NULL );

     CPXNETchgobjsen( env, prob, CPX_MIN );
     CPXNETcopybase( env, prob, colStatus.data(), rowStatus.data() );


     status = CPXNETprimopt (env, prob);

     status = CPXNETsolution( env, prob, &solstat, &objValue, primal.data(),
         dual.data(), NULL, NULL);


     if(solstat != CPX_STAT_OPTIMAL ){
       //std::cout << "Arg" << std::endl;
     }

     iCount = CPXNETgetitcnt( env, prob);

     CPXNETgetbase( env, prob, colStatus.data(), rowStatus.data() );


     CPXNETfreeprob( env, &prob );
   
   };



   virtual bool isOptimal(){
     return solstat == CPX_STAT_OPTIMAL;
   };



   virtual double getObjectiveValue(){
     return objValue;
   };


   virtual long getIterationCount(){
     return iCount;
   };


   virtual long getNumberOfRows(){
     return dual.size();
   };


   virtual long getNumberOfColumns(){
     return primal.size();
   };

   virtual void setupStandardBasis(){
     for(long i= 0; i< getNumberOfColumns(); i++){
       colStatus[i] = CPX_AT_LOWER;
     }
     for(long i= 0; i< getNumberOfRows(); i++){
       rowStatus[i] =  CPX_BASIC;
     }
   };



   virtual void createLP(long nSource, long nTarget){
     deleteLP();
     ns=nSource;
     nt=nTarget;

     rowStatus.resize( 1, CPX_BASIC );
     dual.resize( 1, 0);

     primal.resize( ns+nt, 0);
     colStatus.resize( ns+nt, CPX_AT_LOWER );

   };

   virtual void addColumns(long n){
     sInd.resize( sInd.size() + n, -1  );
     tInd.resize( tInd.size() + n, -1 );
     coeff.resize( coeff.size() + n, 0 );
     primal.resize( primal.size() + n, 0 );
     colStatus.resize( colStatus.size() + n, CPX_AT_LOWER );
     colLB.resize( colLB.size() + n, 0 );
     colUB.resize( colUB.size() + n, CPX_INFBOUND );
   };



   virtual void addRows(long n){
     mass.resize( mass.size() + n , 0);
     rowStatus.resize( rowStatus.size() + n, CPX_BASIC );
     dual.resize( dual.size() + n, 0);
   };



   virtual double getRowDual(long row){
     return dual[row];
   };



   virtual double getColumnPrimal(long col){
     return primal[col];
   };


   virtual void setRowBounds(long i, double m){
     mass[i] = m;
   };


   virtual double getRowBounds(long i){
     return mass[i];
   };

   virtual void setColumnBounds( long col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };

   virtual void setColumnBoundsLower(long col, double lb){
     colLB[col] = lb;
     colUB[col] = CPX_INFBOUND;
   };


   virtual Status getColumnStatus(long col){
     //std::cout << "get: " << colStatus[col] << std::endl;
     return convertFromCPLEX( colStatus[col] );

   };

   virtual Status getRowStatus(long row){
     return convertFromCPLEX( rowStatus[row] );
   };

   virtual void setColumnStatus(long col, Status s){
     //std::cout << s << " -> ";
     //std::cout << convertToCPLEX(s) << std::endl;
     colStatus[col] = convertToCPLEX(s);
   };


   virtual void setRowStatus(long row, Status s){
     if(s == LPSolver::BASIC){
       rowStatus[row] = CPX_BASIC;
     }
     else{
       rowStatus[row] = CPX_AT_LOWER;
     }
   };


   virtual void setColumnObjective(long i, double cost){
     coeff[i] = cost;
   };

   virtual void setColumnCoefficients(long col, long s, long t ){
     tInd[col] = t;
     sInd[col] = s;
   };

   virtual long getColumn(long col, long *ind, double *val){
     ind[0] = sInd[col];
     ind[1] = tInd[col];
     val[0] = 1;
     val[1] = -1;
     return 2;
   };



  private:



   virtual void deleteLP(){
     sInd.clear();
     tInd.clear();
     coeff.clear();
     mass.clear();
     dual.clear();
     primal.clear();
     rowStatus.clear();
     colStatus.clear();
     colLB.clear();
     colUB.clear();
     ns = 0;
   };



   Status convertFromCPLEX(long s){

     switch(s){
       case CPX_BASIC:
         return LPSolver::BASIC;
       case CPX_AT_UPPER:
         return LPSolver::UPPER;
       case CPX_AT_LOWER:
         return LPSolver::LOWER;
       case CPX_FREE_SUPER:
         return LPSolver::FREE;
     }
     return LPSolver::END;
   };




   long convertToCPLEX(Status s){

     switch(s){
       case LPSolver::BASIC:
         return CPX_BASIC;
       case LPSolver::UPPER:
         return CPX_AT_UPPER;
       case LPSolver::LOWER:
         return CPX_AT_LOWER;
       case LPSolver::SUPERBASIC:
       case LPSolver::FREE:
         return CPX_FREE_SUPER;
       case LPSolver::FIXED:
       case LPSolver::INF:
       case LPSolver::UNKNOWN:
       case LPSolver::END:
         return -1;
     }
     return -1;

   };



};


#endif
