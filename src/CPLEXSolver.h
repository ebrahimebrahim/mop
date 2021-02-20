

#ifndef CPLEXSOLVER_H
#define CPLEXSOLVER_H

#include "LPSolver.h"

#include <ilcplex/cplex.h>
#include <stdlib.h>



class CPLEXSolver : public LPSolver{

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

    int nVariables;
    int solstat;
    long iCount;
    double objValue;


  public:

    CPLEXSolver( long optimizer = CPX_ALG_AUTOMATIC) {

      nVariables=0;
      int status=-1;
      env = CPXopenCPLEX (&status);
        
      CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 0 );

#ifdef VERBOSE
      CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
      CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 1 );
#endif
      
      CPXsetdblparam (env, CPX_PARAM_NETEPOPT, 1e-10);
      CPXsetintparam (env, CPX_PARAM_LPMETHOD, optimizer);

      solstat = CPX_STAT_ABORT_USER;
    };


    ~CPLEXSolver(){
      deleteLP();
      CPXcloseCPLEX (&env);
    };



   virtual void solveLP(){

     int status;
     CPXLPptr prob = CPXcreateprob (env, &status, "ot");

     CPXnewcols(env, prob, coeff.size(), coeff.data(), colLB.data(), colUB.data(), NULL, NULL);

     CPXnewrows(env, prob, mass.size(), mass.data(), NULL, NULL, NULL);

     std::vector<double> ones(sInd.size(),  1);
     std::vector<double> nones(tInd.size(), -1);
     std::vector<int> columnIndex(tInd.size());
     for(int i =0; i<columnIndex.size(); i++){
       columnIndex[i] = i;
     }
     CPXchgcoeflist( env, prob, columnIndex.size(), sInd.data(), columnIndex.data(),
         ones.data() );
     CPXchgcoeflist( env, prob, columnIndex.size(), tInd.data(), columnIndex.data(),
         nones.data() );
     CPXcopybase (env, prob, colStatus.data(), rowStatus.data());


     status = CPXlpopt (env, prob);

     status = CPXsolution (env, prob, &solstat, &objValue, primal.data(),
         dual.data(), NULL, NULL);

#ifdef VERBOSE
     std::cout << "CPX status: " << status << std::endl;
     std::cout << "CPX solution status: " << solstat << std::endl;
#endif

     iCount = CPXgetitcnt(env, prob);

     CPXgetbase (env, prob, colStatus.data(), rowStatus.data());

     CPXfreeprob(env, &prob);

#ifdef VERBOSE
     std::cout << "Obj: " << objValue << std::endl;
     std::cout << "#iter: " <<iCount << std::endl << std::endl;
#endif
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
     nVariables = nSource + nTarget;
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



   virtual void setRowBounds(long row, double m){
     mass[row] = m;
   };


   virtual double getRowBounds(long row){
     return mass[row];
   };

   virtual void setColumnObjective(long col, double cost){
     coeff[col] = cost;
   };


   virtual void setColumnCoefficients(long col, long s, long t){
     tInd[col] = t;
     sInd[col] = s;
   };


   virtual void setColumnBounds(long col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };

   virtual void setColumnBoundsLower(long col, double lb){
     colLB[col] = lb;
     colUB[col] = CPX_INFBOUND;
   };




   virtual long getColumn(long col, long *ind, double *val){
     ind[0] = sInd[col];
     ind[1] = tInd[col];
     val[0] = colLB[col];
     val[1] = colUB[col];
     return 2;
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
       case LPSolver::END:
       case LPSolver::UNKNOWN:
         return -1;
     }
     return -1;

   };



};


#endif
