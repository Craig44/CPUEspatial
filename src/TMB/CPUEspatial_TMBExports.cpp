// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_CPUEspatial_TMBExports
#include <TMB.hpp>
#include "HelperFuns.hpp"
#include "SpatialTemporalCPUE.hpp"
#include "SpatialTemporalCPUENN.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "SpatialTemporalCPUE") {
    return SpatialTemporalCPUE(this);
  } else if(model == "SpatialTemporalCPUENN") {
    return SpatialTemporalCPUENN(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}