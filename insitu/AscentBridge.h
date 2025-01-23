#pragma once
#include <mpi.h>
#include <string>
#include "palabos3D.h"
#include "palabos3D.hh"
#include <mutex>
#include <ascent.hpp>

using namespace plb; 

class AscentBridge
{
  protected:
    AscentBridge();  // Constructor

    // Private static instance pointer
    static std::unique_ptr<AscentBridge> singleton_;

  public:

    /**
     * Singletons should not be cloneable.
     */
    AscentBridge(const AscentBridge &other) = delete;
    /**
     * Singletons should not be assignable.
     */
    AscentBridge& operator=(const AscentBridge &) = delete;

    ~AscentBridge(); // Destructor

    static AscentBridge& getInstance();  // Method to get the singleton instance

    void Initialize(MPI_Comm world);
    
    void Publish(double **x, double **v, long ntimestep, int nghost, 
                int nlocal, int **anglelist, int nanglelist,
              TensorField3D<double, 3> velocityDoubleArray, 
              TensorField3D<double, 3> vorticityDoubleArray, 
              ScalarField3D<double> velocityNormDoubleArray,
              int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth); 
    
    void Finalize();

  private:
    ascent::Ascent mAscent;
};

