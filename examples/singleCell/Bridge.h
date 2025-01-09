#pragma once
#include <mpi.h>
#include <string>
#include "palabos3D.h"
#include "palabos3D.hh"

#include <ascent.hpp>

using namespace plb; 

class Bridge
{
  protected:
    Bridge();  // Constructor

    // Private static instance pointer
    static std::unique_ptr<Bridge> singleton_;

  public:

    /**
     * Singletons should not be cloneable.
     */
    Bridge(const Bridge &other) = delete;
    /**
     * Singletons should not be assignable.
     */
    Bridge& operator=(const Bridge &) = delete;

    ~Bridge(); // Destructor

    static Bridge& getInstance();  // Method to get the singleton instance

    void Initialize(MPI_Comm world);
    void Publish();
    void Publish(double **x, double **v, long ntimestep, int nghost, 
                int nlocal, int **anglelist, int nanglelist,
              TensorField3D<double, 3> velocityDoubleArray, 
              TensorField3D<double, 3> vorticityDoubleArray, 
              ScalarField3D<double> velocityNormDoubleArray,
              int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth); 
    void Execute(long ntimestep); 
    void Finalize();

  private:
    ascent::Ascent mAscent;
};

