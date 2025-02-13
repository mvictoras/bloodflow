/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "palabos3D.h"
#include "palabos3D.hh"

#include "ibm3D.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>

//#include "opts.h"

#include "mpi.h"
#include "lammps.h"
#include "domain.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"
#include "update.h"
#include "neighbor.h"

#include "latticeDecomposition.h"
#ifdef ENABLE_ASCENT

#include <random>
#include <sstream>
#include <conduit.hpp>
#include "insitu/AscentBridge.h"

#endif

using namespace plb;
using namespace std;




typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu
    T deltaP = - (alpha * nu);
    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }

    sum *= ((T)4 * alpha * a *a /std::pow(pi,(T)3));
    sum += (alpha / (T)2 * (x * x - a*a / (T)4));
    
    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class ShearTopVelocity {
public:
    ShearTopVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T(); 
        u[1] = T(); 
        u[2] = parameters.getLatticeU();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class ShearBottomVelocity {
public:
    ShearBottomVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class WallDomain3D: public DomainFunctional3D{
public:
	WallDomain3D(Array<plint,2> coor_, plint r_): coor(coor_),r(r_){}
	virtual bool operator ()(plint iX, plint iY, plint iZ) const {
          plint result = 0;
              T rSqr = util::sqr(iX-coor[0]) + util::sqr(iY-coor[1]);
              if (rSqr >= r*r) result = 1;
          return result; 
        }
	virtual WallDomain3D<T> * clone() const {
	  return new WallDomain3D<T>(*this);
	}
	private:
        Array<plint,2> coor;
        plint r; 
};

void squarePoiseuilleSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D top    = Box3D(0,    nx-1, ny-1, ny-1, 0, nz-1);
    Box3D bottom = Box3D(0,    nx-1, 0,    0,    0, nz-1);
    
    //Box3D inlet  = Box3D(0,    nx-1, 1,    ny-2, 0,    0);
    //Box3D outlet = Box3D(0,    nx-1, 1,    ny-2, nz-1, nz-1);
    
    Box3D left   = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D right  = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);
    // shear flow top bottom surface
    /*
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow );
    
    setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet, boundary::outflow );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::outflow );
    */
    // channel flow
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet);
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet);
/*    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );
    
    //setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    //setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    
    setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
 */
 	Array<plint,2> coor(nx/2,ny/2);
	plint r = nx/2;
     defineDynamics(lattice,lattice.getBoundingBox(),new WallDomain3D<plint>(coor,r), new BounceBack<T,DESCRIPTOR>);	

    //initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,0.0,0.0));

    lattice.initialize();
}

/// This functional defines a data processor for the instantiation
///   of bounce-back nodes following the cosine shape clots
// 1/29/2025: clot profile: y = A*cos(2*pi/L*(z-zc))
template <typename T, template <typename U> class Descriptor>
class DynamicBoundaryFunctional : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    //DynamicBoundaryFunctional(plint xc_, plint yc_, plint radius_, T rn_, plint it_, plint zl_, plint clotLoc_) : xc(xc_), yc(yc_), radius(radius_), rn(rn_), it(it_), zl(zl_), clotLoc(clotLoc_) { }
    DynamicBoundaryFunctional(plint A_, plint L_, plint Zc_, T omega_) : A(A_), L(L_), Zc(Zc_), omega(omega_) { }
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
    {
	    /*
        // clotLoc - location of clot 
        //BounceBackNodes<T> bbDomain(N, radius);
        // limiting the iteration of Z for greater efficiency, and unlocking the Z location of clot
        int zrange = 10; // NEEDS FINISHING - Nazariy 12/21/2022
        int zmin = domain.z0 + clotLoc - zrange/2; //default zmin
        int zmax = domain.z0 + clotLoc + zrange/2; //default zmax 
        if(clotLoc < (domain.z0 + zrange/2)){ // if clot location is before the start of the domain + center of specified range 
            zmin = domain.z0;
            zmax = domain.z0+zrange;
        } 
        else if(clotLoc < (domain.z0 + zrange/2)){ // if clot location is before the start of the domain + center of specified range 
            zmin = domain.z0;
            zmax = domain.z0+zrange;
        } 
        else if(clotLoc > domain.z1){
            zmin = domain.z1;
            zmax = domain.z1-zrange;
        }*/ 

        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
           for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
               for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {  // was "for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {"
			T clotHeight = A*cos(2*pi*(iZ-Zc)/L);
			if (iY < clotHeight){
                        	lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
			}else{
				lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega));
			}
                  /*  T dist = (iY + relativeOffset.y - yc)*(iY + relativeOffset.y - yc) + 
                             (iX + relativeOffset.x - xc)*(iX + relativeOffset.x - xc); //XXXX change later - Nazariy
                    T xscale = .5 * rn * it; // scaling factor for x (how wide) replace with a slider later.
                    T yscale = .5 * radius * sin(.0314*it);  // scaling factor for y (how close to center) dependent on iteration
                    if(yscale < 0){
                        T yscale = -.5 * radius * sin(.0314*it);
                    }
                    // the representation is y = e ^ (-x^2), but the example is set up in Z direction so here 
                    // the x stands for domain in Z i.e. values from Z axis and y stands for value being subtracted from the radius
                    int cntZ = (iZ+relativeOffset.z-(zl / 2))/xscale; // position of Z with resepect to center of Z        
                    T modRad = yscale * exp(-cntZ*cntZ) + 1; // radius modification parameter
                    if (dist > (radius - modRad)*(radius - modRad)) {
                        lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, DESCRIPTOR>);
                    }
                    else{
                        lattice.attributeDynamics(iX, iY, iZ, new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1./0.95));
                    }*/
                }
            }
        }
    }
    
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dataStructure;
    }
    virtual DynamicBoundaryFunctional<T, Descriptor> *clone() const
    {
        return new DynamicBoundaryFunctional<T, Descriptor>(*this);
    }

private:
   /* plint xc,yc,it,zl,clotLoc; // zl added by Nazariy 7/19 // clotLoc added by Nazariy 12/20/2022
    plint radius;
    T rn; //XXXX added by Nazariy 7/12
	  */
    plint A, L, Zc; //A: amplitude, L: clot size, Zc: the center of clot in z. Assuming flow is in z direction
    T omega;
};

// template <typename T, class Descriptor<typename U>>
// void DynamicBoundaryFunctional<T, Descriptor<U>>::modifyClotLocation(int newClotLocation)
// {
//     // clotLoc = newClotLocation;
// };

/// Automatic instantiation of the bounce-back nodes for the boundary,
///   using a data processor.
void createDynamicBoundaryFromDataProcessor(
    MultiBlockLattice3D<T, DESCRIPTOR> &lattice, plint A, plint L, plint Zc, T omega) // removed 
{
//    plint it = 0; // change back to 0 when finished debugging the bidirectional stu0ff and lammps domain stuff. 
    applyProcessingFunctional(
        new DynamicBoundaryFunctional<T, DESCRIPTOR>(A,L,Zc, omega), lattice.getBoundingBox(),
        lattice);
}

void modifyDynamicBoundaryFromDataProcessor(
    MultiBlockLattice3D<T, DESCRIPTOR> &lattice, plint A, plint L, plint Zc, T omega)
{
    applyProcessingFunctional(
        new DynamicBoundaryFunctional<T, DESCRIPTOR>(A,L,Zc,omega), lattice.getBoundingBox(),
        lattice);
}

T computeRMSerror ( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters )
{
    MultiTensorField3D<T,3> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   SquarePoiseuilleVelocity<T>(parameters, NMAX) );
    MultiTensorField3D<T,3> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
           std::sqrt( computeAverage( *computeNormSqr(
                          *subtract(analyticalVelocity, numericalVelocity)
                     ) ) );
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}
//**************************************Connor Changed 2/17/22
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              Box3D domain, plint iter)
{
    //T dx = parameters.getDeltaX();
    //T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1);
    //vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice,domain), "velocity", 1);
    //vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

LammpsWrapper* wrapper;

#ifdef ENABLE_ASCENT
double random_double(double min, double max) {
    static std::random_device rd;  // Seed for the random number engine
    static std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

void addRandomCell(conduit::Node &params, conduit::Node &output) {
    std::cout << "Inserting a new RBC" << std::endl;
    double pt[] = {random_double(0.0, 20.0), random_double(0.0, 20.0), random_double(0.0, 20.0)};
    std::ostringstream fixDepositString;
    fixDepositString << "fix 3 cells deposit 1 0 1 12345 mol singleRBC region RBC_zone id max gaussian "
                     << pt[0] << " " << pt[1] << " " << pt[2] << " 10 near 2 " << std::endl;

    std::string commandStr = fixDepositString.str();
    std::vector<char> commandVec(commandStr.begin(), commandStr.end());
    commandVec.push_back('\0');

    wrapper->execCommand(commandVec.data());
    fixDepositString.str("");
}
#endif

//**************************************
int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    
/*
    if (argc != 2) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : N\n";
        exit(1);
    }*/

    //const plint N = atoi(argv[1]);
    const plint N = 1;// atoi(argv[1]);
    const T Re = 5e-3;
    //const plint Nref = 50;
    //const T uMaxRef = 0.01;
    const T uMax = 0.00075;//uMaxRef /(T)N * (T)Nref; // Needed to avoid compressibility errors

    //using namespace opts;

#ifdef ENABLE_ASCENT
    AscentBridge::getInstance().Initialize(global::mpi().getGlobalCommunicator());
    ascent::register_callback("addRandomCell", addRandomCell);
#endif

    /*Options ops(argc, argv);
    ops
    #ifdef ENABLE_SENSEI
    >> Option('f', "config", config_file, "Sensei analysis configuration xml (required)")
    #endif
    */
    
    const T maxT = atoi(argv[2]);//6.6e4; //(T)0.01;
    plint iSave = atoi(argv[3]);//10;//2000;//10;
    //plint iCheck = 10*iSave;
    

    wrapper = new LammpsWrapper(argv,global::mpi().getGlobalCommunicator()); // replaced with MPI_COMM_WORLD
    // LammpsWrapper wrapper(argv,MPI_COMM_WORLD);

    char * inlmp = argv[1];

    wrapper->execFile(inlmp);

    LAMMPS_NS::Domain *domain = wrapper->lmp->domain;
    // Get the simulation box boundaries
    double xlo = domain->boxlo[0];
    double xhi = domain->boxhi[0];
    double ylo = domain->boxlo[1];
    double yhi = domain->boxhi[1];
    double zlo = domain->boxlo[2];
    double zhi = domain->boxhi[2];

    // Print the boundaries
    printf("Domain boundaries:\n");
    printf("xlo: %f, xhi: %f\n", xlo, xhi);
    printf("ylo: %f, yhi: %f\n", ylo, yhi);
    printf("zlo: %f, zhi: %f\n", zlo, zhi);


    const int nx = static_cast<int>(xhi - xlo);
    const int ny = static_cast<int>(yhi - ylo);
    const int nz = static_cast<int>(zhi - zlo);
    IncomprFlowParam<T> parameters(
            uMax,
            Re,
            N,
            nx,        // lx
            ny,        // ly
            nz         // lz
    );
    writeLogFile(parameters, "3D square Poiseuille");
    
    
   
    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    plint mysize = global::mpi().getSize();
    plint localdomain[mysize][6]; //First index: Rank value Second Index: extents (0: xlo 1: xhi 2: ylo 3: yhi 4: zlo 5: zhi)

    pcout<<"Nx,Ny,Nz "<<parameters.getNx()<<" "<<parameters.getNy()<<" "<<parameters.getNz()<<endl;
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper->lmp, localdomain); //XXX sending a double array to store extents of each processor for palabos data (Connor Murphy 2/19/22)

    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 3;


    MultiBlockManagement3D management = MultiBlockManagement3D(blockStructure, threadAttribution, envelopeWidth);//XXX New Change
    
    

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (management,
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    
    //********************************************* XXX Possible alternative for passing local extents to computeVelocity function
    /*
    ThreadAttribution const & orgThreadAttribution = management.getThreadAttribution();
    std::vector<plint> localBlocks = blockStructure.getLocalBlocks(orgThreadAttribution);
    std::map<plint, Box3D> bulksMap = blockStructure.getBulks();
    std::vector<Box3D> bulks;
    plint blockId;
    Box3D bulk;
    auto it=bulksMap.begin();
    cout << " bulkSmap size " << bulksMap.size() << endl;
    for(;it!=bulksMap.end();++it)
    {
        bulk = it->second;
        blockId = it->first;
        cout<<"block id= " << blockId<< " proc " << myrank <<" bulk: Nx " << bulk.getNx() << " Ny " << bulk.getNy() << " Nz " << bulk.getNz() << endl;
    }
    */
    //**********************************************
    
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    //pcout<<"51 works"<<endl;

    /*
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );*/

    // Use periodic boundary conditions.
    lattice.periodicity().toggle(2,true);
   

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);

    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
    //coupling between lammps and palabos
    Array<T,3> force(0,0,1e-6);
    setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
    //LAMMPS

   /* plint xc, yc, radius, iterationCAS, zLength, clotLoc;

    xc = nx/2; // X center 
    yc = ny/2; // Y center
    radius = nx / 2; // radius
    T radiusNorm = radius/maxT; // double radius/max iterations
    iterationCAS = 10; // iterations for collideAndStream in the loop 
    zLength = parameters.getNz(); // because domain.z0 gives local value pass this instead.
    clotLoc = zLength;  
  */ 
    plint A, L, Zc;
    A = ny/10;
    L = nz/4;
    Zc = nz/2;  
    //createDynamicBoundaryFromDataProcessor(lattice, A,L,Zc,parameters.getOmega()); // added by NT 7/18/2022

    for (plint iT=0;iT<1e2;iT++){ //(plint iT=0;iT<4e3;iT++){
        lattice.collideAndStream();
    }

    long time = 0; 
 
    // for (plint iT=0;iT<4e3;iT++){
    //     lattice.collideAndStream();
    // }
    T timeduration = T();
    global::timer("mainloop").start();
    plint nlocal;
    long ntimestep;
    int nanglelist;
    int nghost;
    double **x;
    double **v;
   
    int **anglelist;
    // Array<double,3> center(0.,0.,0.);
    //  std::vector<std::double> center;

    plint myrank = global::mpi().getRank();
    MultiTensorField3D<T,3> vel(lattice);
    MultiTensorField3D<double, 3> vort(lattice);
    MultiScalarField3D<double> velNorm(lattice);

    //TensorField3D<T,3> velocityArray = vel.getComponent(myrank);
    //TensorField3D<T,3> vorticityArray = vort.getComponent(myrank);
    //ScalarField3D<T> velocityNormArray = velNorm.getComponent(myrank);
    int test2 = 1;
    int t_count = 0;
    float t = 0.;
    std::stringstream fixDepositString;
    int fixID = 3;
  

    for (plint iT=0; iT<maxT; ++iT) {
        
        // if (iT%iSave ==0 && iT >0){
        //     wrapper->execCommand("fix 3 cells deposit 1 0 1 12345 mol singleRBC region RBC_zone id max gaussian 15 15 10 10 near 4");
        //     // wrapper->execFile("in.deposit");
        //     // wrapper->execCommand("dump 1 cells xyz 1 dump.rbc.xyz");
        // }
        // lammps to calculate force
        // wrapper->execCommand("run 1"); // pre no post no");
        wrapper->execCommand("run 1 pre no post no");

        //Some values are dynamically changing
        nlocal = wrapper->lmp->atom->nlocal;
        ntimestep = wrapper->lmp->update->ntimestep;
        nanglelist = wrapper->lmp->neighbor->nanglelist;
        nghost = wrapper->lmp->atom->nghost;
        x = wrapper->lmp->atom->x;
        v = wrapper->lmp->atom->v;
        anglelist = wrapper->lmp->neighbor->anglelist;
        

        //*************************************
        vel = *computeVelocity(lattice,lattice.getBoundingBox());
        vort = *computeVorticity(vel);
        velNorm = *computeVelocityNorm(lattice,lattice.getBoundingBox());

	    TensorField3D<T,3> velocityArray = vel.getComponent(myrank);
   	    TensorField3D<T,3> vorticityArray = vort.getComponent(myrank);
      	ScalarField3D<T> velocityNormArray = velNorm.getComponent(myrank);
      
        //Box3D domain = Box3D(localdomain[myrank][0]-envelopeWidth,localdomain[myrank][1]+envelopeWidth,localdomain[myrank][2]-envelopeWidth,localdomain[myrank][3]+envelopeWidth,localdomain[myrank][4]-envelopeWidth,localdomain[myrank][5]+envelopeWidth);
        Box3D domain = Box3D(localdomain[myrank][0],localdomain[myrank][1],localdomain[myrank][2],localdomain[myrank][3],localdomain[myrank][4],localdomain[myrank][5]);
        //*************************************
        
        cout<<"Rank: " << myrank <<" local domain Extents: x: " <<domain.x0 << " " << domain.x1 << " y: " << domain.y0 <<" "<<domain.y1<< " z "<<domain.z0<<" "<<domain.z1<<endl;
        cout<<"Rank: " << myrank <<" Vorticity Extents: " <<vorticityArray.getNx() << " " << vorticityArray.getNy() << " " << vorticityArray.getNz()<<endl;
        //cout<<"Rank: " << myrank <<" Velocity Extents: " <<velocityArray.getNx() << " " << velocityArray.getNy() << " " << velocityArray.getNz()<<endl;
        //cout<<"Rank: " << myrank <<" Velocity Norm Extents: " <<velocityNormArray.getNx() << " " << velocityNormArray.getNy() << " " << velocityNormArray.getNz()<<endl;
#ifdef ENABLE_ASCENT
        if (iT%(iSave) ==0 && iT >0){

            AscentBridge::getInstance().Publish(x, v, ntimestep, nghost ,nlocal, anglelist, nanglelist,
                                velocityArray, vorticityArray, velocityNormArray, 
                                nx, ny, nz, domain, envelopeWidth);
            if(iT == 5) {
                std::cout << "Inserting a new RBC" << std::endl;
                int pt[] = {10, 10, 10};
                fixDepositString << "fix 3 cells deposit 1 0 1 12345 mol singleRBC region RBC_zone id max gaussian "<<pt[0]<<" "<<pt[1]<<" "<< pt[2] << " 1 near 2 "<<endl;
                std::cout << "Deposit string: " << fixDepositString.str() << std::endl;
                //fix 3 cells deposit 1 0 1 12345 mol singleRBC region RBC_zone id max gaussian 10 10 10 10 near 2 # vz 10 20 
                wrapper->execCommand(fixDepositString);
                //wrapper->execCommand("fix 3 cells deposit 1 0 1 12345 mol singleRBC region RBC_zone id max gaussian 10 10 5 10 near 2 ");// this is working, 7/6/2023 TISHCHENKO
                fixDepositString.str("");
            }
        }
#endif        

        // Clear and spread fluid force
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
        ////------------ classical ibm coupling-------------//
        spreadForce3D(lattice,*wrapper);
        ///--------------redefine a new domain--------------// NT 12/20
        if(iT == 6) {
            std::cout << "New clot" << std::endl;
            //clotLoc = 1;
            //radius = 15;
            //modifyDynamicBoundaryFromDataProcessor(lattice, A, L, Zc, parameters.getOmega());
        }
        ////// Lattice Boltzmann iteration step.

        // for(int iteration=0; iteration<iterationCAS; iteration++){
        //     lattice.collideAndStream();
        // }
        // if the modifyDynamicBoundaryFromDataProcessor is used the above is the "proper" way of doing it 12/21
        
        lattice.collideAndStream();
        ////// Interpolate and update solid position
        interpolateVelocity3D(lattice,*wrapper);
        //-----force FSI ibm coupling-------------//
        //forceCoupling3D(lattice,wrapper);
        //lattice.collideAndStream();
        //writeVTK(lattice, domainBox, iT);
	
    }
    wrapper->execCommand("dump 2 cells xyz 1 dump2.rbc.xyz");
    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
#ifdef ENABLE_ASCENT
    AscentBridge::getInstance().Finalize();
#endif
}


