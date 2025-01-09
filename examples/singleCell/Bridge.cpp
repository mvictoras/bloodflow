#include "Bridge.h"
#include <iostream>

using namespace std;
using namespace plb; 

std::unique_ptr<Bridge> Bridge::singleton_;

Bridge::Bridge() {
}

Bridge::~Bridge() {
}

// Define the static method to get the singleton instance
Bridge& Bridge::getInstance() {
  static std::once_flag flag;
  std::call_once(flag, [](){
      singleton_.reset(new Bridge());
      });
  return *singleton_;
}

void Bridge::Initialize(MPI_Comm world) {
  conduit::Node ascent_opts;
  ascent_opts["mpi_comm"] = MPI_Comm_c2f(world);
  mAscent.open(ascent_opts);
}
void Bridge::Publish  () {
  cout << "Publishing" << endl;
}

void Bridge::Publish  (double **x, double **v, long ntimestep, int nghost, 
                      int nlocal, int **anglelist, int nanglelist, 
                      TensorField3D<double, 3> velocityArray,
                      TensorField3D<double, 3> vorticityArray,
                      ScalarField3D<double> velocityNormArray,
                      int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth) {

  
  //GlobalDataAdaptor->AddLAMMPSData(x, ntimestep, nghost, nlocal, anglelist, nanglelist);

  conduit::Node mesh;
  //conduit::Node &volume_coords = mesh["coordsets/volume_coords"];
  
  const int num_tris = nanglelist;
  const int conn_len = num_tris * 3;

/*
  data["coordsets/coords/type"] = "explicit";
  
  data["coordsets/coords/values/x"].set(conduit::DataType::int32(num_tris));
  data["coordsets/coords/values/y"].set(conduit::DataType::int32(num_tris));
  data["coordsets/coords/values/z"].set(conduit::DataType::int32(num_tris));

  data["topologies/mesh/type"]           = "unstructured";
  data["topologies/mesh/coordset"]       = "coords";
  data["topologies/mesh/elements/shape"] = "tri";
  data["topologies/mesh/elements/connectivity"].set(conduit::DataType::int32(conn_len));

  // Get pointers to fill in
  int32_t *x_ptr = data["coordsets/coords/values/x"].value();
  int32_t *y_ptr = data["coordsets/coords/values/y"].value();
  int32_t *z_ptr = data["coordsets/coords/values/z"].value();
  
  // Fill the connectivity array
  for(int i = 0; i < num_tris; i++) {
      x_ptr[i] = anglelist[i][0];
      y_ptr[i] = anglelist[i][1];
      z_ptr[i] = anglelist[i][2];
  }
*/
vector<double> velocity_x, velocity_y, velocity_z;
int nlx = velocityArray.getNx(); 
int nly = velocityArray.getNy();
int nlz = velocityArray.getNz();

std::cout << nlx << " " << nly << " " << nlz << std::endl;
mesh["coordsets/volume_coords/type"] = "uniform";
mesh["coordsets/volume_coords/dims/i"] = nlx;
mesh["coordsets/volume_coords/dims/j"] = nly;
mesh["coordsets/volume_coords/dims/k"] = nlz;
//mesh["coordsets/volume_coords/origin/x"] = 0.0;
//mesh["coordsets/volume_coords/origin/y"] = 0.0;
//mesh["coordsets/volume_coords/origin/z"] = 0.0;
//mesh["coordsets/volume_coords/spacing/dx"] = 1.0;
//mesh["coordsets/volume_coords/spacing/dy"] = 1.0;
//mesh["coordsets/volume_coords/spacing/dz"] = 1.0;
mesh["topologies/volume_topo/type"] = "uniform";
mesh["topologies/volume_topo/coordset"] = "volume_coords";

  Array<double,3> vel(0.,0.,0.); //jifu 5/31/2022
  Array<double,3> vor(0.,0.,0.);
  for (int k=0; k<nlz; k++)
  {
    for (int j=0; j<nly; j++)
    {
     for (int i=0; i<nlx; i++)
      {
        vel = velocityArray.get(i,j,k); //jifu 5/31/2022 Array<double 3> vel ->vel
        vor = vorticityArray.get(i,j,k);
        double norm = velocityNormArray.get(i,j,k);
        int index = (j) * (nlx) + (i) + (k) * (nlx) * (nly);
        velocity_x.push_back(vel[0]);
        velocity_y.push_back(vel[1]);
        velocity_z.push_back(vel[2]);
        //velocityDoubleArray->SetTuple3(index,vel[0],vel[1],vel[2]);
        //vorticityDoubleArray->SetTuple3(index,vor[0],vor[1],vor[2]);
        //velocityNormDoubleArray->SetTuple1(index,norm);
      }
    }
  }

  std::cout << "VEL:" << velocity_x.size() << std::endl;
  // Now add velocity as separate arrays
  mesh["fields/velocity/association"] = "vertex";
  mesh["fields/velocity/topology"]    = "volume_topo";
  mesh["fields/velocity/values/x"].set(velocity_x);
  mesh["fields/velocity/values/y"].set(velocity_y);
  mesh["fields/velocity/values/z"].set(velocity_z);

    //vel["values/v"].set(vy, NX*NY*NZ);
    //vel["values/w"].set(vz, NX*NY*NZ);
  //data["fields/velocity_x/association"]  = "vertex";
  //data["fields/velocity_x/topology"]     = "mesh";
  //data["fields/velocity_x/values"].set(velocity_x);

/*
  vector<double> velocity_x;

  int nlx = velocityArray.getNx(); 
  int nly = velocityArray.getNy();
  int nlz = velocityArray.getNz();

  volume_coords["dims/i"] = nlx;
  volume_coords["dims/j"] = nly;
  volume_coords["dims/k"] = nlz;

  // Optional origin and spacing
  volume_coords["origin/x"] = 0.0;
  volume_coords["origin/y"] = 0.0;
  volume_coords["origin/z"] = 0.0;
  volume_coords["spacing/dx"] = 1.0;
  volume_coords["spacing/dy"] = 1.0;
  volume_coords["spacing/dz"] = 1.0;

  conduit::Node &topo = mesh["topologies/volume_topo"];
  topo["type"] = "uniform";
  topo["coordset"] = "volume_coords";

  Array<double,3> vel(0.,0.,0.); //jifu 5/31/2022
  Array<double,3> vor(0.,0.,0.);
  for (int k=0; k<nlz; k++)
  {
    for (int j=0; j<nly; j++)
    {
     for (int i=0; i<nlx; i++)
      {
        vel = velocityArray.get(i,j,k); //jifu 5/31/2022 Array<double 3> vel ->vel
        vor = vorticityArray.get(i,j,k);
        double norm = velocityNormArray.get(i,j,k);
        int index = (j) * (nlx) + (i) + (k) * (nlx) * (nly);
        velocity_x.push_back(vel[0]);
        //velocityDoubleArray->SetTuple3(index,vel[0],vel[1],vel[2]);
        //vorticityDoubleArray->SetTuple3(index,vor[0],vor[1],vor[2]);
        //velocityNormDoubleArray->SetTuple1(index,norm);
      }
    }
  }

  std::cout << "VEL:" << velocity_x.size() << std::endl;
    // Now add velocity as separate arrays
    conduit::Node &velnode = mesh["fields/velocity_x"];
    velnode["association"] = "vertex";
    velnode["topology"]    = "volume_topo";
    velnode["values"].set(velocity_x);
    //vel["values/v"].set(vy, NX*NY*NZ);
    //vel["values/w"].set(vz, NX*NY*NZ);
*/
  //data["fields/velocity_x/association"]  = "vertex";
  //data["fields/velocity_x/topology"]     = "mesh";
  //data["fields/velocity_x/values"].set(velocity_x);


  int nvals = nlocal+nghost;
  mesh["coordsets/surface_coords/type"] = "explicit";
  mesh["coordsets/surface_coords/values/x"].set(conduit::DataType::float64(nvals));
  mesh["coordsets/surface_coords/values/y"].set(conduit::DataType::float64(nvals));
  mesh["coordsets/surface_coords/values/z"].set(conduit::DataType::float64(nvals));

  mesh["topologies/rbc/type"]           = "unstructured";
  mesh["topologies/rbc/coordset"]       = "surface_coords";
  //mesh["topologies/rbc/elements/shape"] = "point";
  mesh["topologies/rbc/elements/shape"] = "tri";
  //mesh["topologies/rbc/elements/connectivity"].set(conduit::DataType::int32(nlocal));
  mesh["topologies/rbc/elements/connectivity"].set(conduit::DataType::int32(nanglelist * 3));

  // Get pointers to fill in
  double_t *x_ptr = mesh["coordsets/surface_coords/values/x"].value();
  double_t *y_ptr = mesh["coordsets/surface_coords/values/y"].value();
  double_t *z_ptr = mesh["coordsets/surface_coords/values/z"].value();
  int *conn = mesh["topologies/rbc/elements/connectivity"].value();
  
  for(int i=0; i<nvals; ++i) {
    x_ptr[i] = x[i][0];
    y_ptr[i] = x[i][1];
    z_ptr[i] = x[i][2];

  }

  // Fill the connectivity array
  for(int i = 0, j =0; i < nanglelist; i++) {
      conn[j++] = anglelist[i][0];
      conn[j++] = anglelist[i][1];
      conn[j++] = anglelist[i][2];
  }

  // print out results

  mesh["fields/velocity_magnitude"]["association"] = "vertex";   // data per particle (vertex)
  mesh["fields/velocity_magnitude"]["topology"]    = "rbc";
  mesh["fields/velocity_magnitude"]["values"].set(conduit::DataType::float64(nlocal));
  double *velMag = mesh["fields/velocity_magnitude"]["values"].value();

  for(int i=0; i<nlocal; ++i)
  {
      double vx = v[i][0];
      double vy = v[i][1];
      double vz = v[i][2];
      velMag[i] = std::sqrt(vx*vx + vy*vy + vz*vz);
  }
  
  //std::cout << mesh.to_yaml() << std::endl;

  conduit::Node verify_info;
if(!conduit::blueprint::mesh::verify(mesh, verify_info))
{
    std::cout << "Verify failed!" << std::endl;
    verify_info.print();
    exit(0);
}
  //std::cout << mesh.to_yaml() << std::endl;
  mAscent.publish(mesh);

  conduit::Node actions;
  mAscent.execute(actions);
  /*vtkDoubleArray *velocityDoubleArray = vtkDoubleArray::New();
  vtkDoubleArray *vorticityDoubleArray = vtkDoubleArray::New();
  vtkDoubleArray *velocityNormDoubleArray = vtkDoubleArray::New();
  */

//XXXNew local values added with domainBox 2/23/22*****
  //int nlx = velocityArray.getNx(); 
  //int nly = velocityArray.getNy();
  //int nlz = velocityArray.getNz();
  //plint myrank = global::mpi().getRank();

/*****************************************************
  velocityDoubleArray->SetNumberOfComponents(3);
  velocityDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz)); 

  vorticityDoubleArray->SetNumberOfComponents(3);
  vorticityDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz));

   velocityNormDoubleArray->SetNumberOfComponents(1);
   velocityNormDoubleArray->SetNumberOfTuples((nlx) * (nly) * (nlz));

   //plint EW = envelopeWidth;

//XXX Need to convert this to zero copy: FUTURE WORK

        Array<double,3> vel(0.,0.,0.); //jifu 5/31/2022
        Array<double,3> vor(0.,0.,0.);
  for (int k=0; k<nlz; k++)
  {
    for (int j=0; j<nly; j++)
    {
     for (int i=0; i<nlx; i++)
      {
        vel = velocityArray.get(i,j,k); //jifu 5/31/2022 Array<double 3> vel ->vel
        vor = vorticityArray.get(i,j,k);
        double norm = velocityNormArray.get(i,j,k);
        int index = (j) * (nlx) + (i) + (k) * (nlx) * (nly);
        velocityDoubleArray->SetTuple3(index,vel[0],vel[1],vel[2]);
        vorticityDoubleArray->SetTuple3(index,vor[0],vor[1],vor[2]);
        velocityNormDoubleArray->SetTuple1(index,norm);
      }
    }
  }
 GlobalDataAdaptor->AddPalabosData(velocityDoubleArray, vorticityDoubleArray, velocityNormDoubleArray, nx, ny, nz, domainBox, envelopeWidth);
*/
}
void Bridge::Execute(long ntimestep)
{
  //conduit::Node actions;
  //mAscent.execute(actions);
}
void Bridge::Finalize() {
  mAscent.close();
}

