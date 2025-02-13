#include "AscentBridge.h"
#include <iostream>

using namespace std;
using namespace plb; 

std::unique_ptr<AscentBridge> AscentBridge::singleton_;

AscentBridge::AscentBridge() {
}

AscentBridge::~AscentBridge() {
}

// Define the static method to get the singleton instance
AscentBridge& AscentBridge::getInstance() {
  static std::once_flag flag;
  std::call_once(flag, [](){
      singleton_.reset(new AscentBridge());
      });
  return *singleton_;
}

void AscentBridge::Initialize(MPI_Comm world) {
  conduit::Node ascent_opts;
  ascent_opts["mpi_comm"] = MPI_Comm_c2f(world);
  mAscent.open(ascent_opts);
}

void AscentBridge::Publish  (double **x, double **v, long ntimestep, int nghost, 
                      int nlocal, int **anglelist, int nanglelist, 
                      TensorField3D<double, 3> velocityArray,
                      TensorField3D<double, 3> vorticityArray,
                      ScalarField3D<double> velocityNormArray,
                      int nx, int ny, int nz, Box3D domainBox, plint envelopeWidth) {

  

  conduit::Node mesh;
  const int num_tris = nanglelist;
  const int conn_len = num_tris * 3;
  vector<double> velocity_x, velocity_y, velocity_z;
  vector<double> vorticity_x, vorticity_y, vorticity_z;
  int nlx = velocityArray.getNx() - 2*envelopeWidth; //jifu: 2/12/2025 
  int nly = velocityArray.getNy() - 2*envelopeWidth;
  int nlz = velocityArray.getNz() - 2*envelopeWidth;

  mesh["state/cycle"] = ntimestep;
  mesh["coordsets/fluid_coords/type"] = "uniform";
  mesh["coordsets/fluid_coords/dims/i"] = nlx+1; 
  mesh["coordsets/fluid_coords/dims/j"] = nly+1;
  mesh["coordsets/fluid_coords/dims/k"] = nlz+1;
  mesh["topologies/fluid_topo/type"] = "uniform";
  mesh["topologies/fluid_topo/coordset"] = "fluid_coords";
  mesh["coordsets/fluid_coords/origin/x"] = domainBox.x0; //+ envelopeWidth;
  mesh["coordsets/fluid_coords/origin/y"] = domainBox.y0;// + envelopeWidth;
  mesh["coordsets/fluid_coords/origin/z"] = domainBox.z0;// + envelopeWidth;
  mesh["coordsets/fluid_coords/spacing/dx"] = 1;
  mesh["coordsets/fluid_coords/spacing/dy"] = 1;
  mesh["coordsets/fluid_coords/spacing/dz"] = 1;

  Array<double,3> vel(0.,0.,0.); //jifu 5/31/2022
  Array<double,3> vor(0.,0.,0.);
  for (int k=0; k<nlz; k++)
  {
    for (int j= 0; j<nly ; j++)
    {
     for (int i=0 ; i<nlx; i++)
      {
        vel = velocityArray.get(i+envelopeWidth,j+envelopeWidth,k+envelopeWidth); //jifu 5/31/2022 Array<double 3> vel ->vel
        vor = vorticityArray.get(i+envelopeWidth,j+envelopeWidth,k+envelopeWidth);
        double norm = velocityNormArray.get(i+envelopeWidth,j+envelopeWidth,k+envelopeWidth);
        int index = (j) * (nlx) + (i) + (k) * (nlx) * (nly);
        velocity_x.push_back(vel[0]);
        velocity_y.push_back(vel[1]);
        velocity_z.push_back(vel[2]);
        vorticity_x.push_back(vor[0]);
        vorticity_y.push_back(vor[1]);
        vorticity_z.push_back(vor[2]);
        //velocityNormDoubleArray->SetTuple1(index,norm);
      }
    }
  }

  // Now add velocity as separate arrays
  mesh["fields/fluid_velocity/association"] = "element";//"vertex";
  mesh["fields/fluid_velocity/topology"]    = "fluid_topo";
  mesh["fields/fluid_velocity/values/x"].set(velocity_x);
  mesh["fields/fluid_velocity/values/y"].set(velocity_y);
  mesh["fields/fluid_velocity/values/z"].set(velocity_z);

  mesh["fields/fluid_vorticity/association"] = "element";//"vertex";
  mesh["fields/fluid_vorticity/topology"]    = "fluid_topo";
  mesh["fields/fluid_vorticity/values/x"].set(vorticity_x);
  mesh["fields/fluid_vorticity/values/y"].set(vorticity_y);
  mesh["fields/fluid_vorticity/values/z"].set(vorticity_z);

  int nvals = nlocal;//+nghost;
  mesh["coordsets/particle_coords/type"] = "explicit";
  mesh["coordsets/particle_coords/values/x"].set(conduit::DataType::float64(nvals));
  mesh["coordsets/particle_coords/values/y"].set(conduit::DataType::float64(nvals));
  mesh["coordsets/particle_coords/values/z"].set(conduit::DataType::float64(nvals));

  mesh["topologies/particle_topo/type"]           = "unstructured";
  mesh["topologies/particle_topo/coordset"]       = "particle_coords";
  mesh["topologies/particle_topo/elements/shape"] = "tri";
  mesh["topologies/particle_topo/elements/connectivity"].set(conduit::DataType::int32(nanglelist * 3));

  // Get pointers to fill in
  double_t *x_ptr = mesh["coordsets/particle_coords/values/x"].value();
  double_t *y_ptr = mesh["coordsets/particle_coords/values/y"].value();
  double_t *z_ptr = mesh["coordsets/particle_coords/values/z"].value();
  int *conn = mesh["topologies/particle_topo/elements/connectivity"].value();
  
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

  mesh["fields/particle_velocity_magnitude"]["association"] = "vertex";   // data per particle (vertex)
  mesh["fields/particle_velocity_magnitude"]["topology"]    = "particle_topo";
  mesh["fields/particle_velocity_magnitude"]["values"].set(conduit::DataType::float64(nlocal));
  double *velMag = mesh["fields/particle_velocity_magnitude"]["values"].value();

  for(int i=0; i<nlocal; ++i)
  {
      double vx = v[i][0];
      double vy = v[i][1];
      double vz = v[i][2];
      velMag[i] = std::sqrt(vx*vx + vy*vy + vz*vz);
  }
  
  conduit::Node verify_info;
  if(!conduit::blueprint::mesh::verify(mesh, verify_info))
  {
      std::cerr << "Mesh verification failed!" << std::endl;
      verify_info.print();
      exit(EXIT_FAILURE);
  }
  
  mAscent.publish(mesh);

  conduit::Node actions;
  mAscent.execute(actions);
  
}

void AscentBridge::Finalize() {
  mAscent.close();
}

