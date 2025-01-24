import conduit
import ascent.mpi

import numpy as np
from multiprocessing.managers import BaseManager

from mpi4py import MPI


class QueueManager(BaseManager):
    pass


def main():
    # obtain a mpi4py mpi comm object
    comm = MPI.Comm.f2py(ascent_mpi_comm_id())

    # get task id and number of total tasks
    task_id = comm.Get_rank()
    num_tasks = comm.Get_size()

    # run Trame tasks
    interactive = np.array([False], bool)
    update_data = None
    if task_id == 0:
        update_data = executeMainTask(task_id, num_tasks, comm)
    else:
        executeDependentTask(task_id, num_tasks, comm)
    
    # broadcast updates to all ranks
    update_data = comm.bcast(update_data, root=0)

    # TODO: pass updates to Ascent callback
    #update_node = conduit.Node()
    #update_node['task_id'] = task_id
    #output = conduit.Node()
    #ascent.mpi.execute_callback('steeringCallback', update_node, output)


def executeMainTask(task_id, num_tasks, comm):
    interactive = np.array([False], bool)
    update_data = {}    

    # attempt to connect to Trame queue manager
    QueueManager.register('get_data_queue')
    QueueManager.register('get_signal_queue')
    mgr = QueueManager(address=('127.0.0.1', 8000), authkey=b'ascent-trame')
    try:
        mgr.connect()
        interactive[0] = True
    except:
        mgr = None

    # broadcast to all processes whether Trame is currently running
    comm.Bcast((interactive, 1, MPI.BOOL), root=0)

    if interactive[0]:
        # get access to Trame's queues
        queue_data = mgr.get_data_queue()
        queue_signal = mgr.get_signal_queue()

        # get published blueprint data
        mesh_data = ascent_data().child(0)

        # TODO: repartition data -> gather on main process (0)
        #result = repartitionMeshData(task_id, num_tasks, comm)

        # extract red blood cell mesh polydata
        x_coords = mesh_data['coordsets/particle_coords/values/x']
        y_coords = mesh_data['coordsets/particle_coords/values/y']
        z_coords = mesh_data['coordsets/particle_coords/values/z']
        vertices = np.column_stack([x_coords, y_coords, z_coords])

        faces = np.reshape(mesh_data['topologies/particle_topo/elements/connectivity'], (-1, 3))

        # send rbc data to Trame
        queue_data.put({'vertices': vertices, 'faces': faces})
        
        # get steering updates from Trame
        update_data = queue_signal.get()

    return update_data


def executeDependentTask(task_id, num_tasks, comm):
    interactive = np.array([False], bool)

    # receive whether session is interactive or not from main task
    comm.Bcast((interactive, 1, MPI.BOOL), root=0)

    if interactive[0]:
        # TODO: repartition data -> gather on main process (0)
        #repartitionMeshData(task_id, num_tasks, comm)
        pass


if __name__ == "__main__":
    main()

"""
import sys

from hmac import new
import os
import zipfile

from multiprocessing.managers import BaseManager

import ascent.mpi
import conduit
import numpy as np
import yaml

from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter

from mpi4py import MPI


class QueueManager(BaseManager):
    pass


def main():
    # obtain a mpi4py mpi comm object
    comm = MPI.Comm.f2py(ascent_mpi_comm_id())

    # get task id and number of total tasks
    task_id = comm.Get_rank()
    num_tasks = comm.Get_size()

    ###
    # TEST
    ###
    # get published blueprint data
    mesh_data = ascent_data().child(0)    

    x_coords = mesh_data['coordsets/particle_coords/values/x']
    y_coords = mesh_data['coordsets/particle_coords/values/y']
    z_coords = mesh_data['coordsets/particle_coords/values/z']
    vertices = np.column_stack([x_coords, y_coords, z_coords])

    faces = np.reshape(mesh_data['topologies/particle_topo/elements/connectivity'], (-1, 3))

    poly_data = createTriangleVtkPolyData(vertices, faces)

    writer = vtkXMLPolyDataWriter()
    writer.SetFileName('ascent_rbc.vtp')
    writer.SetInputData(poly_data)
    writer.Write()
    ###
    # END: TEST
    ###

    # run Trame tasks
    interactive = np.array([False], bool)
    update_data = None
    if task_id == 0:
        update_data = executeMainTask(comm)
    else:
        executeDependentTask(comm)

    # broadcast updates to all ranks
    update_data = comm.bcast(update_data, root=0)


# Realized that the client is pretty much always going to run where the server is, making this not necessary
def create_zip_file(source_dir, zip_name="cinema_data.zip"):
    with zipfile.ZipFile(zip_name, "w") as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                zipf.write(file_path, arcname=os.path.relpath(file_path, source_dir))
    return zip_name


def executeMainTask(comm):
    interactive = np.array([False], bool)
    update_data = {}

    # attempt to connect to Trame queue manager
    QueueManager.register("get_data_queue")
    QueueManager.register("get_signal_queue")
    mgr = QueueManager(address=("127.0.0.1", 8000), authkey=b"ascent-trame")
    try:
        mgr.connect()
        interactive[0] = True
    except Exception as e:
        print(f"Failed to connect to Trame queue manager: {e}")
        mgr = None

    # broadcast to all processes whether Trame is currently running
    comm.Bcast((interactive, 1, MPI.BOOL), root=0)

    if interactive[0]:
        queue_data = mgr.get_data_queue()
        queue_signal = mgr.get_signal_queue()

        zip_file = create_zip_file("cinema_databases")
        if not os.path.exists(zip_file):
            print(f"Zip file creation failed: {zip_file} not found.")
            return
        print(f"Zip file created at: {zip_file}")

        try:
            with open(zip_file, "rb") as f:
                zip_data = f.read()
                queue_data.put({"zip_content": zip_data})
                os.remove(zip_file)
            print("Zip file sent to server.")
        except Exception as e:
            print(f"Failed to send zip file: {e}")

        update_data = queue_signal.get()

        update_node = conduit.Node()
        output_node = conduit.Node()

        if "new_cell" in update_data:
            ascent.mpi.execute_callback("addRandomCell", update_node, output_node)

    return update_data


def executeDependentTask(comm):
    interactive = np.array([False], bool)

    # receive whether session is interactive or not from main task
    comm.Bcast((interactive, 1, MPI.BOOL), root=0)

    if interactive[0]:
        pass


def createTriangleVtkPolyData(vertices, faces):
    points = vtkPoints()
    for v in vertices:
        points.InsertNextPoint(v)

    triangles = vtkCellArray()
    for f in faces:
        triangles.InsertNextCell(3)
        for i in f:
            triangles.InsertCellPoint(i)

    poly_data = vtkPolyData()
    poly_data.SetPoints(points)
    poly_data.SetPolys(triangles)

    return poly_data


if __name__ == "__main__":
    main()

"""

