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

    # pass updates to Ascent callback
    update_node = conduit.Node()
    update_node['task_id'] = task_id
    for key, item in update_data.items():
        update_node[key] = item
    output = conduit.Node()
    ascent.mpi.execute_callback('updateVesselStenosis', update_node, output)


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
        result = gatherRbcMeshData(task_id, num_tasks, comm)

        # extract red blood cell mesh polydata
        #x_coords = mesh_data['coordsets/particle_coords/values/x']
        #y_coords = mesh_data['coordsets/particle_coords/values/y']
        #z_coords = mesh_data['coordsets/particle_coords/values/z']
        #vertices = np.column_stack([x_coords, y_coords, z_coords])

        #faces = np.reshape(mesh_data['topologies/particle_topo/elements/connectivity'], (-1, 3))

        # send rbc data to Trame
        #queue_data.put({'vertices': vertices, 'faces': faces})
        queue_data.put(result)        

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
        gatherRbcMeshData(task_id, num_tasks, comm)
        pass


def gatherRbcMeshData(task_id, num_tasks, comm):
    result = None

    # get published blueprint data
    mesh_data = ascent_data().child(0)

    # extract red blood cell mesh polydata
    x_coords = mesh_data['coordsets/particle_coords/values/x']
    y_coords = mesh_data['coordsets/particle_coords/values/y']
    z_coords = mesh_data['coordsets/particle_coords/values/z']
    vertices = np.column_stack([x_coords, y_coords, z_coords])

    faces = np.reshape(mesh_data['topologies/particle_topo/elements/connectivity'], (-1, 3))

    all_vertices = comm.gather(vertices, root=0)
    all_faces = comm.gather(faces, root=0)
    
    if task_id == 0:
        vertex_counts = np.array([x.shape[0] for x in all_vertices], dtype=np.uint32)
        vertex_offsets = np.cumsum(vertex_counts, dtype=np.uint32)
        for i in range(len(all_faces)):
            # filter out all faces whose verts are out of range
            mask = np.apply_along_axis(np.any, 1, all_faces[i] >= vertex_counts[i])
            all_faces[i] = all_faces[i][~mask]   
            # add offset into global array
            if i >= 1:
                np.add(all_faces[i], vertex_offsets[i - 1], out=all_faces[i])
        
        all_vertices = np.concatenate(all_vertices)
        all_faces = np.concatenate(all_faces)

        result = {'vertices': all_vertices, 'faces': all_faces}
 
        f = open('sim_data.txt', 'w', encoding='utf-8')
        f.write(f'Vertex Counts\n{vertex_counts}\n')
        f.write(f'Vertex Offsets\n{vertex_offsets}\n')
        #f.write(f'Face Counts\n{face_count}\n')
        f.write(f'Vertex Data\n{all_vertices}\n')
        f.write(f'Face Data\n{all_faces}\n')
        f.close()

    return result


if __name__ == "__main__":
    main()

