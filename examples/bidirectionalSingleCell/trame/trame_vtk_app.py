import asyncio
import numpy as np

import queue
from multiprocessing import Process, Queue
from multiprocessing.managers import BaseManager

from vtkmodules.vtkCommonCore import vtkPoints, vtkDoubleArray, vtkUnsignedCharArray
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkRenderingCore import vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkPolyDataMapper, vtkActor 
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleSwitch  # noqa
from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkFiltersCore import vtkTubeFilter, vtkPolyDataNormals
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter

from trame.app import get_server, asynchronous
from trame.widgets import vuetify, vtk as vtk_widgets
from trame.ui.vuetify import SinglePageLayout


class QueueManager(BaseManager):
    pass


def main():
    # create queues for Trame state and updates
    state_queue = Queue()
    update_queue = Queue()
    
    # start Trame app in new thread
    trame_thread = Process(target=runTrameServer, args=(state_queue, update_queue))
    trame_thread.daemon = True
    trame_thread.start()

    # create queues from Ascent data
    queue_data = Queue()
    queue_signal = Queue()

    # start Queue Manager in new thread
    queue_mgr_thread = Process(target=runQueueManager, args=(queue_data, queue_signal))
    queue_mgr_thread.daemon = True
    queue_mgr_thread.start()
   
    # wait for data coming from Ascent 
    data = None
    while data != '':
        print('waiting on data... ', end='')
        sim_data = queue_data.get()
        print(f'received!')
        
        state_queue.put(sim_data)
        updates = update_queue.get()

        queue_signal.put(updates)


def runTrameServer(state_queue, update_queue):
    # set up Trame application
    server = get_server(client_type="vue2")
    state = server.state
    ctrl = server.controller
    
    # create VTK render data
    vtk_data = {}
    vtk_data['renderer'] = vtkRenderer()
    vtk_data['renderer'].SetBackground(0.2, 0.3, 0.4)
    vtk_data['renderWindow'] = vtkRenderWindow()
    vtk_data['renderWindow'].AddRenderer(vtk_data['renderer'])

    vtk_data['renderWindowInteractor'] = vtkRenderWindowInteractor()
    vtk_data['renderWindowInteractor'].SetRenderWindow(vtk_data['renderWindow'])
    vtk_data['renderWindowInteractor'].GetInteractorStyle().SetCurrentStyleToTrackballCamera()

    num_vessel_points = 32
    z_step = (80.0 - 0.0) / (num_vessel_points - 1)
    vtk_data['vessel_points'] = vtkPoints()
    for i in range(num_vessel_points):
        z = 0.0 + (i * z_step)
        vtk_data['vessel_points'].InsertNextPoint([16.0, 16.0, z])
    vtk_data['vessel_center_line'] = vtkCellArray()
    vtk_data['vessel_center_line'].InsertNextCell(num_vessel_points)
    for i in range(num_vessel_points):
        vtk_data['vessel_center_line'].InsertCellPoint(i)
    vtk_data['vessel_polydata'] = vtkPolyData()
    vtk_data['vessel_polydata'].SetPoints(vtk_data['vessel_points'])
    vtk_data['vessel_polydata'].SetLines(vtk_data['vessel_center_line'])
    vtk_data['vessel_radius'] = vtkDoubleArray()
    vtk_data['vessel_radius'].SetName('TubeRadius')
    vtk_data['vessel_radius'].SetNumberOfTuples(num_vessel_points)
    vtk_data['vessel_colors'] = vtkUnsignedCharArray()
    vtk_data['vessel_colors'].SetName('TubeColors')
    vtk_data['vessel_colors'].SetNumberOfTuples(num_vessel_points)
    vtk_data['vessel_colors'].SetNumberOfComponents(3)
    for i in range(num_vessel_points):
        vtk_data['vessel_radius'].SetTuple1(i, 18.0)
        vtk_data['vessel_colors'].InsertTuple3(i, 255, 255, 255)
    vtk_data['vessel_polydata'].GetPointData().AddArray(vtk_data['vessel_radius'])
    vtk_data['vessel_polydata'].GetPointData().AddArray(vtk_data['vessel_colors'])
    vtk_data['vessel_polydata'].GetPointData().SetActiveScalars('TubeRadius')
    vtk_data['vessel_tube'] = vtkTubeFilter()
    vtk_data['vessel_tube'].SetInputData(vtk_data['vessel_polydata'])
    vtk_data['vessel_tube'].SetNumberOfSides(16)
    vtk_data['vessel_tube'].SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    vtk_data['vessel_mapper'] = vtkPolyDataMapper()
    vtk_data['vessel_actor'] = vtkActor()
    vtk_data['vessel_mapper'].SetInputConnection(vtk_data['vessel_tube'].GetOutputPort())
    vtk_data['vessel_mapper'].ScalarVisibilityOn()
    vtk_data['vessel_mapper'].SetScalarModeToUsePointFieldData()
    vtk_data['vessel_mapper'].SelectColorArray('TubeColors')
    vtk_data['vessel_actor'].SetMapper(vtk_data['vessel_mapper'])
    vtk_data['vessel_actor'].GetProperty().FrontfaceCullingOn()
    vtk_data['renderer'].AddActor(vtk_data['vessel_actor'])

    vertices = [[7.0, 7.0, 18.5], [25.0, 7.0, 36.1], [16.0, 25.0, 27.3]] # array of 3D points
    faces = [[0, 1, 2]] # array of triangles - each triangle has 3 vertex indices
    vtk_data['polyData'] = createTriangleVtkPolyData(vertices, faces)
    vtk_data['mapper'] = vtkPolyDataMapper()
    vtk_data['actor'] = vtkActor()
    vtk_data['mapper'].SetInputData(vtk_data['polyData'])
    vtk_data['actor'].SetMapper(vtk_data['mapper'])
    vtk_data['renderer'].AddActor(vtk_data['actor'])
    vtk_data['renderer'].ResetCamera()

    vtk_data['trame_view'] = None

    # start loop to check for state updates
    @ctrl.add('on_server_ready')
    def initVtkApp(**kwargs):
        print('SERVER READY')
        asynchronous.create_task(checkForStateUpdates(state, state_queue, update_queue, vtk_data))

    # callback for steering enabled change
    def uiStateEnableSteeringUpdate(enable_steering, **kwargs):
        if state.connected:
            state.allow_submit = enable_steering
        if not enable_steering:
            update_queue.put({})

    # callback for clicking submit button
    def submitSteeringOptions():
        steering_data = {}
        update_queue.put(steering_data)

    # register callbacks
    state.change('enable_steering')(uiStateEnableSteeringUpdate)

    # define webpage layout
    state.allow_submit = False
    with SinglePageLayout(server) as layout:
        layout.title.set_text('Ascent-Trame')
        with layout.toolbar:
            vuetify.VDivider(vertical=True, classes='mx-2')
            vuetify.VSwitch(
                label='Enable Steering',
                v_model=('enable_steering', True),
                hide_details=True,
                dense=True
            )
            vuetify.VSpacer()
            vuetify.VBtn(
                'Submit',
                color='primary',
                disabled=('!allow_submit',),
                click=submitSteeringOptions
            )
        with layout.content:
            with vuetify.VContainer(fluid=True, classes='pa-0 fill-height'):
                vtk_data['trame_view'] = vtk_widgets.VtkLocalView(vtk_data['renderWindow'], ref='view')
                vtk_data['trame_view'].reset_camera()
    # start Trame server
    server.start()


async def checkForStateUpdates(state, state_queue, update_queue, vtk_data):
    while True:
        try:
            state_data = state_queue.get(block=False)
           
            state.connected = True
            if state.enable_steering:
                state.allow_submit = True

            vtk_data['renderer'].RemoveActor(vtk_data['actor'])

            vtk_data['polyData'] = createTriangleVtkPolyData(state_data['vertices'], state_data['faces'])

            # maybe?
            #vtk_data['polyDataSmooth'] = vtkPolyDataNormals()
            #vtk_data['polyDataSmooth'].SetInputData(vtk_data['polyData'])

            vtk_data['mapper'] = vtkPolyDataMapper()
            vtk_data['actor'] = vtkActor()
            vtk_data['mapper'].SetInputData(vtk_data['polyData'])
            #vtk_data['mapper'].SetInputConnection(vtk_data['polyDataSmooth'].GetOutputPort())
            vtk_data['actor'].SetMapper(vtk_data['mapper'])
            vtk_data['actor'].GetProperty().SetColor(0.68, 0.05, 0.05)
            vtk_data['renderer'].AddActor(vtk_data['actor'])
           
            vtk_data['trame_view'].update()

            if not state.enable_steering:
                update_queue.put({})
        except queue.Empty:
            pass
        except Exception as e:
            print(f'ERROR: {type(e)}')
            print(e)
        await asyncio.sleep(0)


def runQueueManager(queue_data, queue_signal):
    # register queues with Queue Manager
    QueueManager.register('get_data_queue', callable=lambda:queue_data)
    QueueManager.register('get_signal_queue', callable=lambda:queue_signal)
    
    # create Queue Manager
    mgr = QueueManager(address=('127.0.0.1', 8000), authkey=b'ascent-trame')
    
    # start Queue Manager server
    server = mgr.get_server()
    server.serve_forever()


def createTriangleVtkPolyData(vertices, faces):
    points = vtkPoints()
    for v in vertices:
        points.InsertNextPoint(v)

    triangles = vtkCellArray()
    for f in faces:
        triangles.InsertNextCell(3)
        for i in f:
            triangles.InsertCellPoint(i)

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(triangles)

    return polydata

if __name__ == '__main__':
    main()

