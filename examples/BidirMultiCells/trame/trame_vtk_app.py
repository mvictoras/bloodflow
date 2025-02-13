import asyncio
import numpy as np

import queue
from multiprocessing import Process, Queue
from multiprocessing.managers import BaseManager

from PIL import Image

from vtkmodules.vtkCommonCore import vtkPoints, vtkDoubleArray, vtkUnsignedCharArray, VTK_UNSIGNED_CHAR
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData, vtkDataObject, vtkImageData
from vtkmodules.vtkRenderingCore import (vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
                                         vtkTexture, vtkPolyDataMapper, vtkActor) 
from vtkmodules.vtkFiltersSources import vtkPlaneSource
from vtkmodules.vtkFiltersCore import vtkTubeFilter, vtkCleanPolyData, vtkPolyDataNormals
from vtkmodules.vtkIOImage import vtkPNGReader
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleSwitch  # noqa

from vtkmodules.util import numpy_support

from trame.app import get_server, asynchronous
#from trame.widgets import vuetify, vtk as vtk_widgets
from trame.widgets import vuetify, vtklocal as vtk_widgets
from trame.ui.vuetify import SinglePageLayout

DIMS = {'x': 90.0, 'y': 90.0, 'z': 150.0}

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
    z_step = DIMS['z'] / (num_vessel_points - 1)
    vtk_data['vessel_points'] = vtkPoints()
    for i in range(num_vessel_points):
        z = i * z_step
        vtk_data['vessel_points'].InsertNextPoint([0.5 * DIMS['x'], 0.5 * DIMS['y'], z])
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
        vtk_data['vessel_radius'].SetTuple1(i, 0.5 * DIMS['x'])
        vtk_data['vessel_colors'].InsertTuple3(i, 255, 255, 255)
    vtk_data['vessel_polydata'].GetPointData().AddArray(vtk_data['vessel_radius'])
    vtk_data['vessel_polydata'].GetPointData().AddArray(vtk_data['vessel_colors'])
    vtk_data['vessel_polydata'].GetPointData().SetActiveScalars('TubeRadius')
    vtk_data['vessel_tube'] = vtkTubeFilter()
    vtk_data['vessel_tube'].SetInputData(vtk_data['vessel_polydata'])
    vtk_data['vessel_tube'].SetNumberOfSides(24)
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

    vtk_data['fluid_yz_plane'] = vtkPlaneSource()
    vtk_data['fluid_yz_plane'].SetOrigin(0.5 * DIMS['x'], 0.0, DIMS['z'])
    vtk_data['fluid_yz_plane'].SetPoint1(0.5 * DIMS['x'], 0.0, 0.0)
    vtk_data['fluid_yz_plane'].SetPoint2(0.5 * DIMS['x'], DIMS['y'], DIMS['z'])
    vtk_data['fluid_yz_mapper'] = vtkPolyDataMapper()
    vtk_data['fluid_yz_actor'] = vtkActor()
    vtk_data['fluid_yz_mapper'].SetInputConnection(vtk_data['fluid_yz_plane'].GetOutputPort())
    vtk_data['fluid_yz_actor'].SetMapper(vtk_data['fluid_yz_mapper'])
    vtk_data['renderer'].AddActor(vtk_data['fluid_yz_actor'])

    vtk_data['fluid_xz_plane'] = vtkPlaneSource()
    vtk_data['fluid_xz_plane'].SetOrigin(0.0, 0.5 * DIMS['y'], 0.0)
    vtk_data['fluid_xz_plane'].SetPoint1(0.0, 0.5 * DIMS['y'], DIMS['z'])
    vtk_data['fluid_xz_plane'].SetPoint2(DIMS['x'], 0.5 * DIMS['y'], 0.0)
    vtk_data['fluid_xz_mapper'] = vtkPolyDataMapper()
    vtk_data['fluid_xz_actor'] = vtkActor()
    vtk_data['fluid_xz_mapper'].SetInputConnection(vtk_data['fluid_xz_plane'].GetOutputPort())
    vtk_data['fluid_xz_actor'].SetMapper(vtk_data['fluid_xz_mapper'])
    vtk_data['renderer'].AddActor(vtk_data['fluid_xz_actor'])

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

    # callback for slice opacity slider
    def uiStateUpdateSliceOpacity(slice_opacity, **kwargs):
        vtk_data['fluid_yz_actor'].GetProperty().SetOpacity(slice_opacity)
        vtk_data['fluid_xz_actor'].GetProperty().SetOpacity(slice_opacity)
        vtk_data['trame_view'].update()

    # callback for clicking submit button
    def submitSteeringOptions():
        steering_data = {}
        update_queue.put(steering_data)

    # register callbacks
    state.change('enable_steering')(uiStateEnableSteeringUpdate)
    state.change('slice_opacity')(uiStateUpdateSliceOpacity)

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
            vuetify.VSlider(
                label='Slice Opacity',
                v_model=('slice_opacity', 0.5),
                min=0.0,
                max=1.0,
                step=0.05,
                hide_details=True,
                dense=True
            )
            vuetify.VCol(
                '{{slice_opacity.toFixed(2)}}'
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
                #vtk_data['trame_view'] = vtk_widgets.VtkLocalView(vtk_data['renderWindow'], ref='view')
                vtk_data['trame_view'] = vtk_widgets.LocalView(vtk_data['renderWindow'], ref='view')
                vtk_data['trame_view'].reset_camera()

    # start Trame server
    server.start()


async def checkForStateUpdates(state, state_queue, update_queue, vtk_data):
    while True:
        try:
            state_data = state_queue.get(block=False)
           
            state.connected = True
            
            print(f'connected: steering={state.enable_steering}')
            if state.enable_steering:
                state.allow_submit = True

            vtk_data['renderer'].RemoveActor(vtk_data['actor'])
            vtk_data['renderer'].RemoveActor(vtk_data['fluid_yz_actor'])
            vtk_data['renderer'].RemoveActor(vtk_data['fluid_xz_actor'])

            vtk_data['polyData'] = createTriangleVtkPolyData(state_data['vertices'], state_data['faces'])
            #vtk_data['polyDataClean'] = vtkCleanPolyData()
            #vtk_data['polyDataClean'].SetInputData(vtk_data['polyData'])

            # maybe?
            #vtk_data['polyDataSmooth'] = vtkPolyDataNormals()
            #vtk_data['polyDataSmooth'].SetInputConnection(vtk_data['polyDataClean'].GetOutputPort())

            vtk_data['mapper'] = vtkPolyDataMapper()
            vtk_data['actor'] = vtkActor()
            vtk_data['mapper'].SetInputData(vtk_data['polyData'])
            #vtk_data['mapper'].SetInputConnection(vtk_data['polyDataClean'].GetOutputPort())
            #vtk_data['mapper'].SetInputConnection(vtk_data['polyDataSmooth'].GetOutputPort())
            vtk_data['actor'].SetMapper(vtk_data['mapper'])
            vtk_data['actor'].GetProperty().SetColor(0.68, 0.05, 0.05)
            vtk_data['renderer'].AddActor(vtk_data['actor'])

            vtk_data['image_reader_yz'] = vtkPNGReader()
            vtk_data['image_reader_yz'].SetFileName('../fluid_vel_mag_yz.png')
            vtk_data['image_reader_yz'].Update()

            vtk_data['fluid_yz_tex'] = vtkTexture()
            vtk_data['fluid_yz_tex'].SetInputConnection(vtk_data['image_reader_yz'].GetOutputPort())
            vtk_data['fluid_yz_tex'].InterpolateOn()

            #fluid_yz_tex = createVtkTextureFromImage('../fluid_vel_mag_yz.png')
            vtk_data['fluid_yz_actor'] = vtkActor()
            vtk_data['fluid_yz_actor'].SetMapper(vtk_data['fluid_yz_mapper'])
            vtk_data['fluid_yz_actor'].SetTexture(vtk_data['fluid_yz_tex'])
            vtk_data['fluid_yz_actor'].GetProperty().SetOpacity(state.slice_opacity)
            vtk_data['renderer'].AddActor(vtk_data['fluid_yz_actor'])

            vtk_data['image_reader_xz'] = vtkPNGReader()
            vtk_data['image_reader_xz'].SetFileName('../fluid_vel_mag_xz.png')
            vtk_data['image_reader_xz'].Update()

            vtk_data['fluid_xz_tex'] = vtkTexture()
            vtk_data['fluid_xz_tex'].SetInputConnection(vtk_data['image_reader_xz'].GetOutputPort())
            vtk_data['fluid_xz_tex'].InterpolateOn()

            #fluid_xz_tex = createVtkTextureFromImage('../fluid_vel_mag_xz.png')
            vtk_data['fluid_xz_actor'] = vtkActor()
            vtk_data['fluid_xz_actor'].SetMapper(vtk_data['fluid_xz_mapper'])
            vtk_data['fluid_xz_actor'].SetTexture(vtk_data['fluid_xz_tex'])
            vtk_data['fluid_xz_actor'].GetProperty().SetOpacity(state.slice_opacity)
            vtk_data['renderer'].AddActor(vtk_data['fluid_xz_actor'])

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

    max_sqr_dist = 75.0 * 75.0
    triangles = vtkCellArray()
    for f in faces:
        if f[0] >= len(vertices) or f[1] >= len(vertices) or f[2] >= len(vertices):
            #print(f'bad face: {f}')
            pass
        else:
            d01 = sqrDistance3d(vertices[f[0]], vertices[f[1]])
            d02 = sqrDistance3d(vertices[f[0]], vertices[f[2]])
            d12 = sqrDistance3d(vertices[f[1]], vertices[f[2]])
            if d01 < max_sqr_dist and d02 < max_sqr_dist and d12 < max_sqr_dist:
                triangles.InsertNextCell(3)
                triangles.InsertCellPoint(f[0])
                triangles.InsertCellPoint(f[1])
                triangles.InsertCellPoint(f[2])

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(triangles)

    return polydata


def createVtkTextureFromImage(filename):
    print(f'loading PNG: {filename}')
    image_reader = vtkPNGReader()
    image_reader.SetFileName(filename)
    image_reader.Update()   
 
    texture = vtkTexture()
    texture.SetInputConnection(image_reader.GetOutputPort())
    texture.InterpolateOn()
    texture.Update()

    return texture

def sqrDistance3d(p0, p1):
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    dz = p1[2] - p0[2]
    return (dx * dx) + (dy * dy) + (dz * dz)

if __name__ == '__main__':
    main()

