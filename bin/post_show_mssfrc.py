###import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

total_proc= 16
#Render={}
#Disp={}
#get active source.
#ActiveReader= GetActiveSource()

# get active view
#renderView1 = GetActiveViewOrCreate('RenderView')

for i in range(total_proc):
   name='H5PartReader'+str(i+1)
      
   # find source
   Render= FindSource(name)

   # set active source
   SetActiveSource(Render)

   # get active view
   renderView1 = GetActiveViewOrCreate('RenderView')

   # get display properties
   Disp = GetDisplayProperties(Render, view=renderView1)

   # set scalar coloring
   ColorBy(Disp, ('POINTS', 'mssfrc'))

   # rescale color and/or opacity maps used to include current data range
   Disp.RescaleTransferFunctionToDataRange(True)

   # show color bar/color legend
   Disp.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Vx'
mssfrcLUT = GetColorTransferFunction('mssfrc')

# get opacity transfer function/opacity map for 'mssfrc'
mssfrcPWF = GetOpacityTransferFunction('mssfrc')
