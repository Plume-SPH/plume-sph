#import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
total_proc=16
maindirectory='/eng/home/zhixuanc/Desktop/rohit-zhixuanc/Documents/C_plusplus/Resutls/'
subdirectory='ccr/new1009/short1/'
Disp={}
Plot={}
for i in range(0, 32):
   if i<10:
      name='pvplot00'+str(i)+'.h5part'
   elif i<100:
      name='pvplot0'+str(i)+'.h5part'
   elif i<1000:
      name='pvplot'+str(i)+'.h5part'
   
   if i<10:
      plotname='pvplot00'+str(i)+'h5part'
   elif i<100:
      plotname='pvplot0'+str(i)+'h5part'
   elif i<1000:
      plotname='pvplot'+str(i)+'h5part'   
# create a new 'H5PartReader'
   Plot[plotname] = H5PartReader(FileName=maindirectory+subdirectory+name)
   Plot[plotname].Xarray = 'x'
   Plot[plotname].Yarray = 'y'
   Plot[plotname].Zarray = 'z'
   Plot[plotname].PointArrays = ['ID', 'Involved', 'Rho', 'Vx', 'Vy', 'Vz', 'bctp', 'engr', 'guest', 'mssfrc', 'phase', 'prss', 'x', 'y', 'z']

   # get animation scene
   animationScene1 = GetAnimationScene()

   # update animation scene based on data timesteps
   animationScene1.UpdateAnimationUsingDataTimeSteps()

   # show color bar/color legend
   #h5PartReader2Display.SetScalarBarVisibility(renderView1, True)

   # update animation scene based on data timesteps
   #animationScene1.UpdateAnimationUsingDataTimeSteps()

   # set active source
   SetActiveSource(Plot[plotname])

   # get active view
   renderView1 = GetActiveViewOrCreate('RenderView')
   # show data in view
   displayname=plotname+'Display'
   Disp[displayname] = Show(Plot[plotname], renderView1)

   # reset view to fit data
   renderView1.ResetCamera()

   # show color bar/color legend
   #Disp[displayname].SetScalarBarVisibility(renderView1, True)

   # set scalar coloring
   ColorBy(Disp[displayname], ('POINTS', 'mssfrc'))

   # rescale color and/or opacity maps used to include current data range
   Disp[displayname].RescaleTransferFunctionToDataRange(True)

   # show color bar/color legend
   Disp[displayname].SetScalarBarVisibility(renderView1, True)

   # get color transfer function/color map for 'mssfrc'
   mssfrcLUT = GetColorTransferFunction('mssfrc')

   # get opacity transfer function/opacity map for 'mssfrc'
   mssfrcPWF = GetOpacityTransferFunction('mssfrc')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
