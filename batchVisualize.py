def batchVis(c1File,particleFile,step,saveAs):
    """Renders a bijel top down view in paraview and saves a screenshot."""
    import paraview.simple as pv
    # visualize a vtk file 
    c1 = pv.LegacyVTKReader(FileNames=c1File)
    p = pv.LegacyVTKReader(FileNames=particleFile)
    renderView1 = pv.GetActiveViewOrCreate('RenderView')
    renderView1.ViewSize = [1298, 860]
    renderView1.Background = [1.0, 1.0, 1.0]
    renderView1.InteractionMode = '2D'
    pDisplay = pv.Show(p, renderView1)
    c1Display = pv.Show(c1, renderView1)

    # create particle glyphs
    glyph = pv.Glyph(Input=p,GlyphType="Sphere")
    glyph.ScaleFactor = 1.0
    glyph.GlyphMode = 'All Points'
    glyph.GlyphType.Radius = 1.0
    glyph.GlyphType.ThetaResolution = 20
    glyph.GlyphType.PhiResolution = 20
    glyph.Scalars = ['POINTS','radius']
    glyph.Vectors = ['POINTS','None']
    glyph.ScaleMode = 'scalar'

    # show data in view
    glyphDisplay = pv.Show(glyph, renderView1)
    pv.ColorBy(glyphDisplay, None)
    pv.SetActiveSource(c1)
    pv.ColorBy(c1Display, ('POINTS', 'c1'))
    c1Display.RescaleTransferFunctionToDataRange(True)
    c1Display.SetRepresentationType('Volume')

    # make box outline
    # box = pv.Box()
    # box.XLength = 128.0
    # box.YLength = 128.0
    # box.ZLength = 64.0
    # box.Center = [64.0, 64.0, 32.0]
    # boxDisplay = pv.Show(box, renderView1)
    # boxDisplay.SetRepresentationType('Outline')
    # boxDisplay.AmbientColor = [0.0, 0.0, 0.0]

    # set coloring of c1
    c1LUT = pv.GetColorTransferFunction('c1')
    c1LUT.RGBPoints = [0.006000000052154064, 0.231373, 0.298039, 0.752941, 0.5120000033639371, 0.865003, 0.865003, 0.865003, 1.0180000066757202, 0.705882, 0.0156863, 0.14902]
    c1LUT.ColorSpace = 'Diverging'
    c1LUT.BelowRangeColor = [0.0, 0.0, 0.0]
    c1LUT.AboveRangeColor = [1.0, 1.0, 1.0]
    c1LUT.NanColor = [1.0, 1.0, 0.0]
    c1LUT.Discretize = 1
    c1LUT.NumberOfTableValues = 256
    c1LUT.ScalarRangeInitialized = 1.0
    c1LUT.AllowDuplicateScalars = 1

    c1PWF = pv.GetOpacityTransferFunction('c1')
    c1PWF.Points = [0.0, 0.05, 0.5, 0.0, 0.3, 0.05, 0.5, 0.0, 0.4, 0.5, 0.5, 0.0, 0.6, 0.5, 0.5, 0.0, 0.7, 0.05, 0.5, 0.0, 1., 0.05, 0.5, 0.0]

    # annotate time step in rendering
    # text = pv.Text
    # text.Text = 'Step '+str(step)
    # textDisplay = pv.Show(text,renderView1)
    # textDisplay.Color = [0.0, 0.0, 0.0]
    # textDisplay.WindowLocation = 'UpperCenter'

    # reset view to fit data
    renderView1.ResetCamera()
    # pv.Render()

    # save screen shot
    viewLayout1 = pv.GetLayout()
    print(saveAs)
    pv.SaveScreenshot(saveAs, layout=viewLayout1, magnification=1, quality=100)

    # clean up
    # pv.Delete(box)
    pv.Delete(glyph)
    pv.Delete(p)
    pv.Delete(c1)
    del c1
    del p
    del glyph
    # del box


def getLargestStepVtk(tagName,path):
    """Returns vtk file with largest step associated with tagName."""
    files = os.listdir(path)
    for f in sorted(files):
        tagCheck = f.split('.')[0].split('_')[0]
        if tagName != tagCheck:
            files.remove(f)
    # sort files according to step (i.e. assumes file names are formatted as tagname_step.vtk)
    files.sort(key=lambda x: int(x.split('.')[0].split('_')[1]))
    # return the file with largest step value
    selectedVtk = files[len(files)-1]
    return selectedVtk



# ------------------------------------------------------
# main program
# ------------------------------------------------------

import os
import errno

# studyDir = '/home/joseph/Desktop/bijel_studies/film_geometry_e-field_tests'
# studyDir = '/home/joseph/Desktop/bijel_studies/film_geometry_e-field_tests1'
# studyDir = '/home/joseph/Desktop/bijel_studies/radius_e-field_study_run3'
# studyDir = '/home/joseph/Desktop/bijel_studies/r5_thickness_study'
# studyDir = '/home/joseph/Desktop/bijel_studies/radius_e-field_study/simData/radius_e-field_study_extremes'
# studyDir = '/home/joseph/Desktop/bijel_studies/radius_repulsive_capStr_study/simData/radius_repulsive_capStr_study'
# studyDir = '/home/joseph/Desktop/bijel_studies/radius_e-field_study/simData/r3_e-field_study'
# studyDir = '/home/joseph/Desktop/bijel_studies/radius_e-field_study/simData/radius_e-field_no_eCH_study_run1'
# studyDir = '/media/joseph/external_drive_1/simulation_data/electric_field_studies/radius_e-field_study/simData/radius_e-field_study_run4'
# studyDir = '/media/joseph/external_drive_1/simulation_data/electric_field_studies/thickness_radius_study/simData/r3_thickness_study_run1'
# studyDir = '/media/joseph/external_drive_1/simulation_data/electric_field_studies/thickness_radius_study/simData/r5_thickness_study_run1'
studyDir = '/media/joseph/external_drive_1/simulation_data/electric_field_studies/blend_ratio_and_particle_volume_fraction_studies/simData/highKz_blendRatio_pvf_C_study_run1'

# create the images directory for saving visualizations in 
saveDir = studyDir+'/images'
try:
    os.makedirs(saveDir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
    else:
        os.system('rm -rf '+saveDir+'/*')

# get list of simulation directories from parStudyInfo.txt file (generated by parstubuilder)
simList = []
with open(studyDir+'/parStudyInfo.txt') as fin:
    for line in fin:
        if 'Unique parameter set' in line:
            break
    for line in fin:
        simList.append(line.rstrip())
fin.close()

# find the data files to visualize for each simulation
c1_visList = []
particle_visList = []
steps = []
for pth in sorted(simList):
    dataPath = studyDir+'/'+pth+'/vtkoutput'
    # getLargestStepVtk function used below can be replaced with an alternative
    # get concentration vtk file
    c1File = getLargestStepVtk('c1',dataPath)       
    # get particle vtk file
    parFile = getLargestStepVtk('particles',dataPath)
    # get step from the vtk file name
    step = parFile.split('.')[0].split('_')[1]
    # save this info in some arrays for later processing
    steps.append(step)
    c1_visList.append(c1File)
    particle_visList.append(parFile)

# visualize the data in paraview and save it. batchVis can be replaced with another
# paraview script.
for i, f in enumerate(c1_visList):
    c1Path = studyDir+'/'+simList[i]+'/vtkoutput/'+f
    parPath = studyDir+'/'+simList[i]+'/vtkoutput/'+particle_visList[i]
    savePath = saveDir+'/'+simList[i]+'.png'
    batchVis(c1Path,parPath,steps[i],savePath)
