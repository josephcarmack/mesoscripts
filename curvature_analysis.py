
def calc_curvature(filename,arr_name,gauss_filter,mean_filter,sample_rate):
    import paraview.simple as ps
    import numpy as np
    import paraview as pv
    from vtk.util.numpy_support import vtk_to_numpy

    # have paraview open the vtk data file
    reader = ps.OpenDataFile(filename)
    reader.UpdatePipeline()
    sys_data = pv.servermanager.Fetch(reader)
    nx,ny,nz = sys_data.GetDimensions()
    dx,dy,dz = sys_data.GetSpacing()

    # downsample the data (makes for a smoother contour surface)
    ds = ps.ExtractSubset(reader)
    ds.SampleRateI = sample_rate 
    ds.SampleRateJ = sample_rate
    ds.SampleRateK = sample_rate
    ds.VOI[1]=nx-1
    ds.VOI[3]=ny-1
    ds.VOI[5]=nz-1
    ds.IncludeBoundary = 1
    ds.UpdatePipeline()

    # have paraview apply a contour surface at a concentration value of cont_val
    contour = ps.Contour(ds)
    cont_val = 0.0						 # this might change depending on order parameter
    contour.ContourBy = ['POINTS',arr_name] # Viz is the name of the vtk array
    contour.Isosurfaces = [cont_val]	 
    contour.SetPropertyWithName('ComputeNormals',0)
    contour.UpdatePipeline()

    # have paraview calculate the curvature (default is Gaussian)
    curvature = ps.Curvature(contour)

    # filter the curvatures
    threshold = ps.Threshold(curvature)
    threshold.Scalars = ['POINTS','Gauss_Curvature']
    threshold.ThresholdRange = [-gauss_filter,gauss_filter]
    threshold.AllScalars = 1
    threshold.UpdatePipeline()
    gauss_data = pv.servermanager.Fetch(threshold)
    size_gauss = gauss_data.GetPointData().GetArray(0).GetSize()

    # integrate the surface area and curvature
    summed_curv = ps.IntegrateVariables(threshold)
    summed_curv_data = pv.servermanager.Fetch(summed_curv)
    g_surf_area = summed_curv_data.GetCellData().GetArray(0).GetValue(0)
    gauss_curv = summed_curv_data.GetPointData().GetArray(0).GetValue(0)/g_surf_area

    # have paraview recalculate the mean curvature
    curvature.SetPropertyWithName('CurvatureType','Mean')
    curvature.UpdatePipeline()
    threshold = ps.Threshold(curvature)
    threshold.Scalars = ['POINTS','Mean_Curvature']
    threshold.ThresholdRange = [-mean_filter,mean_filter]
    threshold.UpdatePipeline()
    mean_data = pv.servermanager.Fetch(threshold)
    size_mean = mean_data.GetPointData().GetArray(0).GetSize()
    summed_curv = ps.IntegrateVariables(threshold)
    summed_curv.UpdatePipeline()
    summed_curv_data = pv.servermanager.Fetch(summed_curv)
    m_surf_area = summed_curv_data.GetCellData().GetArray(0).GetValue(0)
    mean_curv = summed_curv_data.GetPointData().GetArray(0).GetValue(0)/m_surf_area

    # calculate the surface area to volume ratio
    vol = nx*dx*ny*dy*nz*dz
    g_surf_to_vol = g_surf_area/vol
    m_surf_to_vol = m_surf_area/vol

    # calculate percent of data used in threshold
    curvature_data = pv.servermanager.Fetch(curvature)
    size_curv = curvature_data.GetPointData().GetArray(0).GetSize()

    # calculate fractions of data used
    gc_frac = float(size_gauss)/float(size_curv)
    mc_frac = float(size_mean)/float(size_curv)

    # delete paraview sources and filters
    ps.Delete(summed_curv)
    ps.Delete(threshold)
    ps.Delete(curvature)
    ps.Delete(contour)
    ps.Delete(ds)
    ps.Delete(reader)
    del(sys_data)
    del(gauss_data)
    del(mean_data)
    del(curvature_data)
    del(summed_curv_data)

    return gauss_curv,mean_curv,gc_frac,mc_frac,g_surf_area,m_surf_area,vol

def calc_curvature_batch(data_directory,save_directory,output_file,arr_name,gauss_filter,mean_filter,sample_rate):
    import numpy as np
    import os
    from sys import argv
    import fnmatch

    # get list of all the vtk files in the batch
    filelist = os.listdir(data_directory)
    filelist = fnmatch.filter(filelist,'fluid*')
    filelist.sort(key=lambda x: int(x.split('.')[0].split('_')[1]))

    # create arrays for storing the data (gaussian and mean curvatures,
    # and the fractions of data used in the curvature calculations for both)
    points=len(filelist)
    gc = np.zeros(points)
    mc = np.zeros(points)
    gc_frac = np.zeros(points)
    mc_frac = np.zeros(points)
    gsurf = np.zeros(points)
    msurf = np.zeros(points)
    vol = np.zeros(points)

    # calculate gaussian and mean curvatures for each data file
    for i,filename in enumerate(filelist):
        if (i==0):
            continue
        try:
            a,b,c,d,e,f,g = calc_curvature(data_directory+'/'+filename,arr_name,gauss_filter,mean_filter,sample_rate)
            gc[i] = a
            mc[i] = b
            gc_frac[i] = c
            mc_frac[i] = d
            gsurf[i] = e
            msurf[i] = f
            vol[i] = g
        except ZeroDivisionError:
            print('division by zero. Curvatures will not be sensible!')
        except OSError:
            print('No vtk data for this system!')

    # save data
    skip = int(filelist[points-1].split('.')[0].split('_')[1])/(points-1)
    t = np.arange(0,points*skip,skip) 

    # save the ads data for replotting or manipulation later
    np.savez(save_directory+'/'+output_file+'.npz',gc,mc,gc_frac,mc_frac,t,gsurf,msurf,vol)

    return 0

def get_curvatures(vtk_file,vtk_file_dir,output_file,output_file_dir,arr_name,sample_rate,gauss_filter,mean_filter):
    import paraview.simple as ps
    import numpy as np
    import paraview as pv
    from vtk.util.numpy_support import vtk_to_numpy

    # have paraview open the vtk data file
    reader = ps.OpenDataFile(vtk_file_dir + vtk_file)
    sys_data = pv.servermanager.Fetch(reader)
    nx,ny,nz = sys_data.GetDimensions()
    dx,dy,dz = sys_data.GetSpacing()

    # downsample the data (makes for a smoother contour surface)
    ds = ps.ExtractSubset(reader)
    ds.SampleRateI = sample_rate 
    ds.SampleRateJ = sample_rate
    ds.SampleRateK = sample_rate
    ds.VOI[1]=nx-1
    ds.VOI[3]=ny-1
    ds.VOI[5]=nz-1
    ds.IncludeBoundary = 1
    ds.UpdatePipeline()

    # have paraview apply a contour surface at a concentration value of cont_val
    contour = ps.Contour(ds)
    cont_val = 0.0						 # this might change depending on order parameter
    contour.ContourBy = ['POINTS',arr_name] # Viz is the name of the vtk array
    contour.Isosurfaces = [cont_val]	 
    contour.SetPropertyWithName('ComputeNormals',0)
    contour.UpdatePipeline()

    # have paraview calculate the curvature (default is Gaussian)
    curvature = ps.Curvature(contour)

    # filter the curvatures
    threshold = ps.Threshold(curvature)
    threshold.Scalars = ['POINTS','Gauss_Curvature']
    threshold.ThresholdRange = [-gauss_filter,gauss_filter]
    threshold.AllScalars = 1
    threshold.UpdatePipeline()
    gauss_data = pv.servermanager.Fetch(threshold)
    size_gauss = gauss_data.GetPointData().GetArray(0).GetSize()

    # convert vtk array to numpy array
    gauss_curv = vtk_to_numpy(gauss_data.GetPointData().GetArray(0))

    # integrate the surface area and curvature
    summed_curv = ps.IntegrateVariables(threshold)
    summed_curv_data = pv.servermanager.Fetch(summed_curv)
    g_surf_area = summed_curv_data.GetCellData().GetArray(0).GetValue(0)

    # have paraview recalculate the mean curvature
    curvature.SetPropertyWithName('CurvatureType','Mean')
    curvature.UpdatePipeline()
    threshold = ps.Threshold(curvature)
    threshold.Scalars = ['POINTS','Mean_Curvature']
    threshold.ThresholdRange = [-mean_filter,mean_filter]
    threshold.UpdatePipeline()
    summed_curv = ps.IntegrateVariables(threshold)
    summed_curv.UpdatePipeline()
    summed_curv_data = pv.servermanager.Fetch(summed_curv)
    m_surf_area = summed_curv_data.GetCellData().GetArray(0).GetValue(0)
    mean_data = pv.servermanager.Fetch(threshold)

    # convert vtk array to numpy array
    mean_curv = vtk_to_numpy(mean_data.GetPointData().GetArray(0))

    # calculate the surface area to volume ratio
    g_surf_to_vol = g_surf_area/(nx*dx*ny*dy*nz*dz)
    m_surf_to_vol = m_surf_area/(nx*dx*ny*dy*nz*dz)

    # calculate percent of data used in threshold
    curvature_data = pv.servermanager.Fetch(curvature)
    size_curv = curvature_data.GetPointData().GetArray(0).GetSize()

    # save the numpy arrays for later manipulation
    np.savez(output_file_dir+output_file,gauss_curv,mean_curv,g_surf_to_vol,g_surf_area,m_surf_to_vol,m_surf_area)

    # delete paraview sources and filters
    ps.Delete(summed_curv)
    ps.Delete(contour)
    ps.Delete(ds)
    ps.Delete(reader)
    del(sys_data)
    del(summed_curv_data)

    return 0

def calc_surf_to_vol(filename,arr_name,sample_rate):
    import paraview.simple as ps
    import numpy as np
    import paraview as pv
    from vtk.util.numpy_support import vtk_to_numpy

    # have paraview open the vtk data file
    reader = ps.OpenDataFile(filename)
    sys_data = pv.servermanager.Fetch(reader)
    nx,ny,nz = sys_data.GetDimensions()
    dx,dy,dz = sys_data.GetSpacing()

    # downsample the data (makes for a smoother contour surface)
    ds = ps.ExtractSubset()
    ds.SampleRateI = sample_rate 
    ds.SampleRateJ = sample_rate
    ds.SampleRateK = sample_rate
    ds.VOI[1]=nx-1
    ds.VOI[3]=ny-1
    ds.VOI[4]=1     # leave off bottom layer for CHBDThinFilm
    ds.VOI[5]=nz-2  # leave off top layer for CHBDThinFilm
    ds.IncludeBoundary = 1
    ds.UpdatePipeline()

    # have paraview apply a contour surface at a concentration value of cont_val
    contour = ps.Contour()
    cont_val = 0.5	 # this might change depending on order parameter
    contour.ContourBy = ['POINTS',arr_name] # Viz is the name of the vtk array
    contour.Isosurfaces = [cont_val]	 
    contour.SetPropertyWithName('ComputeNormals',0)
    contour.UpdatePipeline()

    # integrate the surface area and curvature
    summed_curv = ps.IntegrateVariables()
    summed_curv_data = pv.servermanager.Fetch(summed_curv)
    surf_area = summed_curv_data.GetCellData().GetArray(0).GetValue(0)

    # calculate the surface area to volume ratio
    volume = nx*dx*ny*dy*nz*dz
    surf_to_vol = surf_area/volume

    # delete paraview sources and filters
    ps.Delete(reader)
    ps.Delete(ds)
    ps.Delete(contour)
    ps.Delete(summed_curv)
    del(sys_data)
    del(summed_curv_data)

    return surf_to_vol,surf_area,volume

def calc_stvr_batch(data_directory,save_directory,output_file,arr_name,sample_rate):
    import numpy as np
    import os
    from sys import argv
    import fnmatch

    # get list of all the vtk files in the batch
    filelist = os.listdir(data_directory)
    filelist = fnmatch.filter(filelist,'c2*')
    filelist.sort(key=lambda x: int(x.split('.')[0].split('_')[1]))

    # create arrays for storing the data (gaussian and mean curvatures,
    # and the fractions of data used in the curvature calculations for both)
    points=len(filelist)
    stvr = np.zeros(points)
    surf_area = np.zeros(points)
    volume = np.zeros(points)

    # calculate gaussian and mean curvatures for each data file
    for i,filename in enumerate(filelist):
        try:
            a,b,c = calc_surf_to_vol(data_directory+'/'+filename,arr_name,sample_rate)
            stvr[i] = a
            surf_area[i] = b
            volume[i] = c
        except ZeroDivisionError:
            print('division by zero. Curvatures will not be sensible!')
        except OSError as err:
            print('No vtk data for this system! '+str(err))
        except AttributeError as ae:
            print('\tencountered an attribute error: '+str(ae))
            print('\toccured while performing calc_surf_to_vol on file:')
            print('\t\t'+filename)

    # save data
    skip = int(filelist[points-1].split('.')[0].split('_')[1])/(points-1)
    t = np.arange(0,points*skip,skip) 

    # save the ads data for replotting or manipulation later
    np.savez(save_directory+'/'+output_file+'.npz',stvr,surf_area,volume,t)

    return 0
