from osgeo import gdal
import fiona
import numpy as np
import pickle
from pyproj import Proj
import os

def extract_thickness(gdir_rgi,data):
	#gdir_rgi : str
	#glacier RGI id
	#data : str
	#data source (farinotti or millan)
	
	#Reading the geotiff with gdal
    gdir_region=gdir_rgi[6:8]

    gdir_subreg=gdir_rgi[6:11]
    print(type(gdir_region))
    a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Farinotti/RGI60-'+gdir_region+'/'+gdir_rgi+'_thickness.tif')
    
	#projection info 
    proj=a.GetProjection()
    temp=proj[-8::]
    projstr=temp[0:5]
    zone=projstr[3:]
    print(zone)
    orient=projstr[2]

    if (orient=='6'): #326
        orientation='north'
    else: #327
        orientation='south'
        
	#projection definition
    myProj = Proj("+proj=utm +zone="+zone+" +"+orientation+" +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    
    #crop Millan data region from Farinotti regions
    if (data=='millan'): 
        ulx, xres, xskew, uly, yskew, yres = a.GetGeoTransform()
        lrx = ulx + (a.RasterXSize * xres)
        lry = uly + (a.RasterYSize * yres)

        os.system('gdalwarp -t_srs EPSG:'+str(projstr)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' -tr '+str(xres)+' '+str(yres)+' -r bilinear /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region+'/THICKNESS_RGI-'+gdir_region+'_2021July09.tif /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region+'/'+gdir_rgi+'_thickness_m_from_f.tif')
		
        a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region+'/'+gdir_rgi+'_thickness_m_from_f.tif')
	
    b=a.GetRasterBand(1)
    H=b.ReadAsArray()
    
	#Read a shapefile coordinates
    path_to_shp='/home/lucillegimenes/oggm-workflow-b/per_glacier/RGI60-'+gdir_region+'/RGI60-'+gdir_subreg+'/'+gdir_rgi+'/outflow.shp'
    print(path_to_shp)
    
    lines=fiona.open(path_to_shp)
    nb_lines=len(lines)  #several flowlines
    all_hflowline=[0]*nb_lines #pre allocate
    path_to_shp='/home/lucillegimenes/oggm-workflow-b/per_glacier/RGI60-'+gdir_region+'/RGI60-'+gdir_subreg+'/'+gdir_rgi+'/outflow.shp'
    print(path_to_shp)
    
    lines=fiona.open(path_to_shp)
    nb_lines=len(lines)  #several flowlines
    all_hflowline=[0]*nb_lines #pre allocate

	#Useful information on the geotif image
    geotransform=a.GetGeoTransform()
    xmin=geotransform[0] #x coordinate of the lower left corner of the geotiff (x origin)
    ymax=geotransform[3] #y coordinate of the upper left corner of the image (y origin)
    posting = geotransform[1] #pixelwidth
    size=np.shape(H)
    ny= size[0]
    nx=size[1]
    print(nb_lines)
    
    for n in range(0,nb_lines) : #loop on the number of flowlines
    
        coord=lines[n]['geometry']['coordinates']
	#Useful information on the geotif image
    geotransform=a.GetGeoTransform()
    xmin=geotransform[0] #x coordinate of the lower left corner of the geotiff (x origin)
    ymax=geotransform[3] #y coordinate of the upper left corner of the image (y origin)
    posting = geotransform[1] #pixelwidth
    size=np.shape(H)
    ny= size[0]
    nx=size[1]
    print(nb_lines)
    
    for n in range(0,nb_lines) : #loop on the number of flowlines
    
        coord=lines[n]['geometry']['coordinates']
        print(len(coord))
        xps1=np.zeros(len(coord))
        yps1=np.zeros(len(coord))
        
        for i in range(0,len(coord)):
			#converting coordinates in lat & lon to coordinates in UTM 32N projection
            xps1[i], yps1[i] = myProj(lines[n]['geometry']['coordinates'][i][0],lines[n]['geometry']['coordinates'][i][1])
        
        newxps=np.copy(xps1)
        newyps=np.copy(yps1)
		#print(lines[n]['geometry']['coordinates'][i][0],lines[n]['geometry']['coordinates'][i][1])
		
        xnew=(newxps-xmin)/posting 
        ynew=(ymax-newyps)/posting 
        xnew=xnew.astype(int) #conversion to integer
        ynew=ynew.astype(int)
        print(xnew,ynew)
		    
		#This steps remove flowline points outside of the image (maybe not necessarry)
        w=np.where(xnew<nx) #nx is the number of pixel in the x direction, ie the number of columns, ny is the number of pixel in y direction ie the number of lines
        xnew=xnew[w]
        ynew=ynew[w]
        xps=newxps[w]
        yps=newyps[w]
        w1=np.where(ynew<ny)
        xnew=xnew[w1]
        ynew=ynew[w1]
        xps=xps[w1]
        yps=yps[w1]
        #print(xnew,ynew)
        #Extract the thickness value along the flowline
        hflowline=H[ynew,xnew] 
        all_hflowline[n]=hflowline

    path_to_pkl='/home/lucillegimenes/oggm-workflow-b/per_glacier/RGI60-'+gdir_region+'/RGI60-'+gdir_subreg+'/'+gdir_rgi+'/db_'+data+'.pkl'
    with open(path_to_pkl, 'wb') as db_file:
        pickle.dump(all_hflowline, file=db_file)

    
