from osgeo import gdal
import fiona
import numpy as np
import pickle
from pyproj import Proj, CRS, Transformer
import os
import ogr
from osgeo import osr

import gzip
import xarray as xr


def extract_thickness(gdirs,data):
    #This fonction will crop Millan tif file (if it doesnt already exist) according the the glaciers coordinates
    #Then it will extract the thicknesses on the flowlines and put it in a .pkl file
	#gdir_rgi : list of str
	#glacier RGI ids
	#data : str
	#data source (farinotti or millan)
    for gdir_rgi in gdirs:
    
        gdir_region=gdir_rgi[6:8]
        gdir_subreg=gdir_rgi[6:11]
        
        #Reading the corresponding Farinotti geotif
        a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Farinotti/RGI60-'+gdir_region+'/'+gdir_rgi+'_thickness.tif')
        
    	#projection info 
        proj=a.GetProjection()
        temp=proj[-8::]
        projstr=temp[0:5]
        zone=projstr[3:]
        orient=projstr[2]
        right_file=''
    
        if (orient=='6'): #326
            orientation='north'
        else: #327
            orientation='south'
            
    	#projection definition (from lon,lat to right UTM)
        myProj = Proj("+proj=utm +zone="+zone+" +"+orientation+" +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
        
        #crop Millan data region from Farinotti geotif
        if (data=='millan'): 
            
             #Adjustment for the 13-15 regions
            if (gdir_region in ['13','14','15']):
                gdir_region_m='13-15'
            else :
                gdir_region_m=gdir_region
                    
            #Checking if the millan glacier geotif already exists or not
            exist = ''+gdir_rgi+'_thickness_m_from_f.tif'
            list_exist = [_ for _ in os.listdir('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/') if _.endswith(exist)]
            
            if (list_exist == []):
                #Coordinates of Farinotti file
                ulx, xres, xskew, uly, yskew, yres = a.GetGeoTransform()
                lrx = ulx + (a.RasterXSize * xres)
                lry = uly + (a.RasterYSize * yres)
                
                
                #Selecting all the tif files in the region directory
                file_ext='.tif'
                list_file = [_ for _ in os.listdir('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/') if _.endswith(file_ext)]
                print(list_file)
                
                for file_m in list_file :
                    tr=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/'+file_m+'')
                    tr_proj=tr.GetProjection()
                    tr_temp=tr_proj[-8::]
                    tr.projstr=tr_temp[0:5]
                    print(file_m)
                    print(tr.projstr)
                
                    #reading the files coordinates and projection
                    tr.ulx, tr.xres, tr.xskew, tr.uly, tr.yskew, tr.yres = tr.GetGeoTransform()
                    tr.lrx = tr.ulx + (tr.RasterXSize * tr.xres)
                    tr.lry = tr.uly + (tr.RasterYSize * tr.yres)
                    print(tr.ulx,tr.lry,tr.lrx,tr.uly)
                    
                    if (projstr!=tr.projstr):
                        #converting the coord in the same projection as the farinotti file if not the case
                        transformation=Transformer.from_crs(int(tr.projstr),int(projstr))
                        tr.xmin,tr.ymin,tr.xmax,tr.ymax = transformation.transform(tr.ulx,tr.lry,tr.lrx,tr.uly)
                    else :
                        tr.xmin,tr.ymin,tr.xmax,tr.ymax = tr.ulx, tr.lry, tr.lrx, tr.uly
        
                    print(tr.xmin,tr.ymin,tr.xmax,tr.ymax)
                    if (tr.xmin <= ulx) and (tr.ymin <= lry) and (tr.xmax >= lrx) and (tr.ymax >= uly): #the right file contains the farinotti file
                        right_file=file_m
                        
                    #cropping the millan file to create a new one centered on the glacier    
                    os.system('gdalwarp -t_srs EPSG:'+str(projstr)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' -tr '+str(xres)+' '+str(yres)+' -r bilinear /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/'+right_file+' /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
        		
                    #checking if theright file was chosen
                    test_a = gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
                    test_b = test_a.GetRasterBand(1)
                    test_H= test_b.ReadAsArray()
                    if (int(np.sum(test_H))!=0):
                        continue 
                        #if there is ice in the file chosen, we now assumed it's the right one and we don't have to go through the others
                    else:
                        #if it's not the right one, we have to remove it 
                        os.remove('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
                    
                
            a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
            
        b=a.GetRasterBand(1)
        H=b.ReadAsArray()

        
    	#Reads the shapefile coordinates
        path_to_shp='/home/lucillegimenes/oggm-workflow-b/per_glacier/RGI60-'+gdir_region+'/RGI60-'+gdir_subreg+'/'+gdir_rgi+'/outflow.shp'
        
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
            print(len(coord))
            xps1=np.zeros(len(coord))
            yps1=np.zeros(len(coord))
            
            for i in range(0,len(coord)):
    			#converting coordinates in lat & lon to coordinates in relevant UTM projection
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
            
        if (int(np.sum(H))==0):
            print('Problème avec choix ou découpage du geotif !')
            

               
    
def extract_tif_millan(gdirs):
    #This fonction only crop Millan tif file (if it doesnt already exist) according the the glaciers coordinates
    # MAYBE TO PUT AS AN OPTION IN extract_thickness ????????????
	#gdir_rgi : list of str
	#glacier RGI ids
	#data : str
	#data source (farinotti or millan)
    for gdir_rgi in gdirs:
    
        gdir_region=gdir_rgi[6:8]
        gdir_subreg=gdir_rgi[6:11]
        
        #Reading the corresponding Farinotti geotif
        a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Farinotti/RGI60-'+gdir_region+'/'+gdir_rgi+'_thickness.tif')
        
    	#projection info 
        proj=a.GetProjection()
        temp=proj[-8::]
        projstr=temp[0:5]
        zone=projstr[3:]
        orient=projstr[2]
        right_file=''
    
        if (orient=='6'): #326
            orientation='north'
        else: #327
            orientation='south'
            
    	#projection definition (from lon,lat to right UTM)
        myProj = Proj("+proj=utm +zone="+zone+" +"+orientation+" +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
        
        #crop Millan data region from Farinotti geotif
            
             #Adjustment for the 13-15 regions
        if (gdir_region in ['13','14','15']):
            gdir_region_m='13-15'
        else :
            gdir_region_m=gdir_region
                    
            #Checking if the millan glacier geotif already exists or not
        exist = ''+gdir_rgi+'_thickness_m_from_f.tif'
        list_exist = [_ for _ in os.listdir('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/') if _.endswith(exist)]
            
        if (list_exist == []):
                #Coordinates of Farinotti file
            ulx, xres, xskew, uly, yskew, yres = a.GetGeoTransform()
            lrx = ulx + (a.RasterXSize * xres)
            lry = uly + (a.RasterYSize * yres)
            print(ulx,lry,lrx,uly)
                
                #Selecting all the tif files in the region directory
            file_ext='.tif'
            list_file = [_ for _ in os.listdir('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/') if _.endswith(file_ext)]
            print(list_file)
                
            for file_m in list_file :
                tr=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/'+file_m+'')
                tr_proj=tr.GetProjection()
                tr_temp=tr_proj[-8::]
                tr.projstr=tr_temp[0:5]
                print(file_m)
                print(tr.projstr)
                
                    #reading the files coordinates and projection
                tr.ulx, tr.xres, tr.xskew, tr.uly, tr.yskew, tr.yres = tr.GetGeoTransform()
                tr.lrx = tr.ulx + (tr.RasterXSize * tr.xres)
                tr.lry = tr.uly + (tr.RasterYSize * tr.yres)
                print(tr.ulx,tr.lry,tr.lrx,tr.uly)
                    
                if (projstr!=tr.projstr):
                        #converting the coord in the same projection as the farinotti file if not the case
                    transformation=Transformer.from_crs(int(tr.projstr),int(projstr))
                    tr.xmin,tr.ymin,tr.xmax,tr.ymax = transformation.transform(tr.ulx,tr.lry,tr.lrx,tr.uly)
                else :
                    tr.xmin,tr.ymin,tr.xmax,tr.ymax = tr.ulx, tr.lry, tr.lrx, tr.uly
        
                print(tr.xmin,tr.ymin,tr.xmax,tr.ymax)
                if (tr.xmin <= ulx) and (tr.ymin <= lry) and (tr.xmax >= lrx) and (tr.ymax >= uly): #the right file contains the farinotti file
                    right_file=file_m
                        
                    #cropping the millan file to create a new one centered on the glacier    
                os.system('gdalwarp -t_srs EPSG:'+str(projstr)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' -tr '+str(xres)+' '+str(yres)+' -r bilinear /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/'+right_file+' /home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
        		
                    #checking if theright file was chosen
                test_a = gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
                test_b = test_a.GetRasterBand(1)
                test_H= test_b.ReadAsArray()
                if (int(np.sum(test_H))!=0):
                    continue 
                        #if there is ice in the file chosen, we now assumed it's the right one and we don't have to go through the others
                else:
                        #if it's not the right one, we have to remove it 
                    os.remove('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
                    




def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array,projnum):
    #create a raster file from an array
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(projnum)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()                 
                    
    
    
def extract_new_millan(gdirs,elev_file,year):
    #This fonction will 1. crop the hugonnet file according to the glaciers coordinates
    #2. Create a new tif file of the thicknesses for 2000 (or 2010) using the elevation change from hugonnet
	#gdir_rgi : list of str
	#glacier RGI ids
	#elev_file : str
    #path of the elevation change file .tif
    #returns a tif file updated with elevation change from hugonnet
    #year : str
    #year of starting time (2000 or 2010)
    for gdir_rgi in gdirs:
    
        gdir_region=gdir_rgi[6:8]
        gdir_subreg=gdir_rgi[6:11]
        
              #Adjustment for the 13-15 regions
        if (gdir_region in ['13','14','15']):
            gdir_region_m='13-15'
        else :
            gdir_region_m=gdir_region
        
        #Reading the corresponding Millan geotif
        a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
        
    	#projection info 
        proj=a.GetProjection()
        temp=proj[-8::]
        projstr=temp[0:5]
        zone=projstr[3:]
        orient=projstr[2]

    
        if (orient=='6'): #326
            orientation='north'
        else: #327
            orientation='south'
            
    	#projection definition (from lon,lat to right UTM)
        myProj = Proj("+proj=utm +zone="+zone+" +"+orientation+" +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
        
        #Coordinates of the file
        ulx, xres, xskew, uly, yskew, yres = a.GetGeoTransform()
        lrx = ulx + (a.RasterXSize * xres)
        lry = uly + (a.RasterYSize * yres)
    
        #Step 1 : crop hugonnet file acording the the millan file, so that it's centered on the glacier   
        os.system('gdalwarp -t_srs EPSG:'+str(projstr)+' -te '+str(ulx)+' '+str(lry)+' '+str(lrx)+' '+str(uly)+' -tr '+str(xres)+' '+str(yres)+' -r bilinear '+elev_file+' /home/lucillegimenes/Bureau/Elevation_change/RGI-'+gdir_region_m+'/'+gdir_rgi+'_elev_ch.tif')
        		
        #Step 2 : create a new .tif file
        a=gdal.Open('/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_m_from_f.tif')
        b=a.GetRasterBand(1)
        M=b.ReadAsArray() #array from Millan file
        
        fM=M #future file with thickness from 2000
        hugo=gdal.Open('/home/lucillegimenes/Bureau/Elevation_change/RGI-'+gdir_region_m+'/'+gdir_rgi+'_elev_ch.tif')
        bh=hugo.GetRasterBand(1)
        H=bh.ReadAsArray() #array from Hugonnet file
        
        fM=M-(H*(2019-int(year)))
            
        rasterOrigin=(ulx,uly)
        newRasterfn='/home/lucillegimenes/Bureau/THICKNESS_Millan/RGI-'+gdir_region_m+'/per_glacier/'+gdir_rgi+'_thickness_'+year+'.tif'
        pixelWidth=xres
        pixelHeight=yres
        #creating a new raster file from the array
        projnum=int(projstr)
        array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,fM,projnum)
  