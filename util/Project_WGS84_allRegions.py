
# import system modules
import arcpy,os
import subprocess



regions=['Aus','CAm','Ind','NAf','SAf','SAm','SAs']



CAm_gdal= 'gdal_translate -a_srs "+a=6376192.450 +rf=50000000 +proj=laea +lat_0=17.500000 +lon_0=-86.000000" -a_ullr -3200000.000000 1710549.926758 3203549.804688 -1400000.000000 -of Gtiff' # inputxxx outxxxx'
NAM_gdal= 'gdal_translate -a_srs "+a=6367415.828 +rf=50000000 +proj=laea +lat_0=45.000000 +lon_0=-92.500000" -a_ullr -4200000.000000 2813049.804688 4206049.804688 -2300000.000000 -of Gtiff'
SAM_gdal= 'gdal_translate -a_srs "+a=6375250.017 +rf=50000000 +proj=laea +lat_0=-21.500000 +lon_0=-57.500000" -a_ullr -2900000.000000 4028049.804688 2925049.804688 -4200000.000000 -of Gtiff'
SAf_gdal= 'gdal_translate -a_srs "+a=6376877.497 +rf=50000000 +proj=laea +lat_0=-14.000000 +lon_0=29.000000" -a_ullr -2700000.000000 2702549.804688 2724549.804688 -2900000.000000 -of Gtiff'
Ind_gdal= 'gdal_translate -a_srs "+a=6378094.108 +rf=50000000 +proj=laea +lat_0=-2.500000 +lon_0=129.000000" -a_ullr -4000000.000000 1421549.926758 4005549.804688 -1600000.000000 -of Gtiff'
Aus_gdal= 'gdal_translate -a_srs "+a=6373089.386 +rf=50000000 +proj=laea +lat_0=-29.000000 +lon_0=145.000000" -a_ullr -3900000.000000 2201549.804688 3927549.804688 -2600000.000000 -of Gtiff'
SAs_gdal= 'gdal_translate -a_srs "+a=6376192.450 +rf=50000000 +proj=laea +lat_0=17.500000 +lon_0=95.000000" -a_ullr -3900000.000000 1799549.926758 3927549.804688 -1400000.000000 -of Gtiff'
NAf_gdal= 'gdal_translate -a_srs "+a=6375376.559 +rf=50000000 +proj=laea +lat_0=21.000000 +lon_0=22.500000" -a_ullr -4700000.000000 2824049.804688 4729549.804688 -2200000.000000 -of Gtiff '

regions_gdal=[Aus_gdal,CAm_gdal,Ind_gdal,NAf_gdal,SAf_gdal,SAM_gdal,SAs_gdal]


for i_region in range(0, len(regions)):
    print(regions[i_region])



    region_name=regions[i_region]  

    arcpy.env.workspace = "E:\\ASCAT\\unzip_sir_aIC\\" + region_name+ "_asc"
    outworkspace="E:\\ASCAT\\unzip_sir_aIC\\" + region_name+ "_asc_prj"
    outworkspace_temp="E:\\ASCAT\\unzip_sir_aIC\\" + region_name+ "_asc_prj_temp"

    if not os.path.exists(outworkspace):
        os.mkdir(outworkspace)

    if not os.path.exists(outworkspace_temp):
        os.mkdir(outworkspace_temp)


    fcList = arcpy.ListRasters()

    for raster1 in fcList:
        print raster1


        outRaster1=outworkspace_temp+os.sep+raster1
        outRaster2=outworkspace+os.sep+raster1[:-4]+'.tif'
        if os.path.exists(outRaster2):
            continue

        fullCmd=regions_gdal[i_region]+' '+arcpy.env.workspace+os.sep+raster1+' '+outRaster1 ### change here

        p = subprocess.Popen(fullCmd,shell=True)
        p.wait()


        
        arcpy.ProjectRaster_management(outRaster1, outRaster2,\
                                   r"D:\Amazon\GCS_WGS_1984.prj", "#", "0.044799895",\
                                   "#", "#", "#")

        arcpy.Delete_management(outRaster1)
        arcpy.Delete_management(raster1)
