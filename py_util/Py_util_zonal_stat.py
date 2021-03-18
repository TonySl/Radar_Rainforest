### Shengli Tao, sltao1990@gmail.com

# Import system modules
import arcpy,os
from arcpy import env
from arcpy.sa import *


# set dir
dir_SAm = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\SAm_asc_NoStrip_C20_prj"
dir_CAm = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\CAm_asc_NoStrip_C20_prj"
dir_SAf = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\SAf_asc_NoStrip_C20_prj"
dir_SAs = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\SAs_asc_NoStrip_C20_prj"
dir_Ind = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\Ind_asc_NoStrip_C20_prj"
dir_Aus = r"D:\Drought\ASCAT\monthly_layers\unzip_sir_aIC\Aus_asc_NoStrip_C20_prj"

all_dir=[dir_SAm,dir_CAm,dir_SAf,dir_SAs,dir_Ind,dir_Aus]


# Set finshet
inZoneData_SAm = r"D:\Drought\shps\Americas_fishnet_SAm.shp" # 0.25 degree
inZoneData_CAm = r"D:\Drought\shps\Americas_fishnet_CAm.shp" # 0.25 degree
inZoneData_SAf = r"D:\Drought\shps\Africa_fishnet.shp" # 0.25 degree
inZoneData_SAs = r"D:\Drought\shps\Asia_fishnet_SAs.shp" # 0.25 degree
inZoneData_Ind = r"D:\Drought\shps\Asia_fishnet_Ind.shp" # 0.25 degree
inZoneData_Aus = r"D:\Drought\shps\Asia_fishnet_Aus.shp" # 0.25 degree

all_fish=[inZoneData_SAm,inZoneData_CAm,inZoneData_SAf,inZoneData_SAs,inZoneData_Ind,inZoneData_Aus]

zoneField = "FID_"


for i in range(0,len(all_dir)):
    print(all_dir[i])
    env.workspace = all_dir[i]
    out_workspace = all_dir[i]+"_zonal"

    if not os.path.exists(out_workspace):
        os.makedirs(out_workspace)


    inZoneData=all_fish[i]
    
    rasters=arcpy.ListRasters()

    for inValueRaster in rasters:
        print inValueRaster

    # Check out the ArcGIS Spatial Analyst extension license
        arcpy.CheckOutExtension("Spatial")

    # Execute ZonalStatistics
        
        basename = os.path.basename(inValueRaster)
        

        outname=outname=out_workspace+os.sep+ basename.replace('.tif', '')+".dbf"


        ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster, 
                                     outname, "DATA", "MEAN")
