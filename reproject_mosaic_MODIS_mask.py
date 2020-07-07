import os, osr, glob, datetime, gdal
import numpy as np

##### BOP
# #DESCRIPTION:
# Simply designed for reprojecting and mosaicing "one" specified MODIS SIN or ISIN data
# Sinusodial (SIN) informaion for MODIS data: https://modis-land.gsfc.nasa.gov/GCTP.html
# Two SIN grid systems are using for MODIS land products: MODLAND SIN and Intergerized SIN (ISIN)
#
# #REVISION HISTORY:
# 23 Mar 2020: Hyunglok Kim: Initial Specification
# 
# For MOD/MYD16:
# MOD_Grid_MOD16A2:ET_500m

input_dir = '/home/h/tmp/dataset'
output_dir = '/home/h/GoogleDrive/gdal_study'

starting_date = '2019-01-01'
ending_date = '2019-01-03'

product = 'MOD16A2GF' 
variable = 'ET_500m' # can be checked with: gdal.Open(hdf_file).GetSubDatasets()

# Run domain
Run_domain_lower_left_lat = 25.005
Run_domain_lower_left_lon = -124.995
Run_domain_upper_right_lat = 52.995
Run_domain_upper_right_lon = -67.005
Run_domain_resolution_dx = 0.01 #degree for WGS84
Run_domain_resolution_dy = 0.01 #degree for WGS84

# fill and unwanted values
fill_val=32767 #designated fill value
unwanted_vals=[32766, 32765, 32764, 32763, 32762, 32761] #define other fill value(s)

# subdata_index for MOD16A2/MYD16A2/MOD16A2GF/MYD16A2GF
#0: ET_500m
#1: LE_500m
#2: PET_500m
#3: PLE_500m
#4: ET_QC_500m
subdata_index=0

##### EOP
starting_date = datetime.datetime.strptime(starting_date, '%Y-%m-%d').date()
ending_date = datetime.datetime.strptime(ending_date, '%Y-%m-%d').date()
nof_date = (ending_date - starting_date).days + 1

date_list = []
for i in range(nof_date):
     date_str = starting_date + datetime.timedelta(days=i)
     # Make yyyydoy of the file name of MODxxXxx.Ayyyydoy
     date_list.append(str(date_str.timetuple().tm_year) + str(date_str.timetuple().tm_yday).zfill(3))

del date_str,nof_date,starting_date,ending_date

os.chdir(input_dir)
hdf_files_list = sorted(glob.glob(product+'*.hdf'))

# masking files with unwanted data / change unwanted values to fillvalue
def hdf2tif_masking(hdf_file, out_file_name, subdata_index, fill_val, unwanted_vals):
    src_ds = gdal.Open(hdf_file)
    or_data=gdal.Open(src_ds.GetSubDatasets()[subdata_index][0], gdal.GA_ReadOnly)
    or_array=or_data.ReadAsArray().astype(np.int16)

    for i in range(len(unwanted_vals)):
        or_array[np.where(or_array==unwanted_vals[i])]=fill_val

    out_ds=gdal.GetDriverByName('Gtiff').Create(out_file_name, or_data.RasterXSize, or_data.RasterYSize,1, \
                                                gdal.GDT_Int16,['COMPRESS=LZW', 'TILED=YES'])
    out_ds.SetGeoTransform(or_data.GetGeoTransform())
    out_ds.SetProjection(or_data.GetProjection())
    out_ds.GetRasterBand(1).WriteArray(or_array)
    out_ds.GetRasterBand(1).SetNoDataValue(fill_val)
    out_ds = None

os.system('mkdir '+output_dir+'/tmp')
os.chdir(output_dir+'/tmp')

daily_hdf_files = []
tmp_file_group=''
for i in range(len(date_list)):
    t_daily_files=[j for j in hdf_files_list if product+'.A'+date_list[i] in j] #serach if certain date's files exist

    if t_daily_files: #create projected temp tif files
         
         for k in range(len(t_daily_files)):
             hdf_file=input_dir+'/'+t_daily_files[k]
             tmp_file=output_dir+'/tmp/tmp.tif'
             hdf2tif_masking(hdf_file,tmp_file,subdata_index,fill_val,unwanted_vals)
             out_file_name = output_dir+'/tmp/tmp_'+str(k)+'.tif'
             #out_file_name = output_dir+'/tmp/'+t_daily_files[k][:-3]+'.tif'
             pass2gdalwarp='gdalwarp -q -overwrite -nosrcalpha -et 0.1 -t_srs EPSG:4326 -tr '+\
                           str(Run_domain_resolution_dx)+' '+str(Run_domain_resolution_dy)+' '+tmp_file+' '+\
                           out_file_name
             os.system(pass2gdalwarp)
             tmp_file_group=tmp_file_group+' tmp_'+str(k)+'.tif'
       
         pass2gdalwarp='gdalwarp -q -overwrite -nosrcalpha -et 0.1 '+tmp_file_group+' '+output_dir+'/tmp/tmp_mosaic.tif'
         os.system(pass2gdalwarp)

         pass2gdalwarp='gdal_translate -q -projwin '+str(Run_domain_lower_left_lon)+' '+str(Run_domain_upper_right_lat)+\
                       ' '+str(Run_domain_upper_right_lon)+' '+str(Run_domain_lower_left_lat)+' '+output_dir+'/tmp/tmp_mosaic.tif '+\
                       output_dir+'/'+product+'.A'+date_list[i]+'_'+variable+'.tif'
         os.system(pass2gdalwarp)

os.system('rm '+output_dir+'/tmp/tmp*')