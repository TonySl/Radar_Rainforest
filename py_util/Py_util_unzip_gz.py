### Shengli Tao, sltao1990@gmail.com

import gzip
import glob
import os


root_quick=r"E:\ASCAT"
root_out=r'E:\ASCAT\unzip_sir_aIC\SAm'


if not os.path.exists(root_out):
    os.mkdir(root_out)

for iyear in range(2007,2019):
##    print str(iyear)
    root_2_days=root_quick+os.sep+str(iyear)+os.sep+'SAm'
##    root_2_days=root_quick #+os.sep+'sir\qaeh\NAm'
##    print root_2_days
    
    for roots, dirs, files in os.walk(root_2_days):

        for file1 in files:
            if file1.endswith("sir.gz") and file1.startswith("mafa-a-"):  
##                print roots+os.sep+file1
                
                inputname=roots+os.sep+file1  

##                outname=inputname[:-3]
                
                outname=root_out+os.sep+file1[:-3]
                
                if os.path.exists(outname):
                    continue
                
                print outname
                f=gzip.open(inputname, 'rb')
                file_content = f.read()
                with open(outname,'wb') as sir:
                    sir.write(file_content)

##            if file1.endswith("sir.gz") and file1.startswith("mafa-I-"):  ## incident angle std
##                print roots+os.sep+file1
##                
##                inputname=roots+os.sep+file1  
##
####                outname=inputname[:-3]
##                
##                outname=root_out+os.sep+file1[:-3]
##                print outname
##                if os.path.exists(outname):
##                    continue
##
##                f=gzip.open(inputname, 'rb')
##                file_content = f.read()
##                with open(outname,'wb') as sir:
##                    sir.write(file_content)

                    
            if file1.endswith("sir.gz") and file1.startswith("mafa-C-"):  ## angle count file
##                print roots+os.sep+file1
                
                inputname=roots+os.sep+file1 

##                outname=inputname[:-3]
                
                outname=root_out+os.sep+file1[:-3]
                
                if os.path.exists(outname):
                    continue
                
                print outname
                f=gzip.open(inputname, 'rb')
                file_content = f.read()
                with open(outname,'wb') as sir:
                    sir.write(file_content)
