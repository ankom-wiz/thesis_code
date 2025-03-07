#!/home/ankom/gimathesis/thesis_code/pygima/bin/python3
# the Shebang statement above is only needed when this script is called from bash 

import os
import requests
import tarfile
import requests

def main():
    """
    Main function to be called
    """
    dataroot='data'
    nmeadir=os.path.join(dataroot,'nmea')
    if not os.path.exists(nmeadir):
        os.makedirs(nmeadir)
    nmeasources={
        "jinja_nmea":
        {"url":"https://surfdrive.surf.nl/files/index.php/s/mDyKTExFFcAZlZI/download","file":"Jinja_Spring_summer2024.tgz"},
    }
    
    for name,val in nmeasources.items():
    
        fout=os.path.join(nmeadir,val['file'])
        if not os.path.exists(fout):
            print(f"Downloading NMEA data {fout} for {name}")
            r=requests.get(val['url'])
    
            with open(fout,'wb') as fid:
                fid.write(r.content)
    #extract ijinja data in subdir
    

    # Extract and prepare dataq for analysis
    jinja=nmeasources['jinja_nmea']
    archive=os.path.join(nmeadir,jinja['file'])

    jinja_nmeadir=os.path.join(nmeadir,"jinja")
    if not os.path.exists(jinja_nmeadir):
        os.makedirs(jinja_nmeadir)
    
    with tarfile.open(archive,'r:gz') as tf:
        for nmeaf in tf.getmembers():
            if nmeaf.isreg():
                nmeaf.name = os.path.basename(nmeaf.name) 
                print(f"extracting {jinja_nmeadir}/{nmeaf.name}")
                tf.extract(nmeaf,jinja_nmeadir) # extract



if  __name__ == "__main__":
    # This part is executed when this script is called as a script (e.g. python thisscript.py)

    # call the main function as defined above
    main()
