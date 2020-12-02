# pangeo-AMPS is a repo for working with AMPS output (i.e. raw WRF output files) in a cloud environment.

Very early stages - currently just notebooks for creating and uploading zarr files directly to ldeo-glaciology's Google Cloud Storage bucket

---
### AMPS grids

In Stage 1, ldeo-glaciology will be working primarily with two AMPS domains, the Ross sector (domain 3) and the Antarctic Peninsula (domain 6).  We anticipate the larger pan-Antarctic domain (d02) in the near future.  

![AMPS grids](https://www2.mmm.ucar.edu/rt/amps/information/configuration/maps_2017101012/d1_colorfill_nests.png)

(reminder, WRF uses a staggered grid to aid in spatial discretization, in this case the Arikawa-C grid)
Domain | Current Res (km) | dims
------------ | ------------- | -------------
Ross (d03) | 2.67 | 675 x 1035
AP (d06) | 2.67 | 651 x 570

More info on these domains can be found at Kevin's [AMPS site](https://www2.mmm.ucar.edu/rt/amps/information/configuration/maps_2017101012/maps.html)

---

### Moving AMPS output from NCAR to the bucket

1. Log in to dataaccess node.

2. Retrieve data from archive
    
    2.1 HPSS
    
    List files on HPSS
    
    `hsi ls -l /AMPSRT/`
    
    either manually `` or write/execute script to pull from tape archive
    
    `./test_hsiget_batch.sh > & log.txt &`
    
    check what you did:
    
    `ll wrfout_d06_2019090* | du`
    `ll wrfout_d06_2019090* | du -h`
    
3. Transfer to Google Bucket

    Install gsutil
    
        I just follow this [link](https://cloud.google.com/storage/docs/gsutil_install)
    
    `curl https://sdk.cloud.google.com | bash`
    
    Add sdk directory to your path (cshell example), and restart shell:
    
    `setenv PATH $PATH\:$HOME/bin\:$HOME/google-cloud-sdk/bin\:.`
    
    `exec $SHELL -l`
    
    Initialize
    
    `gcloud init`
    
    3.3 Transfer using rsync
    
    Options: -m for multithreaded, -r recursive, -d directory, -n for a dry-run (USEFUL!)
    
    `gsutil -m rsync -r WRF_24 gs://ldeo-glaciology/AMPS/WRF_24/ > & log.txt &`
    
    or be a little more crafty (-x to exclude a pattern)...
    
    `gsutil -m rsync -n -x "wrfout_d06*$\" . gs://ldeo-glaciology/AMPS/WRF_24/ > & log.txt &`
    
4. Check the bucket for new files
    
    `gsutil ls gs://ldeo-glaciology/AMPS`
    
    

