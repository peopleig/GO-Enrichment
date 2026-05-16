# GO-Enrichment

The usage of `step2_scan_motif.py` took place because `step2_extract_and_dreg.sh` was taking too long, almost 3 hours  
So I sought to parallelize the process, as the `.sh` script was only using a single core  
Hence, there's two files for step2, to show my 2 attempts and why the second one took like 20s to iterate through the list and the first one took 3 hours