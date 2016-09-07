Datasets
-----
The datasets HiSeq and MiSeq are included in Kraken's "Timing data" at:
http://ccb.jhu.edu/software/kraken/

The datasets wgsim200 and wgsim250 can be downloaded at:
https://seafile.rlp.net/f/ab4ac98d55/


Example light database
-----
An example light database containing bacterial species-level targets can be downloaded at:
https://seafile.rlp.net/f/f2c9aa6fd6/

If CuCLARK was installed and the database extracted to `<directoryDB>` the classification can be executed with:
```
./CuCLARK-l -T <directoryDB/targets_0.txt> -D <directoryDB/> -O <fileObjects> -R <fileResults> -n <numberofthreads>
```
