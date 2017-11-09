These are instructions for running things on HPC:
1. To submit scripts: "qsub <script_name>"
2. To delete scripts: "qdel <job_id>"
3. To see status: "qstat -u <HPC_userID>"

A note on scripts:
1. script_{i}.sh: PW analysis on ith DREAM4 dataset of size N, where N is set in master scripts in ./btp/pairwise
2. script_ige_{type}_{i}.sh: IGE analysis on ith DREAM4 dataset of size N, where N is set in master scripts in ./btp/global, and "type" is {'' : regular IGE, 'new': soft IGE, 'pag': for pageranked IGE, 'jug': for iterative pagerank IGE}

Correspondingly, there are also "output" and "errorLog" files to aid in debugging.

Note: scripts assume that ./btp is in home folder. If not, make changes to scripts aptly.