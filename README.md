# Codes for Decoing Attended Spatial Location during Complex Scene Analysis with fNIRS


[![Stable Version](https://img.shields.io/badge/fNIRS-v1.0.0-blue)](https://github.com/NoPenguinsLand/fNIRSCodes_Manuscript)
[![Build Status](https://github.com/NSNC-Lab/fNIRS_DecodingSpatialAttention/actions/workflows/github-actions-demo.yml/badge.svg)](https://github.com/NSNC-Lab/fNIRS_DecodingSpatialAttention/actions)
[![Mirror](https://img.shields.io/badge/Mirror-Bitbucket-blue)](https://bitbucket.org/nopenguinsland/fnirscodes_manuscript_mirror/src/master/)

 For "Decoding Attended Spatial Location during Complex Scene Analysis with fNIRS."

**Authors**: Matthew Ning, Meryem Yücel, Alexander Von Lühmann, David A Boas, Kamal Sen.

**Affliation**: Neurophotonics Center, Department of Biomedical Engineering, Boston University.

**Preprint**: https://www.biorxiv.org/content/10.1101/2022.09.06.506821v2

**Citation**: 

> ##### Ning, M., Yücel, M. A., Von Lühmann, A., Boas, D. A. & Sen, K. Decoding Attended Spatial Location during Complex Scene Analysis with fNIRS. bioRxiv (2022) doi:10.1101/2022.09.06.506821.

**Publication**: Coming soon.

**Dataset**: [Google Drive](https://drive.google.com/drive/folders/1MWttV66DmfCk02pftK2U-wSJxp2okwvN?usp=sharing)

**Abstract**: When analyzing complex scenes, humans often focus their attention on an object at a particular spatial location. The ability to decode the attended spatial location would facilitate brain computer interfaces for complex scene analysis. Here, we investigated functional near-infrared spectroscopy’s (fNIRS) capability to decode audio-visual spatial attention in the presence of competing stimuli from multiple locations.  We targeted dorsal frontoparietal network including frontal eye field (FEF) and intra-parietal sulcus (IPS) as well as superior temporal gyrus/planum temporal (STG/PT). They all were shown in previous fMRI studies to be activated by auditory, visual, or audio-visual spatial tasks. We found that fNIRS provides robust decoding of attended spatial locations for most participants and correlates with behavioral performance. Moreover, we found that FEF makes a large contribution to decoding performance. Surprisingly, performance was significantly above chance level ~1s, well before the peak of the fNIRS response. Our results demonstrate that fNIRS is a promising platform for a compact, wearable technology that could be applied to decode attended spatial location and reveal contributions of specific brain regions during complex scene analysis.

**Running Instructions**:
* To run everything in one click, run the runAllPreprocessing.m script within the SbjLvlProcessing folder.
    * **WARNING**: it'll take 1-3 days depending on your computer. Parallel computing version is available at [Parallel Computing Version](#parallel-computing-instructions)

* To run preprocessing pipeline for Figure 5 and 6, in the SbjLvlProcessing folder, run the runPreprocessing_Main_2Class.m script for 2-class classification and similarly, runPreprocessing_Main_3Class.m script for the 3-class classification.
    * Estimated running time: up to 8-20 hours.

* To run preprocessing pipeline for Figure 3, 4, 7 & 8, in the SbjLvlProcessing folder, run the runPreprocessing_Others.m script.
    * Estimated running time: 1-4 hours.

* To generate Figure 3-8, in the SbjLvlProcessing folder, run the genFigures.m script.

# Parallel Computing Instructions

**List of directories for binaural lab:**
* /usr3/graduate/mhn (smaller storage quota, for small personal project).
* /projectnb/binaural/mhn (larger storage quota, shared among lab directory, for bigger project. where I run parallel computing).

Before running batch file, run below command:

```dos2unix *.sh *.m```

Then check status of file using below command:

```file runBatchfNIRSJob.m```

Should be ASCII text executable

Finally run the following command:

```qsub -pe omp 8 -l h_rt=28:00:00 ./runfNIRSBash.sh```

You can vary the parameter values if you know what you're doing.

Helpful tips:

To view current status of job:

```qstat -u mhn```

To view outputs:
view runfNIRSBash.sh.o(job_number)
view runfNIRSBash.sh.po(job_number)

To view error messages:
runfNIRSBash.sh.e(job_number)
view runfNIRSBash.sh.pe(job_number)

vim tip/commands:
go to top of file: ESC + gg
go to bottom of file: SHIFT + g
go to line: in view mode, enter number
to save and quit: ":wq"
to quit without saving: ":q" or ":q!"
to go to edit mode: type "i"

To run MATLAB without desktop:

```module load matlab
matlab -nodisplay```

To view/manage storage quota:
scc-ondemand.bu.edu

To debug MATLAB without desktop but from command line:
https://www.mathworks.com/help/releases/R2019b/matlab/debugging-code.html

Advanced options for sun grid engine:
To view a list of parallel computing environments and their configurations, type qmon, then when main GUI window pops up,
click Parallel Environment Configuration.

For description of different parallel environments on BU SCC, refer to [4]

To run multicores, use -pe option in qsub command.
To submit array of jobs, use -t option in qsub command.

For full list of options for SGE, refer to [3].

References & Additional Readings:
[1] https://docs.oracle.com/cd/E19957-01/820-0698/6ncdvjcmd/index.html
[2] https://www.bu.edu/tech/support/research/software-and-programming/common-languages/matlab/matlab-batch/
[3] https://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
[4] https://www.bu.edu/tech/support/research/system-usage/running-jobs/parallel-batch/#pe
