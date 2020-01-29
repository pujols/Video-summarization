-- Environment check --
After download this package, you will have 3 folders---codes, data, Evaluation---under /seqDPP


-- Installation --
1) Please download the features from
https://www.dropbox.com/s/bpy41o2zglk4ka4/video_summarization.zip?dl=0
Unzip it and put the "video_summarization" folder under "/seqDPP"

2) Please download the downsampled video frames from
https://drive.google.com/open?id=1tKc-yxMAZKDaTNw_AIpIHTst5TZElIK7
Unzip it and put the "Frames_sampled" folder under "/seqDPP/data"
Please put all the files (e.g, vXX.mat) originally within "/seqDPP/data/Frames_sampled/OVP or /Youtube" into "/seqDPP/data/Frames_sampled"
By doing so, you should have 100 vXX.mat (XX from 11 to 110) under "/seqDPP/data/Frames_sampled"

3) Please download minFunc
https://www.cs.ubc.ca/~schmidtm/Software/minFunc_2012.zip
Unzip it and change the folder name "minFunc_2012" to "minFunc"
Put the "minFunc" folder under "/seqDPP/codes"

4) Please download the user summaries
https://www.dropbox.com/s/ilt1jpclzs2o18v/UserSummary.zip
https://www.dropbox.com/s/7rtbyeo64hk8ot7/newUserSummary.zip
Unzip them and put all the folders (e.g., vXX) originally within "UserSummary" and "newUserSummary" into "/seqDPP/data/OVP_YouTube_cmp"
!!!! Please remove all the .eps files in "v37". !!!!
By doing so, you should have 100 vXX folders (XX from 11 to 110) under "/seqDPP/data/OVP_YouTube_cmp"

5) Please download the VSUMM evaluation code
https://sites.google.com/site/vsummsite/download/CUS.jar?attredirects=0
Please change the name of "CUS.jar" to "CUS1.jar"
Please put "CUS1.jar" under "seqDPP/Evaluation"

6) Please change the line 17 in /seqDPP/Evaluation/comparision_seqDPP.m from "/usr/lib/jvm/jre-1.6.0/bin/java" to the java path of yours.


-- To begin --
Try seqDPP/codes/demo.m


-- Advance --
To speedup, please change L97, L228 in seqDPP_NN_all and L97, L230 in seqDPP_linear_all from "for" to "parfor"


-- Reporting problems --
If you find any problem with the code, please contact the authors via weilunchao760414@gmail.com


