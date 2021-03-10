import ROOT,os
#from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
#from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
year = 2017

inputTreeName   = "Ana/passedEvents"
t2_prefix    	= "/cmsuf/data"
if year==2016: 
	inputDir       = t2_prefix+'/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data_80X_2lskim_M17_Feb02/'	#2016
	outputDir      = t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/data2016/new/"
if year==2017: 
	inputDir       = t2_prefix+"/store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/"			#2017
	outputDir      = t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/data2017/new_dimu/"
	#outputDir	   = t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/no_cuts/" 
if year==2018: 
	inputDir       = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data2018_102X_M19_3l_2018jets/"	#2018
	outputDir      = t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/data2018/new/"
if year==2000:
	inputDir	   = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/Data2018_102X_M19_3l_2018jets/"
	outputDir	   = t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/test/"
	fileNames = ["skimZ1L_Data_Run2017-17Nov2017_noDuplicates.root"]

#outputDir     	= "/raid/raid7/kshi/Zprime/20190724/SkimTree_Run2016_signalregion_Data/"
#outputDir      = "/raid/raid7/kshi/Zprime/20190827/SkimTree_Run2016_MMM_Data/"
#outputDir     	= t2_prefix+"/store/user/t2/users/nikmenendez/skimmed/data2018/"

if year==2016:
	fileNames = [ #2016
	    "DoubleEG.root",
	    "DoubleMuon.root",
	    "MuonEG.root",
	    "SingleElectron.root",
	    "SingleMuon.root",
	    ]
elif year==2017:
	fileNames = [ #2017
	    'SingleElectron_Run2017-17Nov2017.root',
	    'SingleMuon_Run2017-17Nov2017.root',
	    'DoubleMuon_Run2017-17Nov2017.root',
	    'DoubleEG_Run2017-17Nov2017.root',
	    'MuonEG_Run2017-17Nov2017.root',
	    ]
elif year==2018:
	fileNames = [ #2018
	    'Data_Run2018A-17Sep2018_noDuplicates.root',
	    'Data_Run2018B-17Sep2018_noDuplicates.root',
	    'Data_Run2018C-17Sep2018_noDuplicates.root',
	    'Data_Run2018D-PromptReco-v2_noDuplicates.root',
	    ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteWto3lMMMTreeProducer_Data_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteWto3lMMMTreeProducer_Data(
            0.35,
            0.35,
            outputDir,
            fileName,
            True,
            )
    ana.setDebugMode(False)
    ana.loop(inputDir+fileName,inputTreeName)
#sendQuickMail(
#            ["kshi@cern.ch",],
#            "UFHZZLiteAnalyzer finished processing ("+getTimeStamp()+") ",
#            "\n".join([
#                "Input directory: "+inputDir,
#                "Output directory: "+outputDir,
#                "File: "+", ".join(fileNames),
#                ]
#                ),
#            )
