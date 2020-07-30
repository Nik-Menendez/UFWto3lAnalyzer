import ROOT,os
from PyUtils.UFTier2 import t2_prefix
from PyUtils.Shell import makedirs
from Emailer.Utils import sendQuickMail,getTimeStamp

# ____________________________________________________________________________________________________________________________________ ||
#bkgTreeDirT2_Feb21      = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2l_Feb21/"
#bkgTreeDirT2_Aug10      = t2_prefix+"/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_2lskim_Aug10/"
#sigTreeDirT2            = t2_prefix+"/store/user/t2/users/mhl/rootfiles_2017/"
#inputTreeName           = "Ana/passedEvents"
##outputDir               = "/raid/raid7/kshi/Zprime/20190729/SkimTree_Run2016_mem_controlregion_MC/"
#outputDir               = "/raid/raid7/kshi/Zprime/20190827/SkimTree_Run2016_MEM_MC/" 

bkgTreeDirT2            = t2_prefix+"/store/user/t2/users/nikmenendez/2017_MC_bkg/"
sigTreeDirT2            = t2_prefix+"/store/user/t2/users/mhl/rootfiles_2017/"
sigTreeDir              = "/home/kshi/Zprime/Zp_data_Ntuple/"
inputTreeName           = "Ana/passedEvents"
outputDir               = "/home/nikmenendez/Data/Wto3l/test/"

fileNames = [
        #bkgTreeDirT2_Feb21+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Feb21+"TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Feb21+"WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #bkgTreeDirT2_Aug10+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
        #sigTreeDirT2+"WmTo3munu_ZpM15_13TeV_MadGraph5_pythia8-v4_hmei-RunIISummer16MiniAODv2.root",
        #sigTreeDirT2+"WmTo3munu_ZpM45_13TeV_MadGraph5_pythia8-v4_hmei-RunIISummer16MiniAODv2.root",
        #sigTreeDirT2+"WpTo3munu_ZpM15_13TeV_MadGraph5_pythia8-v4_hmei-RunIISummer16MiniAODv2.root",
        #sigTreeDirT2+"WpTo3munu_ZpM45_13TeV_MadGraph5_pythia8-v4_hmei-RunIISummer16MiniAODv2.root",
        bkgTreeDirT2+"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root",
        #bkgTreeDirT2+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root",
        #bkgTreeDirT2+"TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8.root",
        #bkgTreeDirT2+"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root",
        ]

# ____________________________________________________________________________________________________________________________________ ||
ROOT.gSystem.Load("include/LiteWto3lMEMTreeProducer_h.so")

makedirs(outputDir)
for fileName in fileNames:
    ana = ROOT.LiteWto3lMEMTreeProducer(
            0.35,
            0.35,
            outputDir,
            os.path.basename(fileName),
            True,
            )
    ana.setDebugMode(False)
    ana.loop(fileName,inputTreeName)
#sendQuickMail(
#            ["kshi@cern.ch",],
#            "UFHZZLiteAnalyzer finished processing ("+getTimeStamp()+") ",
#            "\n".join([
#                "Input directory: "+", ".join([bkgTreeDirT2_Feb21,bkgTreeDirT2_Aug10,]),
#                "Output directory: "+outputDir,
#                "File: "+", ".join(fileNames),
#                ]
#                ),
#            )
