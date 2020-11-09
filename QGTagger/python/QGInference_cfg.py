import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')

# TODO: put this option in cmsRun scripts
options.register('processMode', 
    default='JetLevel', 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process mode: JetLevel or EventLevel")
# Skip Events.
options.register('skipEvents',
    default=0,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
# Set doECALstitched to 1 to produce JetSeeds and JetFrames.
options.register('doECALstitched',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doECALstitched")
# Set doTracksAtECALstitchedPt to 1 to produce JetSeeds and JetFrames.
options.register('doTracksAtECALstitchedPt',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doTracksAtECALstitchedPt")
# Set doTracksAtECALadjPt to 1 to produce JetSeeds and JetFrames.
options.register('doTracksAtECALadjPt',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doTracksAtECALadjPt")
# Set doHBHEenergy to 1 to produce JetSeeds and JetFrames.
options.register('doHBHEenergy',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "set doHBHEenergy")
# Name of the QGInference model to be used for inference.
options.register('QGModelName',
    default='ResNet.pb',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "QGInference Model name")
options.parseArguments()

process = cms.Process("QGClassifier")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    #input = cms.untracked.int32(10)
    )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      #"file:myOutputFile.root"#SinglePhotonPt50_noPU_AODSIM.root
      )
    , skipEvents = cms.untracked.uint32(0)#options.skipEvents
    )
print (" >> Loaded",len(options.inputFiles),"input files from list.")

process.load("E2eDL.FrameProducers.DetFrameProducer_cfi")
process.load("E2eDL.FrameProducers.JetFrameProducer_cfi")
process.load("E2eDL.QGTagger.QGTagger_cfi")
#process.EGTagger.EGModelName = options.EGModelName
process.JetFrames.jetCollection = cms.string("ak4")
process.JetFrames.minJetPt = cms.double(35.)
process.JetFrames.maxJetEta = cms.double(2.4)
process.JetFrames.doHBHEenergy = options.doHBHEenergy
process.JetFrames.doECALstitched = options.doECALstitched
process.JetFrames.doTracksAtECALstitchedPt = options.doTracksAtECALstitchedPt
process.JetFrames.doTracksAtECALadjPt = options.doTracksAtECALadjPt

process.DetFrames.doHBHEenergy = options.doHBHEenergy
process.DetFrames.doECALstitched = options.doECALstitched
process.DetFrames.doTracksAtECALstitchedPt = options.doTracksAtECALstitchedPt
process.DetFrames.doTracksAtECALadjPt = options.doTracksAtECALadjPt

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('QGPt+QGFrames.root') 
    )
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root")#options.outputFile
    )

process.p = cms.Path(process.DetFrames + process.JetFrames+process.QGTagger)
process.ep=cms.EndPath(process.out)

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
