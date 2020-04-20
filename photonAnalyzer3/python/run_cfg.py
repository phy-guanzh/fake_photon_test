import FWCore.ParameterSet.Config as cms

process = cms.Process("photonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/z/zguan/EWZG/CMSSW_10_2_18/src/photonAna/photonAnalyzer/ZG1.root'
    )
)

process.demo = cms.EDAnalyzer('photonAnalyzer',
		 Photons =  cms.InputTag("slimmedPhotons"),
		 rho = cms.InputTag("fixedGridRhoFastjetAll"),
)


process.p = cms.Path(process.demo)

process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))
