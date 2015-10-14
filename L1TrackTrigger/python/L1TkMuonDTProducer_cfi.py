import FWCore.ParameterSet.Config as cms

L1TkMuonsDT = cms.EDProducer("L1TkMuonDTProducer")


L1TkMuonsMerge =  cms.EDProducer("L1TkMuonMerger",
   TkMuonCollections = cms.VInputTag( cms.InputTag("l1TkMuonsExt",""),
                                      cms.InputTag("l1TkMuonsExtCSC","") ),
   absEtaMin = cms.vdouble( 0. , 1.1),
   absEtaMax = cms.vdouble( 1.1 , 5.0)
)


        # --- or using the Padova muons in the central region:
L1TkMuonsMergeWithDT = cms.EDProducer("L1TkMuonMerger",
   TkMuonCollections = cms.VInputTag(
                                      cms.InputTag("L1TkMuonsDT","DTMatchInwardsTTTrackFullReso"),
                                      cms.InputTag("l1TkMuonsExt",""),
                                      cms.InputTag("l1TkMuonsExtCSC","") ),
   absEtaMin = cms.vdouble( 0. , 1.1, 1.25),
   absEtaMax = cms.vdouble( 1.1 , 1.25, 5.0)
)
