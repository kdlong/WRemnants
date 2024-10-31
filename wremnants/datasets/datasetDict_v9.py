import copy

from utilities import common

lumicsv = f"{common.data_dir}/bylsoutput.csv"
lumijson = (
    f"{common.data_dir}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
)

dataDictV9 = {
    "dataPostVFP": {
        "filepaths": [
            "{BASE_PATH}/SingleMuon/NanoV9Run2016FDataPostVFP_{NANO_PROD_TAG}",
            "{BASE_PATH}/SingleMuon/NanoV9Run2016GDataPostVFP_{NANO_PROD_TAG}",
            "{BASE_PATH}/SingleMuon/NanoV9Run2016HDataPostVFP_{NANO_PROD_TAG}",
        ],
        "group": "Data",
        "lumicsv": lumicsv,
        "lumijson": lumijson,
    },
    "ZmumuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.xsec_ZmmPostVFP,
        "group": "Zmumu",
    },
    "DYJetsToMuMuMass10to50PostVFP": {
        "filepaths": [
            "{BASE_PATH}/DYJetsToMuMu_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.xsec_ZmmMass10to50PostVFP,
        "group": "DYlowMass",
    },
    "ZtautauPostVFP": {
        "filepaths": [
            "{BASE_PATH}/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_H2ErratumFix_PDF_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        # At least one tau->e or mu decay, so everything that's not all other decays
        "xsec": common.xsec_ZmmPostVFP * common.Z_TAU_TO_LEP_RATIO,
        "group": "Ztautau",
    },
    "WplusmunuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WplusJetsToMuNu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.xsec_WpmunuPostVFP,
        "group": "Wmunu",
    },
    "WminusmunuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WminusJetsToMuNu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.xsec_WmmunuPostVFP,
        "group": "Wmunu",
    },
    "WplustaunuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WplusJetsToTauNu_TauToMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.BR_TAUToMU * common.xsec_WpmunuPostVFP,
        "group": "Wtaunu",
    },
    "WminustaunuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WminusJetsToTauNu_TauToMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": common.BR_TAUToMU * common.xsec_WmmunuPostVFP,
        "group": "Wtaunu",
    },
    "TTLeptonicPostVFP": {
        "filepaths": [
            "{BASE_PATH}/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 88.29,
        "group": "Top",
    },
    "TTSemileptonicPostVFP": {
        "filepaths": [
            "{BASE_PATH}/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 366.34,
        "group": "Top",
    },
    "SingleTschanLepDecaysPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 3.609,
        "group": "Top",
    },
    "SingleTtWAntitopPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 19.55,  # 35.85 * (1.0-((1-0.1086*3)*(1-0.1086*3))) = 19.5 pb
        "group": "Top",
    },
    "SingleTtWTopPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 19.55,
        "group": "Top",
    },
    "SingleTtchanAntitopPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 80.0,
        "group": "Top",
    },
    "SingleTtchanTopPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 134.2,
        "group": "Top",
    },
    # inclusive samples, keep for reference
    # 'WWPostVFP' : {
    #                 'filepaths' :
    #                 ["{BASE_PATH}/BKGV9/WW_TuneCP5_13TeV-pythia8"],
    #                 'xsec' : 118.7,
    #                 'group' : "Diboson",
    # },
    # 'WZPostVFP' : {
    #                 'filepaths' :
    #                 ["{BASE_PATH}/BKGV9/WZ_TuneCP5_13TeV-pythia8"],
    #                 'xsec' : 47.026760,  # to check, taken from WZTo1L1Nu2Q dividing by BR: 10.71/(3*0.1086)/(1-3*0.033658-0.2)
    #                 'group' : "Diboson",
    # },
    ##
    "WWTo2L2NuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 12.6,  # 118.7*0.1086*0.1086*9
        "group": "Diboson",
    },
    "WWTo1L1Nu2QPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 52.146,  # 118.7*[(3*0.1086)*(1-3*0.1086)]*2 (2 is because one W or the other can go to Q)
        "group": "Diboson",
    },
    "WZTo3LNuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 4.91,  # 4.42965*1.109, 1.109 is the NLO to NNLO kfactor, for this one would need to make sure about the NLO XS, depends a lot on the dilepton mass cut
        "group": "Diboson",
    },
    "WZTo2Q2LPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 5.4341,  # 4.9*1.109
        "group": "Diboson",
    },
    "WZTo1L1Nu2QPostVFP": {
        "filepaths": [
            "{BASE_PATH}/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 11.781,  # 10.71*1.10
        "group": "Diboson",
    },
    "ZZTo2L2NuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 0.60,
        "group": "Diboson",
    },
    "ZZTo2Q2LPostVFP": {
        "filepaths": [
            "{BASE_PATH}/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 5.1,
        "group": "Diboson",
    },
    "QCDmuEnrichPt15PostVFP": {
        "filepaths": [
            "{BASE_PATH}/QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}"
        ],
        "xsec": 238800,
        "group": "QCD",
    },
    "GGToLLMass5to50PostVFP": {
        "filepaths": [
            "{BASE_PATH}/GGToLL_M-5To50_TuneCP5_13TeV-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": 9.448e02,
        "group": "PhotonInduced",
    },
    "GGToLLPostVFP": {
        "filepaths": [
            "{BASE_PATH}/GGToLL_TuneCP5_13TeV-pythia8/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": 14.93,
        "group": "PhotonInduced",
    },
    "QGToDYQTo2LPostVFP": {
        "filepaths": [
            "{BASE_PATH}/QGToDYQTo2L_TuneCP5_13TeV-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": 1.373,
        "group": "PhotonInduced",
    },
    "QGToWQToLNuPostVFP": {
        "filepaths": [
            "{BASE_PATH}/QGToWQToLNu_TuneCP5_13TeV-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
        ],
        "xsec": 4.224e01,
        "xsec_up": 4.827e01,
        "xsec_dn": 3.588e01,
        "group": "PhotonInduced",
    },
}

# extended version with additional samples (but missing some pdf sets)
dataDictV9extended = copy.deepcopy(dataDictV9)

dataDictV9extended["ZmumuPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)

dataDictV9extended["ZtautauPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/DYJetsToTauTau_M-50_AtLeastOneEorMuDecay_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)

dataDictV9extended["WplusmunuPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)

dataDictV9extended["WminusmunuPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/WminusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)

dataDictV9extended["WplustaunuPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/WplusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)

dataDictV9extended["WminustaunuPostVFP"]["filepaths"].extend(
    [
        "{BASE_PATH}/WminusJetsToTauNu_TauToMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_{NANO_PROD_TAG}",
    ]
)
