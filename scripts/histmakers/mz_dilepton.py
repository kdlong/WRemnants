from utilities import boostHistHelpers as hh, common, logging, differential
from utilities.io_tools import output_tools

parser,initargs = common.common_parser(True)

import ROOT
import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation, muon_calibration, muon_selections, unfolding_tools
from wremnants.histmaker_tools import scale_to_data, aggregate_groups
from wremnants.datasets.dataset_tools import getDatasets
from wremnants.datasets.datagroups import Datagroups
import hist
import lz4.frame
import math
import time
import os

parser.add_argument("--csVarsHist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--axes", type=str, nargs="*", default=["mll", "ptll"], help="")
parser.add_argument("--finePtBinning", action='store_true', help="Use fine binning for ptll")
parser.add_argument("--useDileptonTriggerSelection", action='store_true', help="Use dilepton trigger selection (default uses the Wlike one, with one triggering muon and odd/even event selection to define its charge, staying agnostic to the other)")
parser.add_argument("--noAuxiliaryHistograms", action="store_true", help="Remove auxiliary histograms to save memory (removed by default with --unfolding or --theoryAgnostic)")

parser = common.set_parser_default(parser, "aggregateGroups", ["Diboson", "Top", "Wtaunu", "Wmunu"])
parser = common.set_parser_default(parser, "ewTheoryCorr", ["virtual_ew", "pythiaew_ISR", "horaceqedew_FSR", "horacelophotosmecoffew_FSR",])
parser = common.set_parser_default(parser, "excludeProcs", ["QCD"])

args = parser.parse_args()
isUnfolding = args.analysisMode == "unfolding"
isPoiAsNoi = isUnfolding and args.poiAsNoi

if isUnfolding:
    parser = common.set_parser_default(parser, "genAxes", ["ptVGen", "absYVGen"])
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)


thisAnalysis = ROOT.wrem.AnalysisType.Dilepton if args.useDileptonTriggerSelection else ROOT.wrem.AnalysisType.Wlike
era = args.era
datasets = getDatasets(maxFiles=args.maxFiles,
                       filt=args.filterProcs,
                       excl=args.excludeProcs, 
                       nanoVersion="v9",
                       base_path=args.dataPath,
                       extended = "msht20an3lo" not in args.pdfs,
                       era = era)

# dilepton invariant mass cuts
mass_min = 60
mass_max = 120

ewMassBins = theory_tools.make_ew_binning(mass = 91.1535, width = 2.4932, initialStep=0.010)

dilepton_ptV_binning = common.get_dilepton_ptV_binning(args.finePtBinning)
# available axes for dilepton validation plots
all_axes = {
    # "mll": hist.axis.Regular(60, 60., 120., name = "mll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "mll": hist.axis.Variable([60,70,75,78,80,82,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,100,102,105,110,120], name = "mll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),    "yll": hist.axis.Regular(20, -2.5, 2.5, name = "yll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "absYll": hist.axis.Regular(10, 0., 2.5, name = "absYll", underflow=False, overflow=not args.excludeFlow),
    "ptll": hist.axis.Variable(dilepton_ptV_binning, name = "ptll", underflow=False, overflow=not args.excludeFlow),
    "etaPlus": hist.axis.Variable([-2.4,-1.2,-0.3,0.3,1.2,2.4], name = "etaPlus"),
    "etaMinus": hist.axis.Variable([-2.4,-1.2,-0.3,0.3,1.2,2.4], name = "etaMinus"),
    "etaRegionSign": hist.axis.Regular(3, 0, 3, name = "etaRegionSign"),
    "etaRegionRange": hist.axis.Regular(3, 0, 3, name = "etaRegionRange"),
    "absEtaPlus": hist.axis.Regular(8, 0, 2.4, name = "absEtaPlus"),
    "absEtaMinus": hist.axis.Regular(8, 0, 2.4, name = "absEtaMinus"),
    "etaAbsEta": hist.axis.Variable([-2.4, -2.0, -1.6, -1.4, -1.2, -1.0, -0.6, 0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4], name = "etaAbsEta"),
    "etaSum": hist.axis.Regular(12, -4.8, 4.8, name = "etaSum"),
    "etaDiff": hist.axis.Variable([-4.8, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 4.8], name = "etaDiff"),
    "ptPlus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptPlus"),
    "ptMinus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptMinus"),
    "cosThetaStarll": hist.axis.Regular(20, -1., 1., name = "cosThetaStarll", underflow=False, overflow=False),
    "phiStarll": hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phiStarll"),
    #"charge": hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge") # categorical axes in python bindings always have an overflow bin, so use a regular
    "massVgen": hist.axis.Variable(ewMassBins, name = "massVgen", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMll": hist.axis.Variable(ewMassBins, name = "ewMll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMlly": hist.axis.Variable(ewMassBins, name = "ewMlly", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewLogDeltaM": hist.axis.Regular(100, -10, 4, name = "ewLogDeltaM", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "trigMuons_abseta0" : hist.axis.Regular(3, 0., 2.4, name = "trigMuons_abseta0", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "nonTrigMuons_eta0" : hist.axis.Regular(24, -2.4, 2.4, name = "nonTrigMuons_eta0", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "nonTrigMuons_pt0" : hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "nonTrigMuons_pt0"),
    "nonTrigMuons_charge0" : hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "nonTrigMuons_charge0"),
}

auxiliary_gen_axes = ["massVgen", # preFSR variables
    "ewMll", "ewMlly", "ewLogDeltaM" # ew variables
    ]

for a in args.axes:
    if a not in all_axes.keys():
        logger.error(f" {a} is not a known axes! Supported axes choices are {list(all_axes.keys())}")

nominal_cols = args.axes

if args.csVarsHist:
    nominal_cols += ["cosThetaStarll", "phiStarll"]

nominal_axes = [all_axes[a] for a in nominal_cols] 

if isUnfolding:
    unfolding_axes, unfolding_cols, unfolding_selections = differential.get_dilepton_axes(args.genAxes, common.get_gen_axes(isPoiAsNoi or inclusive, dilepton_ptV_binning), add_out_of_acceptance_axis=isPoiAsNoi)
    if not isPoiAsNoi:
        datasets = unfolding_tools.add_out_of_acceptance(datasets, group = "Zmumu")

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.theory_corrections.make_qcd_uncertainty_helper_by_helicity(is_w_like = True)

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # might never use it really anymore, but let's warn the user that this is obsolete
    logger.warning("Only SF with no uT dependence are implemented, and the treatment for trigger is like Wlike")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     is_w_like = True) 
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     what_analysis = thisAnalysis, isoEfficiencySmoothing=args.isoEfficiencySmoothing, smooth3D=args.smooth3dsf, isoDefinition=args.isolationDefinition)
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)
vertex_helper = wremnants.make_vertex_helper(era = era)

calib_filepaths = common.calib_filepaths
closure_filepaths = common.closure_filepaths
diff_weights_helper = ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths['tflite_file']) if (args.muonScaleVariation == 'smearingWeightsSplines' or args.validationHists) else None
mc_jpsi_crctn_helper, data_jpsi_crctn_helper, mc_jpsi_crctn_unc_helper, data_jpsi_crctn_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, calib_filepaths, make_uncertainty_helper=True)
z_non_closure_parametrized_helper, z_non_closure_binned_helper = muon_calibration.make_Z_non_closure_helpers(args, calib_filepaths, closure_filepaths)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)

closure_unc_helper = wremnants.muon_calibration.make_closure_uncertainty_helper(common.closure_filepaths["parametrized"])
closure_unc_helper_A = wremnants.muon_calibration.make_uniform_closure_uncertainty_helper(0, common.correlated_variation_base_size["A"])
closure_unc_helper_M = wremnants.muon_calibration.make_uniform_closure_uncertainty_helper(2, common.correlated_variation_base_size["M"])

smearing_helper, smearing_uncertainty_helper = (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()

bias_helper = muon_calibration.make_muon_bias_helpers(args) 

theory_corrs = [*args.theoryCorr, *args.ewTheoryCorr]
corr_helpers = theory_corrections.load_corr_helpers([d.name for d in datasets if d.name in common.vprocs], theory_corrs)

def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isWorZ = isW or isZ

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    df = df.Define("isEvenEvent", "event % 2 == 0")

    weightsum = df.SumAndCount("weight")

    axes = nominal_axes
    cols = nominal_cols

    if isUnfolding and dataset.name == "ZmumuPostVFP":
        fidmode = "mz_window" if inclusive else "mz_dilepton"
        df = unfolding_tools.define_gen_level(df, args.genLevel, dataset.name, mode=fidmode)
        fidargs = unfolding_tools.get_fiducial_args(fidmode, pt_min=args.pt[1], pt_max=args.pt[2])

        if hasattr(dataset, "out_of_acceptance"):
            logger.debug("Reject events in fiducial phase space")
            df = unfolding_tools.select_fiducial_space(df, mode=fidmode, selections=[], accept=False, **fidargs)
        else:
            logger.debug("Select events in fiducial phase space")
            df = unfolding_tools.select_fiducial_space(df, mode="mz_dilepton", selections=[], select=not isPoiAsNoi, accept=True, **fidargs)
            unfolding_tools.add_xnorm_histograms(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols)
            if not isPoiAsNoi:
                axes = [*nominal_axes, *unfolding_axes] 
                cols = [*nominal_cols, *unfolding_cols]

    if not args.noAuxiliaryHistograms and isZ:
        # gen level variables before selection
        for obs in auxiliary_gen_axes:
            df_gen = df
            df_gen = df_gen.DefinePerSample("exp_weight", "1.0")

            df_gen = theory_tools.define_theory_weights_and_corrs(df_gen, dataset.name, corr_helpers, args)

            results.append(df_gen.HistoBoost(f"gen_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            df_gen = syst_tools.add_theory_hists(results, df_gen, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"gen_{obs}", for_wmass=False)
    df = df.Filter(muon_selections.hlt_string(era))

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, args.pt[1], args.pt[2], dataset.group, nMuons=2, use_trackerMuons=args.trackerMuons, use_isolation=True, isoDefinition=args.isolationDefinition)

    df = muon_selections.define_trigger_muons(df, dilepton=args.useDileptonTriggerSelection)

    df = muon_selections.select_z_candidate(df, mass_min, mass_max)

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "trigMuons")
    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "nonTrigMuons")

    if args.useDileptonTriggerSelection:
        df = muon_selections.apply_triggermatching_muon(df, dataset, "trigMuons", "nonTrigMuons", era=era)
        df = df.Alias("muonsMinus_pt0", "trigMuons_pt0")
        df = df.Alias("muonsPlus_pt0", "nonTrigMuons_pt0")
        df = df.Alias("muonsMinus_eta0", "trigMuons_eta0")
        df = df.Alias("muonsPlus_eta0", "nonTrigMuons_eta0")
        df = df.Alias("muonsMinus_mom4", "trigMuons_mom4")
        df = df.Alias("muonsPlus_mom4", "nonTrigMuons_mom4")
    else:
        df = muon_selections.apply_triggermatching_muon(df, dataset, "trigMuons", era=era)
        df = df.Define("trigMuon_isNegative",  "trigMuons_charge0 == -1")
        df = df.Define("muonsMinus_pt0",  "trigMuon_isNegative ? trigMuons_pt0 : nonTrigMuons_pt0")
        df = df.Define("muonsPlus_pt0",   "trigMuon_isNegative ? nonTrigMuons_pt0 : trigMuons_pt0")
        df = df.Define("muonsMinus_eta0", "trigMuon_isNegative ? trigMuons_eta0 : nonTrigMuons_eta0")
        df = df.Define("muonsPlus_eta0",  "trigMuon_isNegative ? nonTrigMuons_eta0 : trigMuons_eta0")
        df = df.Define("muonsMinus_mom4", "trigMuon_isNegative ? trigMuons_mom4 : nonTrigMuons_mom4")
        df = df.Define("muonsPlus_mom4",  "trigMuon_isNegative ? nonTrigMuons_mom4 : trigMuons_mom4")
    
    df = df.Define("ptll", "ll_mom4.pt()")
    df = df.Define("yll", "ll_mom4.Rapidity()")
    df = df.Define("absYll", "std::fabs(yll)")
    # "renaming" to write out corresponding axis
    df = df.Alias("ptMinus", "muonsMinus_pt0")
    df = df.Alias("ptPlus", "muonsPlus_pt0")
    df = df.Alias("etaMinus", "muonsMinus_eta0")
    df = df.Alias("etaPlus", "muonsPlus_eta0")
    df = df.Define("absEtaMinus", "std::fabs(etaMinus)")
    df = df.Define("absEtaPlus", "std::fabs(etaPlus)")
    df = df.Define("etaAbsEta", "absEtaMinus > absEtaPlus ? etaMinus : etaPlus")

    df = df.Define("etaRegionRange", "(std::abs(muonsPlus_eta0) > 0.9) + (std::abs(muonsMinus_eta0) > 0.9)") # eta region: 0: barrel-barrel, 1: endcap-barrel, 2: endcap-endcap
    df = df.Define("etaRegionSign", "(muonsPlus_eta0 > 0) + (muonsMinus_eta0 > 0)") # eta region: 0: both muons in negative eta, 1: one muon in negative eta, 2: both muons in positive eta

    df = df.Define("etaSum", "muonsPlus_eta0 + muonsMinus_eta0") 
    df = df.Define("etaDiff", "muonsPlus_eta0 - muonsMinus_eta0") # plus - minus 

    df = df.Define("csSineCosThetaPhill", "wrem::csSineCosThetaPhi(muonsPlus_mom4, muonsMinus_mom4)")
    df = df.Define("cosThetaStarll", "csSineCosThetaPhill.costheta")
    df = df.Define("phiStarll", "std::atan2(csSineCosThetaPhill.sinphi, csSineCosThetaPhill.cosphi)")

    # TODO might need to add an explicit cut on trigMuons_pt0 in case nominal pt range
    # extends below 26 GeV e.g. for calibration test purposes
    df = df.Define("trigMuons_abseta0", "std::fabs(trigMuons_eta0)")

    logger.debug(f"Define weights and store nominal histograms")

    if dataset.is_data:
        results.append(df.HistoBoost("nominal", axes, cols))
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        df = df.Define("weight_vtx", vertex_helper, ["GenVtx_z", "Pileup_nTrueInt"])
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        if era == "2016PostVFP":
            weight_expr = "weight_pu*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom"
        else:
            weight_expr = "weight_pu*L1PreFiringWeight_Muon_Nom*L1PreFiringWeight_ECAL_Nom"

        if not args.noVertexWeight:
            weight_expr += "*weight_vtx"            

        muonVarsForSF = ["pt0", "eta0", "SApt0", "SAeta0", "uT0", "charge0"]
        if args.useDileptonTriggerSelection:
            muonVarsForSF.append("passTrigger0")            
        # careful, first all trig variables, then all nonTrig
        columnsForSF = [f"{t}Muons_{v}" for t in ["trig", "nonTrig"] for v in muonVarsForSF]

        df = muon_selections.define_muon_uT_variable(df, isWorZ, smooth3dsf=args.smooth3dsf, colNamePrefix="trigMuons")
        df = muon_selections.define_muon_uT_variable(df, isWorZ, smooth3dsf=args.smooth3dsf, colNamePrefix="nonTrigMuons")
        if not args.smooth3dsf:
            columnsForSF.remove("trigMuons_uT0")
            columnsForSF.remove("nonTrigMuons_uT0")

        if not args.noScaleFactors:
            # FIXME: add flags for pass_trigger for both leptons
            df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, columnsForSF)
            weight_expr += "*weight_fullMuonSF_withTrackingReco"

        logger.debug(f"Experimental weight defined: {weight_expr}")
        df = df.Define("exp_weight", weight_expr)
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)

        results.append(df.HistoBoost("weight", [hist.axis.Regular(100, -2, 2)], ["nominal_weight"], storage=hist.storage.Double()))
        results.append(df.HistoBoost("nominal", axes, [*cols, "nominal_weight"]))

    if isUnfolding and isPoiAsNoi and dataset.name == "ZmumuPostVFP":
        noiAsPoiHistName = Datagroups.histName("nominal", syst="yieldsUnfolding")
        logger.debug(f"Creating special histogram '{noiAsPoiHistName}' for unfolding to treat POIs as NOIs")
        results.append(df.HistoBoost(noiAsPoiHistName, [*nominal_axes, *unfolding_axes], [*nominal_cols, *unfolding_cols, "nominal_weight"]))       

    #for obs in ["ptll", "mll", "yll", "cosThetaStarll", "phiStarll", "etaPlus", "etaMinus", "ptPlus", "ptMinus"]:
    for obs in ["ptll", "mll", "yll", "etaPlus", "etaMinus", "ptPlus", "ptMinus"]:
        if dataset.is_data:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs]))
        else:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            if isWorZ:
                df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"nominal_{obs}", for_wmass=False)

    if not args.noAuxiliaryHistograms and isZ:
        # gen level variables
        for obs in auxiliary_gen_axes:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"nominal_{obs}", for_wmass=False)

    # test plots
    if args.validationHists and args.useDileptonTriggerSelection:
        df_plusTrig = df.Filter("trigMuons_passTrigger0")
        df_minusTrig = df.Filter("nonTrigMuons_passTrigger0")
        df_bothTrig = df.Filter("trigMuons_passTrigger0 && nonTrigMuons_passTrigger0")
        df_plusTrigOnly = df.Filter("trigMuons_passTrigger0 && !nonTrigMuons_passTrigger0")
        df_minusTrigOnly = df.Filter("nonTrigMuons_passTrigger0 && !trigMuons_passTrigger0")
        for obs in ["etaPlus", "etaMinus", "ptPlus", "ptMinus"]:
            if dataset.is_data:
                results.append(df_plusTrig.HistoBoost(f"nominal_{obs}_plusTrig", [all_axes[obs]], [obs]))
                results.append(df_minusTrig.HistoBoost(f"nominal_{obs}_minusTrig", [all_axes[obs]], [obs]))
                results.append(df_bothTrig.HistoBoost(f"nominal_{obs}_bothTrig", [all_axes[obs]], [obs]))
                results.append(df_plusTrigOnly.HistoBoost(f"nominal_{obs}_plusTrigOnly", [all_axes[obs]], [obs]))
                results.append(df_minusTrigOnly.HistoBoost(f"nominal_{obs}_minusTrigOnly", [all_axes[obs]], [obs]))
            else:
                results.append(df_plusTrig.HistoBoost(f"nominal_{obs}_plusTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_minusTrig.HistoBoost(f"nominal_{obs}_minusTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_bothTrig.HistoBoost(f"nominal_{obs}_bothTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_plusTrigOnly.HistoBoost(f"nominal_{obs}_plusTrigOnly", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_minusTrigOnly.HistoBoost(f"nominal_{obs}_minusTrigOnly", [all_axes[obs]], [obs, "nominal_weight"]))

    if not dataset.is_data and not args.onlyMainHistograms:

        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, axes, cols, what_analysis=thisAnalysis, smooth3D=args.smooth3dsf)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, axes, cols)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isWorZ:

            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, axes, cols, for_wmass=False)

            reco_sel = "vetoMuonsPre"
            require_prompt = "tau" not in dataset.name
            df = muon_calibration.define_genFiltered_recoMuonSel(df, reco_sel, require_prompt)
            reco_sel_GF = muon_calibration.getColName_genFiltered_recoMuonSel(reco_sel, require_prompt)
            df = muon_calibration.define_matched_gen_muons_kinematics(df, reco_sel_GF)
            df = muon_calibration.calculate_matched_gen_muon_kinematics(df, reco_sel_GF)
            df = muon_calibration.define_matched_reco_muon_kinematics(df, reco_sel_GF)

            ####################################################
            # nuisances from the muon momemtum scale calibration 
            if (args.muonCorrData in ["massfit", "lbl_massfit"]):
                input_kinematics = [
                    f"{reco_sel_GF}_recoPt",
                    f"{reco_sel_GF}_recoEta",
                    f"{reco_sel_GF}_recoCharge",
                    f"{reco_sel_GF}_genPt",
                    f"{reco_sel_GF}_genEta",
                    f"{reco_sel_GF}_genCharge"
                ]
                if diff_weights_helper:
                    df = df.Define(f'{reco_sel_GF}_response_weight', diff_weights_helper, [*input_kinematics])
                    input_kinematics.append(f'{reco_sel_GF}_response_weight')

                # muon scale variation from stats. uncertainty on the jpsi massfit
                df = df.Define(
                    "nominal_muonScaleSyst_responseWeights_tensor", data_jpsi_crctn_unc_helper,
                    [*input_kinematics, "nominal_weight"]
                )
                muonScaleSyst_responseWeights = df.HistoBoost(
                    "nominal_muonScaleSyst_responseWeights", axes,
                    [*cols, "nominal_muonScaleSyst_responseWeights_tensor"],
                    tensor_axes = data_jpsi_crctn_unc_helper.tensor_axes, storage=hist.storage.Double()
                )
                results.append(muonScaleSyst_responseWeights)

                df = muon_calibration.add_resolution_uncertainty(df, axes, results, cols, smearing_uncertainty_helper, reco_sel_GF)

                if args.nonClosureScheme in ["A-M-separated", "A-only"]:
                    # add the ad-hoc Z non-closure nuisances from the jpsi massfit to muon scale unc
                    df = df.DefinePerSample("AFlag", "0x01")
                    df = df.Define(
                        "Z_non_closure_parametrized_A", z_non_closure_parametrized_helper,
                        [*input_kinematics, "nominal_weight", "AFlag"]
                    )
                    hist_Z_non_closure_parametrized_A = df.HistoBoost(
                        "nominal_Z_non_closure_parametrized_A",
                        axes, [*cols, "Z_non_closure_parametrized_A"],
                        tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                        storage=hist.storage.Double()
                    )
                    results.append(hist_Z_non_closure_parametrized_A)

                if args.nonClosureScheme in ["A-M-separated", "binned-plus-M", "M-only"]:
                    df = df.DefinePerSample("MFlag", "0x04")
                    df = df.Define(
                        "Z_non_closure_parametrized_M", z_non_closure_parametrized_helper,
                        [*input_kinematics, "nominal_weight", "MFlag"]
                    )
                    hist_Z_non_closure_parametrized_M = df.HistoBoost(
                        "nominal_Z_non_closure_parametrized_M",
                        axes, [*cols, "Z_non_closure_parametrized_M"],
                        tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                        storage=hist.storage.Double()
                    )
                    results.append(hist_Z_non_closure_parametrized_M)

                if args.nonClosureScheme == "A-M-combined":
                    df = df.DefinePerSample("AMFlag", "0x01 | 0x04")
                    df = df.Define("Z_non_closure_parametrized", z_non_closure_parametrized_helper,
                        [
                            *input_kinematics,
                            "nominal_weight",
                            "AMFlag"
                        ]
                    )
                    hist_Z_non_closure_parametrized = df.HistoBoost(
                        "Z_non_closure_parametrized_gaus" if args.muonScaleVariation == 'smearingWeightsGaus' else "nominal_Z_non_closure_parametrized",
                        axes,
                        [*cols, "Z_non_closure_parametrized"],
                        tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                        storage=hist.storage.Double()
                    )
                    results.append(hist_Z_non_closure_parametrized)

                # extra uncertainties from non-closure stats
                df = df.Define("muonScaleClosSyst_responseWeights_tensor_splines", closure_unc_helper,
                    [*input_kinematics, "nominal_weight"]
                )
                nominal_muonScaleClosSyst_responseWeights = df.HistoBoost(
                    "nominal_muonScaleClosSyst_responseWeights", axes,
                    [*cols, "muonScaleClosSyst_responseWeights_tensor_splines"],
                    tensor_axes = closure_unc_helper.tensor_axes,
                    storage = hist.storage.Double()
                )
                results.append(nominal_muonScaleClosSyst_responseWeights)

                # extra uncertainties for A (fully correlated)
                df = df.Define("muonScaleClosASyst_responseWeights_tensor_splines", closure_unc_helper_A,
                    [*input_kinematics, "nominal_weight"]
                )
                nominal_muonScaleClosASyst_responseWeights = df.HistoBoost(
                    "nominal_muonScaleClosASyst_responseWeights", axes,
                    [*cols, "muonScaleClosASyst_responseWeights_tensor_splines"],
                    tensor_axes = closure_unc_helper_A.tensor_axes,
                    storage = hist.storage.Double()
                )
                results.append(nominal_muonScaleClosASyst_responseWeights)

                # extra uncertainties for M (fully correlated)
                df = df.Define("muonScaleClosMSyst_responseWeights_tensor_splines", closure_unc_helper_M,
                    [*input_kinematics, "nominal_weight"]
                )
                nominal_muonScaleClosMSyst_responseWeights = df.HistoBoost(
                    "nominal_muonScaleClosMSyst_responseWeights", axes,
                    [*cols, "muonScaleClosMSyst_responseWeights_tensor_splines"],
                    tensor_axes = closure_unc_helper_M.tensor_axes,
                    storage = hist.storage.Double()
                )
                results.append(nominal_muonScaleClosMSyst_responseWeights)

            ####################################################

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                syst_tools.add_muonscale_hist(
                    results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, axes, cols,
                    muon_eta="trigMuons_eta0") ## FIXME: what muon to choose ?

    if hasattr(dataset, "out_of_acceptance"):
        # Rename dataset to not overwrite the original one
        dataset.name = dataset.name+"OOA"
    
    return results, weightsum

logger.debug(f"Datasets are {[d.name for d in datasets]}")
resultdict = narf.build_and_run(datasets, build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args)
