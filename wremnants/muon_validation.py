import hist

import narf
from utilities import boostHistHelpers as hh
from utilities import common, logging

narf.clingutils.Declare('#include "muon_validation.hpp"')

logger = logging.child_logger(__name__)


# "muon" is for mw; "muons" is for wlike, for which we select one of the trig/nonTrig muons
def define_cvh_muon_kinematics(df):
    df = df.Define("goodMuons_cvh_pt0", "Muon_cvhPt[gooodMuons][0]")
    df = df.Define("goodMuons_cvh_eta0", "Muon_cvhEta[goodMuons][0]")
    df = df.Define("goodMuons_cvh_phi0", "Muon_cvhPhi[goodMuons][0]")
    return df


def define_cvh_muons_kinematics(df):
    df = df.Define("trigMuons_cvh_pt0", "Muon_correctedPt[trigMuons][0]")
    df = df.Define("trigMuons_cvh_eta0", "Muon_correctedEta[trigMuons][0]")
    df = df.Define("trigMuons_cvh_phi0", "Muon_correctedPhi[trigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_pt0", "Muon_correctedPt[nonTrigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_eta0", "Muon_correctedEta[nonTrigMuons][0]")
    df = df.Define("nonTrigMuons_cvh_phi0", "Muon_correctedPhi[nonTrigMuons][0]")
    return df


def define_jpsi_crctd_muons_pt_unc(df, helper):
    df = df.Define(
        "trigMuons_jpsi_crctd_pt_unc",
        helper,
        [
            "trigMuons_cvh_eta",
            "trigMuons_cvh_pt",
            "trigMuons_charge",
            "trigMuons_jpsi_crctd_pt",
        ],
    )
    df = df.Define(
        "nonTrigMuons_jpsi_crctd_pt_unc",
        helper,
        [
            "nonTrigMuons_cvh_eta",
            "nonTrigMuons_cvh_pt",
            "nonTrigMuons_charge",
            "nonTrigMuons_jpsi_crctd_pt",
        ],
    )
    return df


def define_jpsi_crctd_z_mass(df):
    df = df.Define(
        "trigMuons_jpsi_crctd_mom4",
        (
            "ROOT::Math::PtEtaPhiMVector("
            "trigMuons_jpsi_crctd_pt, trigMuons_cvh_eta, trigMuons_cvh_phi, wrem::muon_mass)"
        ),
    )
    df = df.Define(
        "nonTrigMuons_jpsi_crctd_mom4",
        (
            "ROOT::Math::PtEtaPhiMVector("
            "nonTrigMuons_jpsi_crctd_pt, nonTrigMuons_cvh_eta, nonTrigMuons_cvh_phi, wrem::muon_mass)"
        ),
    )
    df = df.Define(
        "Z_jpsi_crctd_mom4",
        "ROOT::Math::PxPyPzEVector(trigMuons_jpsi_crctd_mom4)+ROOT::Math::PxPyPzEVector(nonTrigMuons_jpsi_crctd_mom4)",
    )
    df = df.Define("massZ_jpsi_crctd", "Z_jpsi_crctd_mom4.mass()")
    return df


def define_jpsi_crctd_unc_z_mass(df):
    df = df.Define(
        "trigMuons_jpsi_crctd_mom4_unc",
        (
            "ROOT::VecOps::RVec<double> res(trigMuons_jpsi_crctd_pt_unc.size());"
            "for (int i = 0; i < trigMuons_jpsi_crctd_pt_unc.size(); i++) {"
            "    res[i] = ("
            "       ROOT::Math::PtEtaPhiMVector("
            "           trigMuons_jpsi_crctd_pt_unc[i],"
            "           trigMuons_cvh_eta,"
            "           trigMuons_cvh_phi,"
            "           wrem::muon_mass"
            "       )"
            "    );"
            "}"
            "return res;"
        ),
    )
    df = df.Define(
        "nonTrigMuons_jpsi_crctd_mom4_unc",
        (
            "ROOT::VecOps::RVec<double> res(nonTrigMuons_jpsi_crctd_pt_unc.size());"
            "for (int i = 0; i < nonTrigMuons_jpsi_crctd_pt_unc.size(); i++) {"
            "    res[i] = ("
            "        ROOT::Math::PtEtaPhiMVector("
            "            nonTrigMuons_jpsi_crctd_pt_unc[i],"
            "            nonTrigMuons_cvh_eta,"
            "            nonTrigMuons_cvh_phi,"
            "            wrem::muon_mass"
            "        )"
            "    );"
            "}"
            "return res;"
        ),
    )
    df = df.Define(
        "Z_jpsi_crctd_mom4_unc",
        (
            "ROOT::VecOps::RVec<double> res(trigMuons_jpsi_crctd_mom4_unc.size());"
            "for (int i = 0; i < trigMuons_jpsi_crctd_mom4_unc.size(); i++) {"
            "    res[i] = ("
            "        ROOT::Math::PxPyPzEVector(trigMuons_jpsi_crctd_mom4_unc[i]) +"
            "        ROOT::Math::PxPyPzEVector(nonTrigMuons_jpsi_crctd_mom4_unc[i])"
            "    )"
            "}"
            "return res"
        ),
    )
    df = df.Define(
        "massZ_jpsi_crctd_unc",
        (
            "ROOT::VecOps::RVec<double> res(Z_jpsi_crctd_mom4_unc.size());"
            "for (int i = 0; i < Z_jpsi_crctd_mom4_unc.size(); i++) {"
            "    res[i] = Z_jpsi_crctd_mom4_unc[i].mass()"
            "}"
            "return res"
        ),
    )
    return df


def define_cvh_reco_muon_kinematics(df, kinematic_vars=["pt", "eta", "phi", "charge"]):
    for var in kinematic_vars:
        col_name = f"goodMuons_{var.lower()}0_cvh"
        if not (col_name in df.GetColumnNames()):
            df = df.Define(col_name, f"Muon_cvh{var.capitalize()}[goodMuons][0]")
    return df


def define_uncrct_reco_muon_kinematics(
    df, kinematic_vars=["pt", "eta", "phi", "charge"]
):
    for var in kinematic_vars:
        col_name = f"goodMuons_{var.lower()}0_uncrct"
        if not (col_name in df.GetColumnNames()):
            df = df.Define(col_name, f"Muon_{var.lower()}[goodMuons][0]")
    return df


def define_reco_over_gen_cols(df, reco_type, kinematic_vars=["pt", "eta"]):
    kinematic_vars = common.string_to_list(kinematic_vars)
    df = define_cvh_reco_muon_kinematics(df, kinematic_vars)
    df = define_uncrct_reco_muon_kinematics(df, kinematic_vars)
    for var in kinematic_vars:
        reco_col = (
            f"goodMuons_{var.lower()}0"
            if reco_type == "crctd"
            else f"goodMuons_{var.lower()}0_{reco_type}"
        )
        df = df.Define(
            f"goodMuons_{var.lower()}0_{reco_type}_over_gen",
            f"{reco_col}/goodMuons_{var.lower()}0_gen",
        )
    return df


def make_reco_over_gen_hists(df, results):
    nominal_cols_crctd_over_gen = ["goodMuons_pt0_crctd_over_gen"]
    nominal_cols_cvh_over_gen = ["goodMuons_pt0_cvh_over_gen"]
    nominal_cols_uncrct_over_gen = ["goodMuons_pt0_uncrct_over_gen"]
    nominal_cols_gen_smeared_over_gen = [
        "goodMuons_pt0_gen_smeared_over_gen",
    ]
    axis_pt_reco_over_gen = hist.axis.Regular(
        1000, 0.9, 1.1, underflow=True, overflow=True, name="reco_pt_over_gen"
    )
    crctd_over_gen = df.HistoBoost(
        "crctd_over_gen",
        [axis_pt_reco_over_gen],
        [*nominal_cols_crctd_over_gen, "nominal_weight"],
        storage=hist.storage.Double(),
    )
    cvh_over_gen = df.HistoBoost(
        "cvh_over_gen",
        [axis_pt_reco_over_gen],
        [*nominal_cols_cvh_over_gen, "nominal_weight"],
        storage=hist.storage.Double(),
    )
    uncrct_over_gen = df.HistoBoost(
        "uncrct_over_gen",
        [axis_pt_reco_over_gen],
        [*nominal_cols_uncrct_over_gen, "nominal_weight"],
        storage=hist.storage.Double(),
    )
    gen_smeared_over_gen = df.HistoBoost(
        "gen_smeared_over_gen",
        [axis_pt_reco_over_gen],
        [*nominal_cols_gen_smeared_over_gen, "nominal_weight"],
        storage=hist.storage.Double(),
    )

    results.append(crctd_over_gen)
    results.append(cvh_over_gen)
    results.append(uncrct_over_gen)
    results.append(gen_smeared_over_gen)


def define_cols_for_manual_shifts(df):
    df = df.Define(
        "goodMuons_pt0_gen_smeared_scaleUp_mil", "goodMuons_pt0_gen_smeared * 1.001"
    )
    df = df.Define(
        "goodMuons_pt0_gen_smeared_scaleDn_mil", "goodMuons_pt0_gen_smeared / 1.001"
    )
    df = df.Define("goodMuons_pt0_scaleUp_tenthmil", "goodMuons_pt0 * 1.0001")
    df = df.Define("goodMuons_pt0_scaleDn_tenthmil", "goodMuons_pt0 / 1.0001")
    df = df.Define(
        "goodMuons_pt0_gen_smeared_scaleUp_tenthmil",
        "goodMuons_pt0_gen_smeared * 1.0001",
    )
    df = df.Define(
        "goodMuons_pt0_gen_smeared_scaleDn_tenthmil",
        "goodMuons_pt0_gen_smeared / 1.0001",
    )
    return df


def make_hists_for_event_weights_perse(
    df, nominal_axes, nominal_cols, weights_col, nominal_weight_col, method, results
):
    df = df.Define(f"weights_{method}_dn", f"{weights_col}(0,0)/{nominal_weight_col}")
    df = df.Define(f"weights_{method}_up", f"{weights_col}(0,1)/{nominal_weight_col}")
    axis_weights = hist.axis.Regular(
        1000, 0.99, 1.01, underflow=True, overflow=True, name="weights"
    )
    weights_dn = df.HistoBoost(
        f"weights_{method}_dn",
        [*nominal_axes, axis_weights],
        [*nominal_cols, f"weights_{method}_dn"],
        storage=hist.storage.Double(),
    )
    weights_up = df.HistoBoost(
        f"weights_{method}_up",
        [*nominal_axes, axis_weights],
        [*nominal_cols, f"weights_{method}_up"],
        storage=hist.storage.Double(),
    )
    results.append(weights_dn)
    results.append(weights_up)
    return df


# make the histograms for muon momentum scale shifts by
# manually shift the muon pt by a certain proportion
def make_hists_for_manual_scale_shifts(df, axes, cols, cols_gen_smeared, results):
    muonScaleVariationUpMil = df.HistoBoost(
        "nominal_muonScaleVariationUpMil",
        axes,
        [
            cols_gen_smeared[0],
            "goodMuons_pt0_gen_smeared_scaleUp_mil",
            *cols_gen_smeared[2:],
            "nominal_weight",
        ],
        storage=hist.storage.Double(),
    )
    muonScaleVariationDnMil = df.HistoBoost(
        "nominal_muonScaleVariationDnMil",
        axes,
        [
            cols_gen_smeared[0],
            "goodMuons_pt0_gen_smeared_scaleDn_mil",
            *cols_gen_smeared[2:],
            "nominal_weight",
        ],
        storage=hist.storage.Double(),
    )
    muonScaleVariationUpTenthmil = df.HistoBoost(
        "nominal_muonScaleVariationUpTenthmil",
        axes,
        [cols[0], "goodMuons_pt0_scaleUp_tenthmil", *cols[2:], "nominal_weight"],
        storage=hist.storage.Double(),
    )
    muonScaleVariationDnTenthmil = df.HistoBoost(
        "nominal_muonScaleVariationDnTenthmil",
        axes,
        [cols[0], "goodMuons_pt0_scaleDn_tenthmil", *cols[2:], "nominal_weight"],
        storage=hist.storage.Double(),
    )
    muonScaleVariationUpTenthmil_gen_smear = df.HistoBoost(
        "nominal_muonScaleVariationUpTenthmil_gen_smear",
        axes,
        [
            cols_gen_smeared[0],
            "goodMuons_pt0_gen_smeared_scaleUp_tenthmil",
            *cols_gen_smeared[2:],
            "nominal_weight",
        ],
        storage=hist.storage.Double(),
    )
    muonScaleVariationDnTenthmil_gen_smear = df.HistoBoost(
        "nominal_muonScaleVariationDnTenthmil_gen_smear",
        axes,
        [
            cols_gen_smeared[0],
            "goodMuons_pt0_gen_smeared_scaleDn_tenthmil",
            *cols_gen_smeared[2:],
            "nominal_weight",
        ],
        storage=hist.storage.Double(),
    )
    results.append(muonScaleVariationUpMil)
    results.append(muonScaleVariationDnMil)
    results.append(muonScaleVariationUpTenthmil)
    results.append(muonScaleVariationDnTenthmil)
    results.append(muonScaleVariationUpTenthmil_gen_smear)
    results.append(muonScaleVariationDnTenthmil_gen_smear)


def muon_scale_variation_from_manual_shift(
    resultdict,
    procs=["WplusmunuPostVFP", "WminusmunuPostVFP", "ZmumuPostVFP"],
):
    for proc in procs:
        proc_hists = resultdict[proc]["output"]
        manual_shift_hists = [
            proc_hists["nominal_muonScaleVariationDnTenthmil"].get(),
            proc_hists["nominal_muonScaleVariationUpTenthmil"].get(),
        ]
        proc_hists["muonScaleSyst_manualShift"] = narf.ioutils.H5PickleProxy(
            hh.combineUpDownVarHists(*manual_shift_hists)
        )


def make_hists_for_muon_scale_var_weights(
    df, axes, results, cols, nominal_cols_gen_smeared
):
    df = make_hists_for_event_weights_perse(
        df,
        axes,
        cols,
        "muonScaleSyst_responseWeights_tensor_gaus",
        "nominal_weight",
        "gaus",
        results,
    )
    df = make_hists_for_event_weights_perse(
        df,
        axes,
        cols,
        "muonScaleSyst_responseWeights_tensor_splines",
        "nominal_weight",
        "splines",
        results,
    )
    df = make_hists_for_event_weights_perse(
        df,
        axes,
        cols,
        "muonScaleSyst_responseWeights_tensor_massWeights",
        "nominal_weight",
        "massweights",
        results,
    )
    df = define_cols_for_manual_shifts(df)
    make_hists_for_manual_scale_shifts(
        df, axes, cols, nominal_cols_gen_smeared, results
    )
    return df
