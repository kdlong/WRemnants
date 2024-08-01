import hist
import numpy as np
from utilities import boostHistHelpers as hh
from utilities import common, logging
from utilities.fnnls import fnnls
from scipy.optimize import nnls
from scipy.special import comb
import pdb

from scipy.optimize import curve_fit
import uncertainties as unc
from uncertainties import unumpy as unp

logger = logging.child_logger(__name__)

# thresholds, 
abcd_thresholds={
    "pt":[26,28,30],
    "mt":[0,20,40],
    "iso":[0,4,8,12],
    "relIso":[0,0.15,0.3,0.45],
    "relJetLeptonDiff": [0,0.2,0.35,0.5],
    "dxy":[0,0.01,0.02,0.03]
}

def get_selection_edges(axis_name, upper_bound=False):
    # returns edges from pass to fail regions [x0, x1, x2, x3] e.g. x=[x0,x1], dx=[x1,x2], d2x=[x2,x3]
    if axis_name in abcd_thresholds:
        ts = abcd_thresholds[axis_name]
        if axis_name in ["mt", "pt"]:
            # low: failing, high: passing, no upper bound 
            return None, complex(0,ts[2]), complex(0,ts[1]), complex(0,ts[0])
        if axis_name in ["dxy", "iso", "relIso", "relJetLeptonDiff"]:
            # low: passing, high: failing, no upper bound 
            return complex(0,ts[0]), complex(0,ts[1]), complex(0,ts[2]), (complex(0,ts[3]) if upper_bound else None)
    else:
        raise RuntimeError(f"Can not find threshold for abcd axis {axis_name}")


# default abcd_variables to look for
abcd_variables = (("mt", "passMT"), ("relIso", "passIso"), ("iso", "passIso"), ("relJetLeptonDiff", "passIso"), ("dxy", "passDxy"))

def get_parameter_eigenvectors(params, cov, sign=1, force_positive=False):
    # diagonalize and get eigenvalues and eigenvectors
    e, v = np.linalg.eigh(cov) # The column eigenvectors[:, i] is the normalized eigenvector corresponding to the eigenvalue eigenvalues[i], 
    vT = np.transpose(v, axes=(*np.arange(v.ndim-2), v.ndim-1, v.ndim-2)) # transpose v to have row eigenvectors[i, :] for easier computations
    mag = np.sqrt(e)[...,np.newaxis] * vT * sign
    # duplicate params vector to have duplicated rows
    params_brd = np.broadcast_to(params[...,np.newaxis], mag.shape)
    paramsT = np.transpose(params_brd, axes=(*np.arange(params_brd.ndim-2), params_brd.ndim-1, params_brd.ndim-2))
    mask = (paramsT + mag < 0).any(axis=-1)
    if force_positive and mask.sum() > 0:
        logger.info(f"Force {mask.sum()} eigenvectors to be positive")
        # scaling the magnitude of the eigenvector to avoid negative coefficients
        min_idxs = paramsT + mag == np.min(paramsT + mag, axis=-1)[...,np.newaxis] # find indices with minimum value of each eigenvector
        min_paramsT = paramsT[min_idxs]
        min_mag = mag[min_idxs]
        factors = np.reshape(-1*min_paramsT/min_mag, mask.shape)
        factors = np.where(mask, factors, np.ones_like(mask))
        factors = np.broadcast_to(factors[...,np.newaxis], mag.shape)
        mag = factors * mag
    if np.sum(mag.sum(axis=-1) == 0):
        logger.warning(f"Found {np.sum(mag.sum(axis=-1) == 0)} eigenvector shifts where all coefficients are 0") 
    params_var = paramsT + mag
    if np.sum(params_var.sum(axis=-1) == 0):
        logger.warning(f"Found {np.sum(params_var.sum(axis=-1) == 0)} variations where all coefficients are 0") 
    return params_var

def make_eigenvector_predictons(params, cov, func, x1, x2=None, force_positive=False):
    # return alternate values i.e. nominal+/-variation
    params_up = get_parameter_eigenvectors(params, cov, sign=1, force_positive=force_positive)
    y_pred_up = func(x1, params_up) if x2 is None else func(x1, x2, params_up)
    y_pred_up = np.moveaxis(y_pred_up, params.ndim-1, -1) # put parameter variations last
    params_dn = get_parameter_eigenvectors(params, cov, sign=-1, force_positive=force_positive)
    y_pred_dn = func(x1, params_dn) if x2 is None else func(x1, x2, params_dn)
    y_pred_dn = np.moveaxis(y_pred_dn, params.ndim-1, -1) # put parameter variations last
    return np.stack((y_pred_up, y_pred_dn), axis=-1)

def poly(pol, order, order2=None):
    if order2 is None:
        if pol=="power":
            return lambda x, n, p=1: p * x**n
        elif pol=="bernstein":
            return lambda x, n, p=1, o=order: p * comb(o, n) * x**n * (1 - x)**(o - n)
    else:
        return lambda x, x2, n, m, p=1, o1=order, o2=order2: p * poly(pol, o1)(x, n) * poly(pol, o2)(x2, m)

def get_parameter_matrices(x, y, w, order, pol="power"):
    if x.shape != y.shape:
        x = np.broadcast_to(x, y.shape)
    stackX=[] # parameter matrix X 
    stackXTY=[] # and X.T @ Y
    f = poly(pol, order)
    for n in range(order+1):
        p = f(x, n)
        stackX.append( w * p )
        stackXTY.append((w**2 * p * y).sum(axis=(-1)))
    X = np.stack(stackX, axis=-1)
    XTY = np.stack(stackXTY, axis=-1)
    return X, XTY

def get_parameter_matrices_from2D(x, x2, y, w, order, order2=None, pol="power", flatten=False):
    x, x2 = np.broadcast_arrays(x[np.newaxis,...], x2[..., np.newaxis])
    x = np.broadcast_to(x, y.shape)
    x2 = np.broadcast_to(x2, y.shape)
    if order2 is None:
        order2 = [0,]*(order+1)
    elif type(order2) == int:
        order2 = [order2,]*(order+1)
    elif type(order2) == list: 
        if len(order2) < order+1:
            order2.append([0,]*(len(order2)-order+1))
    else:
        raise RuntimeError(f"Input 'order2' requires type 'None', 'int' or 'list'")
    stackX=[]   # parameter matrix X 
    stackXTY=[] # and X.T @ Y
    for n in range(order+1):
        f = poly(pol, order, order2[n])
        for m in range(order2[n]+1):
            p = f(x, x2, n, m)
            stackX.append(w * p)
            stackXTY.append((w**2 * p * y).sum(axis=(-2, -1)))
    X = np.stack(stackX, axis=-1)
    XTY = np.stack(stackXTY, axis=-1)
    if flatten:
        # flatten the 2D array into 1D
        newshape = (*y.shape[:-2],np.product(y.shape[-2:]))
        y = y.reshape(newshape)
        X = X.reshape(*newshape, X.shape[-1])

    return X, y, XTY

def get_regression_function(order1, order2=None, pol="power"):
    if order2 is None:
        def fsum(x, ps, o=order1):
            f = poly(pol, o)
            if hasattr(ps, "ndim") and ps.ndim > 1:
                return sum([f(x, n, ps[...,n,np.newaxis]) for n in range(o+1)])
            else:
                return sum([f(x, n, ps[n]) for n in range(o+1)])
    else:
        def fsum(x1, x2, ps, o1=order1, o2=order2):
            idx=0
            psum = 0
            x1, x2 = np.broadcast_arrays(x1[np.newaxis,...], x2[..., np.newaxis])
            if hasattr(ps, "ndim") and ps.ndim > 1:
                x1 = np.broadcast_to(x1, [*ps.shape[:-1], *x1.shape])
                x2 = np.broadcast_to(x2, [*ps.shape[:-1], *x2.shape]) 
                for n in range(o1+1):
                    f = poly(pol, o1, o2[n])
                    for m in range(o2[n]+1):
                        psum += f(x1, x2, n, m, ps[...,idx,np.newaxis,np.newaxis])
                        idx += 1
            else:
                for n in range(o1+1):
                    f = poly(pol, o1, o2[n])
                    for m in range(o2[n]+1):
                        psum += f(x1, x2, n, m, ps[idx])
                        idx += 1            
            return psum
    return fsum

def solve_leastsquare(X, XTY):
    # compute the transpose of X for the mt and parameter axes 
    XT = np.transpose(X, axes=(*np.arange(X.ndim-2), X.ndim-1, X.ndim-2))
    XTX = XT @ X
    # compute the inverse of the matrix in each bin (reshape to make last two axes contiguous, reshape back after inversion), 
    # this term is also the covariance matrix for the parameters
    XTXinv = np.linalg.inv(XTX.reshape(-1,*XTX.shape[-2:]))
    XTXinv = XTXinv.reshape((*XT.shape[:-2],*XTXinv.shape[-2:])) 
    params = np.einsum('...ij,...j->...i', XTXinv, XTY)
    return params, XTXinv

def solve_nonnegative_leastsquare(X, XTY):
    XT = np.transpose(X, axes=(*np.arange(X.ndim-2), X.ndim-1, X.ndim-2))
    XTX = XT @ X
    XTXinv = np.linalg.inv(XTX.reshape(-1,*XTX.shape[-2:]))
    XTXinv = XTXinv.reshape((*XT.shape[:-2],*XTXinv.shape[-2:])) 
    orig_shape = XTY.shape
    nBins = np.prod(orig_shape[:-1])
    XTY_flat = XTY.reshape(nBins, XTY.shape[-1])
    XTX_flat = XTX.reshape(nBins, XTX.shape[-2], XTX.shape[-1])
    # params = [fnnls(xtx, xty) for xtx, xty in zip(XTX_flat, XTY_flat)] # use fast nnls
    params = [nnls(xtx, xty)[0] for xtx, xty in zip(XTX_flat, XTY_flat)] # use scipy implementation of nnls (may be a bit slower)
    params = np.reshape(params, orig_shape)
    return params, XTXinv

def compute_chi2(y, y_pred, w=None, nparams=1):
    # goodness of fit from parameter matrix 'X', values 'y' and parateters 'params'
    residuals = (y - y_pred)
    if w is not None:
        residuals *= w
    chi2 = np.sum((residuals**2), axis=-1) # per fit chi2
    chi2_total = np.sum((residuals**2)) # chi2 of all fits together

    # Degrees of freedom calculation
    ndf = y.shape[-1] - nparams
    ndf_total = y.size - chi2.size*nparams

    from scipy import stats

    logger.info(f"Total chi2/ndf = {chi2_total}/{ndf_total} = {chi2_total/ndf_total} (p = {stats.chi2.sf(chi2_total, ndf_total)})")
    return chi2, ndf    

def extend_edges(traits, x):
    # extend array for underflow/overflow with distance from difference of two closest values 
    if traits.underflow:
        new_x = x[0] - x[1] + x[0]
        x = np.array([new_x, *x])
        logger.debug(f"Extend array by underflow bin using {new_x}")
    if traits.overflow:
        new_x = x[-1] + x[-1] - x[-2]
        x = np.array([*x, new_x])
        logger.debug(f"Extend array by overflow bin using {new_x}")
    return x

def get_rebinning(edges, axis_name):
    if axis_name == "mt":
        target_edges = common.get_binning_fakes_mt(high_mt_bins=True)
    elif axis_name == "pt":
        target_edges = common.get_binning_fakes_pt(edges[0],edges[-1])
    else:
        raise RuntimeError(f"No automatic rebinning known for axis {axis_name}")

    if len(edges) == len(target_edges) and all(edges == target_edges):
        logger.debug(f"Axis has already in target edges, no rebinning needed")
        return None

    i=0
    rebin = []
    for e in edges:
        if e >= target_edges[i]:
            rebin.append(e)
            i+=1
        if i >= len(target_edges):
            break

    if len(rebin) <= 2:
        logger.warning(f"No automatic rebinning possible for axis {axis_name}")
        return None
    
    logger.debug(f"Rebin axis {axis_name} to {rebin}")
    return rebin

def divide_arrays(num, den, cutoff=1):
    r = num/den
    criteria = abs(den) < cutoff
    if np.sum(criteria) > 0:
        logger.warning(f"Found {np.sum(criteria)} values in denominator below {cutoff}, the ratio will be set to 0 for those")
    r[abs(den) < cutoff] = 0 # if denumerator is close to 0 set ratio to zero to avoid large negative/positive values
    return r

class HistselectorABCD(object):
    def __init__(self, h, name_x=None, name_y=None,
        fakerate_axes=["eta","pt","charge"], 
        smoothing_axis_name="pt", 
        rebin_smoothing_axis="automatic", # can be a list of bin edges, "automatic", or None
        upper_bound_y=None, # using an upper bound on the abcd y-axis (e.g. isolation)
        integrate_x=True, # integrate the abcd x-axis in final histogram (allows simplified procedure e.g. for extrapolation method)   
    ):           
        self.upper_bound_y=upper_bound_y
        self.integrate_x = integrate_x

        self.name_x = name_x
        self.name_y = name_y
        if name_x is None or name_y is None:
            self.set_abcd_axes(h)

        self.axis_x = h.axes[self.name_x]
        self.axis_y = h.axes[self.name_y]

        self.sel_x = None
        self.sel_dx = None
        self.sel_y = None
        self.sel_dy = None
        self.set_selections_x()
        self.set_selections_y()

        if fakerate_axes is not None:
            self.fakerate_axes = fakerate_axes
            self.fakerate_integration_axes = [n for n in h.axes.name if n not in [self.name_x, self.name_y, *fakerate_axes]]
            logger.debug(f"Setting fakerate integration axes to {self.fakerate_integration_axes}")

        self.smoothing_axis_name = smoothing_axis_name
        edges = h.axes[smoothing_axis_name].edges
        if rebin_smoothing_axis == "automatic":
            self.rebin_smoothing_axis = get_rebinning(edges, self.smoothing_axis_name)
        else:
            self.rebin_smoothing_axis = rebin_smoothing_axis

        edges = edges if self.rebin_smoothing_axis is None else self.rebin_smoothing_axis
        edges = extend_edges(h.axes[self.smoothing_axis_name].traits, edges)
        self.smoothing_axis_min = edges[0]
        self.smoothing_axis_max = edges[-1]

    # A
    def get_hist_failX_failY(self, h):
        return h[{self.name_x: self.sel_dx, self.name_y: self.sel_dy}]
    # B
    def get_hist_failX_passY(self, h):
        return h[{self.name_x: self.sel_dx, self.name_y: self.sel_y}]
    # C
    def get_hist_passX_failY(self, h):
        return h[{self.name_x: self.sel_x, self.name_y: self.sel_dy}]
    # D
    def get_hist_passX_passY(self, h):
        return h[{self.name_x: self.sel_x, self.name_y: self.sel_y}]

    def set_abcd_axes(self, h):
        for a, b in abcd_variables:
            if self.name_x is None and a in h.axes.name: 
                self.name_x = a
                continue
            if self.name_x is None and b in h.axes.name: 
                self.name_x = b
                continue            
            if self.name_y is None and a in h.axes.name: 
                self.name_y = a
                continue
            if self.name_y is None and b in h.axes.name: 
                self.name_y = b
                continue
        if self.name_x is None or self.name_y is None:
            raise RuntimeError(f"Can not find ABCD axes among histogram axes {h.axes.name}")
        logger.debug(f"Found ABCD axes {self.name_x} and {self.name_y}")

    # set slices object for selection of signal and sideband regions
    def set_selections_x(self):
        if self.name_x.startswith("pass"):
            self.sel_x = 1
            self.sel_dx = 0
        else:
            x0, x1, x2, x3 = get_selection_edges(self.name_x)
            s = hist.tag.Slicer()
            do = hist.sum if self.integrate_x else None
            self.sel_x = s[x0:x1:do] if x0 is not None and x1.imag > x0.imag else s[x1:x0:do]
            self.sel_dx = s[x1:x3:hist.sum] if x3 is None or x3.imag > x1.imag else s[x3:x1:hist.sum]

    def set_selections_y(self):
        if self.name_y.startswith("pass"):
            self.sel_y = 1
            self.sel_dy = 0
        else:
            y0, y1, y2, y3 = get_selection_edges(self.name_y, upper_bound=self.upper_bound_y)
            s = hist.tag.Slicer()
            self.sel_y = s[y0:y1:hist.sum] if y0 is not None and y1.imag > y0.imag else s[y1:y0:hist.sum]
            self.sel_dy = s[y1:y3:hist.sum] if y3 is None or y3.imag > y1.imag else s[y3:y1:hist.sum]


class SignalSelectorABCD(HistselectorABCD):
    # simple ABCD method
    def __init__(self, h, *args, **kwargs):
        super().__init__(h, *args, **kwargs)

    # signal region selection
    def get_hist(self, h, is_nominal=False):
        return self.get_hist_passX_passY(h)

class FakeSelectorSimpleABCD(HistselectorABCD):
    # simple ABCD method
    def __init__(self, h, *args, 
        smoothing_mode="full", # 'binned', 'fakerate', or 'full'
        smoothing_order_fakerate=2,
        throw_toys=None, #"normal", # None, 'normal' or 'poisson'
        global_scalefactor=1, # apply global correction factor on prediction
        **kwargs
    ):
        """
        :polynomial: choices: ["power","bernstein"]
        """
        super().__init__(h, *args, **kwargs)

        # nominal histogram to be used to transfer variances for systematic variations
        self.h_nominal = None
        self.global_scalefactor = global_scalefactor

        # select appropriate polynomial depending on type of smoothing
        if smoothing_mode == "fakerate":
            self.polynomial = "bernstein"
        else:
            self.polynomial = "power"

        # rebinning doesn't make sense for binned estimation
        if smoothing_mode == "binned":
            self.rebin_smoothing_axis = None

        if hasattr(self, "fakerate_integration_axes"):
            if smoothing_mode == "full" and self.fakerate_integration_axes:
                raise NotImplementedError("Smoothing of full fake prediction is not currently supported together with integration axes.")

        self.throw_toys = throw_toys

        # fakerate factor
        self.smoothing_mode = smoothing_mode
        self.smoothing_order_fakerate = smoothing_order_fakerate
        if self.smoothing_mode != "binned":
            logger.info(f"Fakerate smoothing order is {self.smoothing_order_fakerate}")

        # solve with non negative least squares
        if self.polynomial=="bernstein":
            self.solve = solve_nonnegative_leastsquare
        else:
            self.solve = solve_leastsquare

        # set smooth functions
        self.f_smoothing = None
        if self.smoothing_mode != "binned":
            self.f_smoothing = get_regression_function(self.smoothing_order_fakerate, pol=self.polynomial)

        # histogram with nonclosure corrections
        self.hCorr = None

    def set_correction(self, hQCD, axes_names=False, flow=True):
        # hQCD is QCD MC histogram before selection (should contain variances)
        # axes_names: the axes names to bin the correction in. If empty make an inclusive correction (i.e. a single number)
        hQCD_rebin = hh.rebinHist(hQCD, self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h

        # mirror charge and eta axes to cope with limited QCD MC stat
        hQCD_rebin = hh.mirrorAxes(hQCD_rebin, ["charge", "eta"])

        if self.smoothing_mode == "fakerate" and (len(self.fakerate_integration_axes) or "charge" not in self.fakerate_axes):
            # TODO: should this be done in general?
            sel = {n: hist.sum for n in self.fakerate_integration_axes}
            if "charge" not in self.fakerate_axes:
                sel["charge"] = hist.sum
            hQCD_rebin = hQCD_rebin[sel]                

        # prediction without smoothing
        d, dvar = self.calculate_fullABCD(hQCD_rebin, flow=flow)
        hPred = hist.Hist(
            *hQCD_rebin[{self.name_x: self.sel_x if not self.integrate_x else hist.sum, self.name_y: self.sel_y}].axes, 
            storage=hQCD_rebin.storage_type()
            )
        hPred.values(flow=flow)[...] = d
        if hPred.storage_type == hist.storage.Weight:
            hPred.variances(flow=flow)[...] = dvar
        
        # truth histogram
        hTruth = self.get_hist_passX_passY(hQCD_rebin)

        if any(n not in axes_names for n in hTruth.axes.name if n not in ["charge"]): # charge is anyway (hard coded) mirrored
            # rebin instead of integrate to keep axes to allow for easy post processing
            s = hist.tag.Slicer()
            hTruth = hTruth[{a.name: s[::hist.rebin(a.size)] for a in hTruth.axes if a.name not in axes_names}]
            hPred = hPred[{a.name: s[::hist.rebin(a.size)] for a in hPred.axes if a.name not in axes_names}]

        self.hCorr = hh.divideHists(hTruth, hPred)

        logger.debug(f"Got QCD MC corrections of {self.hCorr.values()}")

        if np.any(~np.isfinite(self.hCorr.values(flow=flow))):
            logger.warning(f"{sum(~np.isfinite(self.hCorr.values(flow=flow)))} Inf or NaN values in QCD MC nonclosure correction")
        if np.any(~np.isfinite(self.hCorr.variances(flow=flow))):
            logger.warning(f"{sum(~np.isfinite(self.hCorr.values(flow=flow)))} Inf or NaN variances in QCD MC nonclosure correction")

    def transfer_variances(self, h, set_nominal=False):
        if set_nominal:
            self.h_nominal = h.copy()
        elif self.h_nominal is not None:
            h = hh.transfer_variances(h, self.h_nominal)
        elif h.storage_type == hist.storage.Weight:
            logger.warning("Nominal histogram is not set but current histogram has variances, use those")
        else:
            raise RuntimeError(f"Failed to transfer variances")
        return h

    def get_hist(self, h, is_nominal=False, variations_smoothing=False, flow=True):
        idx_x = h.axes.name.index(self.name_x)
        if self.smoothing_mode == "fakerate":
            h = self.transfer_variances(h, set_nominal=is_nominal)
            y_frf, y_frf_var = self.compute_fakeratefactor(h, smoothing=True, syst_variations=variations_smoothing)
            c, cvar = self.get_yields_applicationregion(h)
            d = c * y_frf

            if variations_smoothing:
                dvar = c[..., np.newaxis,np.newaxis] * y_frf_var[...,:,:]
            else:
                # only take bin by bin uncertainty from c region
                dvar = y_frf**2 * cvar
        elif self.smoothing_mode == "full":
            h = self.transfer_variances(h, set_nominal=is_nominal)
            d, dvar = self.calculate_fullABCD_smoothed(h, flow=flow, syst_variations=variations_smoothing)
        elif self.smoothing_mode == "binned":
            # no smoothing of rates
            d, dvar = self.calculate_fullABCD(h, flow=flow)
        else:
            raise ValueError("invalid choice of smoothing_mode")

        # set histogram in signal region
        hSignal = hist.Hist(
            *h[{self.name_x: self.sel_x if not self.integrate_x else hist.sum, self.name_y: self.sel_y}].axes, 
            storage=hist.storage.Double() if variations_smoothing else h.storage_type())
        hSignal.values(flow=flow)[...] = d
        if variations_smoothing and self.smoothing_mode != "binned":
            hSignal = self.get_syst_hist(hSignal, d, dvar, flow=flow)
        elif hSignal.storage_type == hist.storage.Weight:
            hSignal.variances(flow=flow)[...] = dvar

        if self.global_scalefactor != 1:
            hSignal = hh.scaleHist(hSignal, self.global_scalefactor)

        return hSignal            

    def get_yields_applicationregion(self, h, flow=True):
        hC = self.get_hist_passX_failY(h)
        c = hC.values(flow=flow)
        if h.storage_type == hist.storage.Weight:
            cvar = hC.variances(flow=flow)
            return c, cvar
        return c, None

    def compute_fakeratefactor(self, h, smoothing=False, syst_variations=False, flow=True, auxiliary_info=False):
        # rebin in smoothing axis to have stable ratios
        sel = {n: hist.sum for n in self.fakerate_integration_axes}
        hNew = hh.rebinHist(h[sel], self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h[sel]

        # select sideband regions
        ha = self.get_hist_failX_failY(hNew)
        hb = self.get_hist_failX_passY(hNew)

        a = ha.values(flow=flow)
        b = hb.values(flow=flow)
        # fakerate factor
        y = divide_arrays(b,a)
        if h.storage_type == hist.storage.Weight:
            avar = ha.variances(flow=flow)
            bvar = hb.variances(flow=flow)
            y_var = bvar/a**2 + b**2*avar/a**4
            y_var[abs(a) < 1]=0

        if self.hCorr:
            # apply QCD MC nonclosure correction and account for variance of correction
            cval = self.hCorr.values(flow=flow)
            y *= cval
            if h.storage_type == hist.storage.Weight:
                cvar = self.hCorr.variances(flow=flow)
                y_var = y**2 * cvar + cval**2 * y_var

        if smoothing:
            x = self.get_bin_centers_smoothing(hNew, flow=True) # the bins where the smoothing is performed (can be different to the bins in h)
            y, y_var = self.smoothen(h, x, y, y_var, syst_variations=syst_variations, auxiliary_info=auxiliary_info, flow=flow)

        # broadcast abcd-x axis and application axes
        slices=[slice(None) if n in ha.axes.name else np.newaxis for n in h[{self.name_x: self.sel_x}].axes.name if n != self.name_y]
        y = y[*slices]
        y_var = y_var[*slices] if y_var is not None else None

        return y, y_var

    def smoothen(self, h, x, y, y_var, syst_variations=False, auxiliary_info=False, flow=True):
        if h.storage_type == hist.storage.Weight:
            # transform with weights
            w = 1/np.sqrt(y_var)
        else:
            logger.warning("Smoothing extended ABCD on histogram without uncertainties, make an unweighted linear squared solution.")
            w = np.ones_like(y)

        # move smoothing axis to last
        axes = [n for n in h.axes.name if n not in [self.name_x, self.name_y, *self.fakerate_integration_axes] ]
        idx_ax_smoothing = axes.index(self.smoothing_axis_name)
        if idx_ax_smoothing != len(axes)-1:
            y = np.moveaxis(y, idx_ax_smoothing, -1)
            w = np.moveaxis(w, idx_ax_smoothing, -1)

        # smooth frf (e.g. in pT)
        X, XTY = get_parameter_matrices(x, y, w, self.smoothing_order_fakerate, pol=self.polynomial)

        params, cov = self.solve(X, XTY)

        # evaluate in range of original histogram
        x_smooth_orig = self.get_bin_centers_smoothing(h, flow=True)
        y_smooth_orig = self.f_smoothing(x_smooth_orig, params)

        if syst_variations:
            y_smooth_var_orig = self.make_eigenvector_predictons_frf(params, cov, x_smooth_orig)
        else: 
            y_smooth_var_orig = None

        # move smoothing axis to original positon again
        if idx_ax_smoothing != len(axes)-1:
            y_smooth_orig = np.moveaxis(y_smooth_orig, -1, idx_ax_smoothing)
            y_smooth_var_orig = np.moveaxis(y_smooth_var_orig, -3, idx_ax_smoothing) if syst_variations else None

        # check for negative rates
        if np.sum(y_smooth_orig<0) > 0:
            logger.warning(f"Found {np.sum(y_smooth_orig<0)} bins with negative fake rate factors")
        if y_smooth_var_orig is not None and np.sum(y_smooth_var_orig<0) > 0:
            logger.warning(f"Found {np.sum(y_smooth_var_orig<0)} bins with negative fake rate factor variations")

        if auxiliary_info:
            y_pred = self.f_smoothing(x, params)
            # flatten
            y_pred = y_pred.reshape(y.shape) 
            w = w.reshape(y.shape)
            chi2, ndf = compute_chi2(y, y_pred, w, nparams=params.shape[-1])
            return y_smooth_orig, y_smooth_var_orig, params, cov, chi2, ndf
        else:
            return y_smooth_orig, y_smooth_var_orig

    def get_syst_hist(self, hNominal, values, alternate, flow=True):
        # return systematic histogram with nominal and nominal+variation on diagonal elements for systematic axes of parameters and up/down variations
        hsyst = hist.Hist(*hNominal.axes, hist.axis.Integer(0, alternate.shape[-2], name="_param", overflow=False, underflow=False), common.down_up_axis, storage=hist.storage.Double())

        # check alternate lower than 0
        # vcheck = alternate < (-1*values[...,np.newaxis,np.newaxis])
        # if np.sum(vcheck) > 0:
        #     logger.warning(f"Found {np.sum(vcheck)} bins with alternate giving negative yields, set these to 0")
        #     alternate[vcheck] = 0 

        hsyst.values(flow=flow)[...] = alternate-values[...,np.newaxis,np.newaxis]

        # decorrelate in fakerate axes
        axes_names = [n for n in self.fakerate_axes if self.smoothing_mode=="binned" or n != self.smoothing_axis_name]
        hsyst = hh.expand_hist_by_duplicate_axes(hsyst, axes_names, [f"_{n}" for n in axes_names])    

        # add nominal hist and broadcast
        hNominal = hh.addHists(hNominal, hsyst)
        return hNominal

    def make_eigenvector_predictons_frf(self, params, cov, x1, x2=None):
        return make_eigenvector_predictons(params, cov, func=self.f_smoothing, x1=x1, x2=x2, force_positive=False)#self.polynomial=="bernstein")

    def get_bin_centers_smoothing(self, h, flow=True):
        return self.get_bin_centers(h, self.smoothing_axis_name, self.smoothing_axis_min, self.smoothing_axis_max, flow=flow)

    def get_bin_centers(self, h, axis_name, xmin=None, xmax=None, flow=True, cap=False):
        # get bin centers for interpolation/smoothing
        x = h.axes[axis_name].centers
        if flow:
            x = extend_edges(h.axes[axis_name].traits, x)

        if self.polynomial=="bernstein":
            # transform bernstein polinomials to [0,1]
            x = (x - xmin) / (xmax - xmin)
            if np.sum(x < 0) or np.sum(x > 1):
                if cap:
                    logger.info("Values outside [0,1] found, those will be capped to [0,1]")
                    x[x < 0] = 0
                    x[x > 1] = 1
                else:
                    raise RuntimeError(f"All values need to be within [0,1] but {np.sum(x < 0)} values smaller 0 ({x[x < 0]}) and {np.sum(x > 1)} larger 1 ({x[x > 1]}) found after transformation with xmin={xmin} and xmax={xmax}")
        return x

    def calculate_fullABCD(self, h, flow=True):
        if len(self.fakerate_integration_axes) > 0:
            logger.warning(f"Binned fake estimation is performed but fakerate integration axes {self.fakerate_integration_axes} are set, the bin-by-bin stat uncertainties are not correct along this axis.")
        if not self.integrate_x:
            logger.warning(f"Binned fake estimation is performed but ABCD x-axis is not integrated, the bin-by-bin stat uncertainties are not correct along this axis.")

        c, cvar = self.get_yields_applicationregion(h)
        frf, frf_var = self.compute_fakeratefactor(h)

        d = c * frf
        if h.storage_type == hist.storage.Weight:
            dvar = frf**2 * cvar + c**2 * frf_var
        
        return d, dvar

    def calculate_fullABCD_smoothed(self, h, flow=True, syst_variations=False):
        hNew = hh.rebinHist(h, self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h

        dsel, dvarsel = self.calculate_fullABCD(hNew, flow=flow)

        if self.hCorr:
            # apply QCD MC nonclosure correction and account for variance of correction
            cval = self.hCorr.values(flow=flow)
            cvar = self.hCorr.variances(flow=flow)
            dsel *= cval
            dvarsel = dsel**2 * cvar + cval**2 * dvarsel

        smoothidx = hNew.axes.name.index(self.smoothing_axis_name)
        smoothing_axis = hNew.axes[self.smoothing_axis_name]
        nax = dsel.ndim

        # underflow and overflow are left unchanged along the smoothing axis
        # so we need to exclude them if they have been otherwise included
        if flow:
            smoothstart = 1 if smoothing_axis.traits.underflow else 0
            smoothstop = -1 if smoothing_axis.traits.overflow else None
            smoothslice = slice(smoothstart, smoothstop)
        else:
            smoothslice = slice(None)

        sel = nax*[slice(None)]
        sel[smoothidx] = smoothslice

        dsel =  dsel[*sel]
        dvarsel = dvarsel[*sel]

        xwidth = hNew.axes[self.smoothing_axis_name].widths
        xwidthtgt = h.axes[self.smoothing_axis_name].widths

        xwidth = xwidth[*smoothidx*[None], :, *(nax - smoothidx - 1)*[None]]
        xwidthtgt = xwidthtgt[*smoothidx*[None], :, *(nax - smoothidx - 1)*[None]]#

        dsel *= 1./xwidth
        dvarsel *= 1./xwidth**2

        logd = np.where(dsel>0., np.log(dsel), 0.)
        logdvar = np.where(dsel>0., dvarsel/dsel**2, 1.)
        x = smoothing_axis.centers

        if syst_variations:
            logd, logdvar = self.smoothen(h, x, logd, logdvar, syst_variations=syst_variations)
        else:
            logd, logdvar, params, cov, chi2, ndf = self.smoothen(h, x, logd, logdvar, syst_variations=syst_variations, auxiliary_info=True)
            logger.debug("ndof", ndf)
            logger.debug("chisq stats", np.mean(chi2), np.std(chi2), np.min(chi2), np.max(chi2))
            # print(f"chisq/ndof = {chi2}/{ndf}")
        # get output shape from original hist axes, but as for result histogram
        d = np.zeros([a.extent if flow else a.shape for a in h[{self.name_x:self.sel_x if not self.integrate_x else hist.sum}].axes if a.name != self.name_y])
        # leave the underflow and overflow unchanged if present
        d[*sel] = np.exp(logd)*xwidthtgt
        if syst_variations:
            dvar = d
            dvar = dvar[..., None, None]*np.ones((*dvar.shape, *logdvar.shape[-2:]), dtype=dvar.dtype)
            # leave the underflow and overflow unchanged if present
            dvar[*sel, :, :] = np.exp(logdvar)*xwidthtgt[..., None, None]**2
        else:
            # with full smoothing all of the statistical uncertainty is included in the
            # explicit variations, so the remaining binned uncertainty is zero
            dvar = np.zeros_like(d)

        return d, dvar

class FakeSelectorSimultaneousABCD(FakeSelectorSimpleABCD):
    # simple ABCD method simultaneous fitting all regions
    def __init__(self, h, *args, **kwargs):
        super().__init__(h, *args, **kwargs)

    def get_hist(self, h):
        if h.storage_type == hist.storage.Weight:
            # setting errors to 0
            h.view(flow=True)[...] = np.stack((h.values(flow=True), np.zeros_like(h.values(flow=True))), axis=-1)

        if self.name_x not in h.axes.name or self.name_y not in h.axes.name:
            raise RuntimeError(f'{self.name_x} and {self.name_y} expected to be found in histogram, but only have axes {h.axes.name}')

        # axes in the correct ordering
        axes = [ax for ax in h.axes.name if ax not in [self.name_x, self.name_y]]
        axes += [self.name_x, self.name_y]

        if set(h.axes.name) != set(axes):
            logger.warning(f"Axes in histogram '{h.axes.name}' are not the same as required '{axes}' or in a different order than expected, try to project")
            h = h.project(*axes)

        # set the expected values in the signal region
        slices = [self.sel_x if n==self.name_x else self.sel_y if n==self.name_y else slice(None) for n in h.axes.name]
        h.values(flow=True)[*slices] = super().get_hist(h).values(flow=True)

        if self.global_scalefactor != 1:
            h = hh.scaleHist(h, self.global_scalefactor)

        return h

class FakeSelectorExtrapolateABCD(FakeSelectorSimpleABCD):
    # extrapolate the fakerate in the abcd x axis by finding an analytic description in the dx region
    def __init__(self, h, *args, 
        extrapolation_order=1,
        rebin_x="automatic", # can be a list of bin edges, "automatic", or None
        **kwargs
    ):
        super().__init__(h, *args, **kwargs)
        self.set_selections_x()

        self.extrapolation_order = extrapolation_order

        if rebin_x == "automatic":
            edges = h.axes[self.name_x].edges
            self.rebin_x = get_rebinning(edges, self.name_x)
        else:
            self.rebin_x = rebin_x

        if self.polynomial == "bernstein":
            logger.warning(f"It is not recommended to use {self.polynomial} polynomials for extrapolation.")

        if self.smoothing_mode != "binned":
            raise NotImplementedError("Smooting fakerate is not implemented")

        self.f_smoothing = get_regression_function(self.extrapolation_order, pol=self.polynomial)

    # set slices object for selection of sideband regions
    def set_selections_x(self):
        x0, x1, x2, x3 = get_selection_edges(self.name_x)
        s = hist.tag.Slicer()
        self.sel_x = s[x0:x1:] if x0 is not None and x1.imag > x0.imag else s[x1:x0:]
        self.sel_dx = s[x1:x3:] if x3 is None or x3.imag > x1.imag else s[x3:x1:]

    def get_hist(self, h, is_nominal=False, variations_smoothing=False, flow=True):
        h = self.transfer_variances(h, set_nominal=is_nominal)
        c, cvar = self.get_yields_applicationregion(h)
        if self.integrate_x:
            # can convert syst variations into bin by bin stat; do for the fist call (nominal histogram) 
            logger.debug("Integrate abcd x-axis, treat uncertainty as bin by bin stat")
            y_frf, y_frf_variation = self.compute_fakeratefactor(h, syst_variations=is_nominal)
            d = c * y_frf
            # sum up abcd-x axis (sum up variations to treat uncertainty fully correlated across abcd-x axis bins)
            idx_x = h[{self.name_y: self.sel_dy}].axes.name.index(self.name_x)
            d = d.sum(axis=idx_x)
            if is_nominal:
                d_variation = c[..., np.newaxis,np.newaxis] * y_frf_variation
                d_variation = d_variation.sum(axis=idx_x)
                dvar = np.sum([(d_variation[...,iparam,0] - d)**2 for iparam in range(d_variation.shape[-2])], axis=0) # sum of squares of up variations
                dvar += (y_frf**2 * cvar).sum(axis=idx_x) # add uncorrelated bin by bin uncertainty from application region
            else:
                dvar = None
        else:
            y_frf, y_frf_var = self.compute_fakeratefactor(h, syst_variations=variations_smoothing)
            if variations_smoothing:
                dvar = c[..., np.newaxis,np.newaxis] * y_frf_var
            else:
                # only take bin by bin uncertainty from c region
                dvar = y_frf**2 * cvar
            d = c * y_frf

        # set histogram in signal region
        axes = h[{self.name_x: self.sel_x if not self.integrate_x else hist.sum, self.name_y: self.sel_y}].axes
        hSignal = hist.Hist(*axes, storage=hist.storage.Double() if variations_smoothing else h.storage_type())
        hSignal.values(flow=flow)[...] = d
        if variations_smoothing:
            hSignal = self.get_syst_hist(hSignal, d, dvar, flow=flow)
        elif dvar is not None and hSignal.storage_type == hist.storage.Weight:
            hSignal.variances(flow=flow)[...] = dvar

        if self.global_scalefactor != 1:
            hSignal = hh.scaleHist(hSignal, self.global_scalefactor)

        return hSignal

    def compute_fakeratefactor(self, h, syst_variations=False, flow=True, auxiliary_info=False):
        sel = {n: hist.sum for n in self.fakerate_integration_axes}
        # rebin in regression axis to have stable ratios
        hNew = hh.rebinHist(h[sel], self.name_x, self.rebin_x) if self.rebin_x is not None else h[sel]

        # the underflow and overflow bins are kept by a simple selection, remove them
        hNew = hh.disableFlow(hNew, self.name_x)

        # select sideband regions
        ha = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_dy}]
        hb = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_y}]

        a = ha.values(flow=flow)
        b = hb.values(flow=flow)

        # fakerate factor
        y = divide_arrays(b,a)
        if h.storage_type == hist.storage.Weight:
            # full variances
            avar = ha.variances(flow=flow)
            bvar = hb.variances(flow=flow)
            y_var = bvar/a**2 + b**2*avar/a**4
            y_var[abs(a) <= 0]=0

        # the bins where the regression is performed (can be different to the bin in h)
        x = self.get_bin_centers(ha, self.name_x, flow=False) 
        if len(x)-1 < self.extrapolation_order:
            raise RuntimeError(f"The extrapolation order is {self.extrapolation_order} but it can not be higher than the number of bins in {self.name_x} of {len(x)} -1")
        # if len(x)-1 == self.extrapolation_order:
        #     logger.info("Extrapolation can be done fully analytical")
        #     slope = (y[...,1] - y[...,0]) / (x[1] - x[0])
        #     offset = y[...,0] - slope * x[0]

        #     x_extrapolation = self.get_bin_centers(h[{self.name_x: self.sel_x}], self.name_x, flow=flow)
        #     y = slope[...,np.newaxis] * x_extrapolation + offset[...,np.newaxis]
        
        y, y_var = self.extrapolate_fakerate(h[{**sel, self.name_x: self.sel_x}], x, y, y_var, syst_variations=syst_variations, auxiliary_info=auxiliary_info, flow=flow)

        # TODO: smoothing
        # if self.smooth_fakerate:

        # broadcast abcd-x axis and application axes
        slices=[slice(None) if n in ha.axes.name else np.newaxis for n in h[{self.name_x: self.sel_x}].axes.name if n != self.name_y and (not self.integrate_x or n != self.name_x)]
        y = y[*slices]
        y_var = y_var[*slices] if y_var is not None else None

        return y, y_var

    def extrapolate_fakerate(self, h, x, y, y_var, syst_variations=False, auxiliary_info=False, flow=True):
        if h.storage_type == hist.storage.Weight:
            # transform with weights
            w = 1/np.sqrt(y_var)
        else:
            logger.warning("Extrapolate ABCD on histogram without uncertainties, make an unweighted linear squared solution.")
            w = np.ones_like(y)

        # move regression axis to last
        axes = [n for n in h.axes.name if n not in [self.name_y,] ]
        idx_ax_regress = axes.index(self.name_x)
        if idx_ax_regress != len(axes)-1:
            y = np.moveaxis(y, idx_ax_regress, -1)
            w = np.moveaxis(w, idx_ax_regress, -1)

        # regression frf (e.g. in mT)
        X, XTY = get_parameter_matrices(x, y, w, self.extrapolation_order, pol=self.polynomial)

        params, cov = self.solve(X, XTY)

        # evaluate in range of application region of original histogram
        x_extrapolation = self.get_bin_centers(h, self.name_x, flow=flow)
        y_extrapolation = self.f_smoothing(x_extrapolation, params)

        if syst_variations:
            y_extrapolation_var = self.make_eigenvector_predictons_frf(params, cov, x_extrapolation)
        else: 
            y_extrapolation_var = None

        if idx_ax_regress != len(axes)-1:
            # move extrapolation axis to original positon again
            y_extrapolation = np.moveaxis(y_extrapolation, -1, idx_ax_regress)
            y_extrapolation_var = np.moveaxis(y_extrapolation_var, -3, idx_ax_regress) if syst_variations else None

        # check for negative rates
        if np.sum(y_extrapolation<0) > 0:
            logger.warning(f"Found {np.sum(y_extrapolation<0)} bins with negative fake rate factors")
        if y_extrapolation_var is not None and np.sum(y_extrapolation_var<0) > 0:
            logger.warning(f"Found {np.sum(y_extrapolation_var<0)} bins with negative fake rate factor variations")

        if auxiliary_info:
            y_pred = self.f_smoothing(x, params)
            # flatten
            y_pred = y_pred.reshape(y.shape) 
            w = w.reshape(y.shape)
            chi2, ndf = compute_chi2(y, y_pred, w, nparams=params.shape[-1])
            return y_extrapolation, y_extrapolation_var, params, cov, chi2, ndf
        else:
            return y_extrapolation, y_extrapolation_var

class FakeSelector1DExtendedABCD(FakeSelectorSimpleABCD):
    # extended ABCD method with 5 control regions as desribed in https://arxiv.org/abs/1906.10831 equation 16
    def __init__(self, h, *args, **kwargs):
        super().__init__(h, *args, **kwargs)
        self.sel_d2x = None    
        self.set_selections_x(integrate_x=self.integrate_x)

    # set slices object for selection of sideband regions
    def set_selections_x(self, integrate_x=True):
        x0, x1, x2, x3 = get_selection_edges(self.name_x)
        s = hist.tag.Slicer()
        do = hist.sum if integrate_x else None
        self.sel_x = s[x0:x1:do] if x0 is not None and x1.imag > x0.imag else s[x1:x0:do]
        self.sel_dx = s[x1:x2:hist.sum] if x2.imag > x1.imag else s[x2:x1:hist.sum]
        self.sel_d2x = s[x2:x3:hist.sum] if x3.imag > x2.imag else s[x3:x2:hist.sum]

    def calculate_fullABCD(self, h, flow=True, syst_variations=False):
        if len(self.fakerate_integration_axes) > 0:
            logger.warning(f"Binned fake estimation is performed but fakerate integration axes {self.fakerate_integration_axes} are set, the bin-by-bin stat uncertainties are not correct along this axis.")
        if not self.integrate_x:
            logger.warning(f"Binned fake estimation is performed but ABCD x-axis is not integrated, the bin-by-bin stat uncertainties are not correct along this axis.")

        sel = {n: hist.sum for n in self.fakerate_integration_axes}

        ha = h[{**sel, self.name_x: self.sel_dx, self.name_y: self.sel_dy}]
        hax = h[{**sel, self.name_x: self.sel_d2x, self.name_y: self.sel_dy}]
        hb = h[{**sel, self.name_x: self.sel_dx, self.name_y: self.sel_y}]
        hbx = h[{**sel, self.name_x: self.sel_d2x, self.name_y: self.sel_y}]

        hc = h[{self.name_x: self.sel_x, self.name_y: self.sel_dy}]

        # calculate extended ABCD method without smoothing and in a single bin in x
        a = ha.values(flow=flow)
        ax = hax.values(flow=flow)
        b = hb.values(flow=flow)
        bx = hbx.values(flow=flow)
        c = hc.values(flow=flow)

        frf = (b/a)**2 * (ax/bx)
        slices=[slice(None) if a in ha.axes else np.newaxis for a in hc.axes]

        frf = frf[*slices]
        d = c * frf

        dvar=None
        if h.storage_type == hist.storage.Weight:
            avar = ha.variances(flow=flow)
            axvar = hax.variances(flow=flow)
            bvar = hb.variances(flow=flow)
            bxvar = hbx.variances(flow=flow)
            cvar = hc.variances(flow=flow)
            dvar = frf**2 * cvar + d**2 * (4 * bvar/b**2 + 4 * avar/a**2 + axvar/ax**2 + bxvar/bx**2)[*slices]
        
        return d, dvar

    def compute_fakeratefactor(self, h, smoothing=False, syst_variations=False, flow=True, auxiliary_info=False):
        # rebin in smoothing axis to have stable ratios
        sel = {n: hist.sum for n in self.fakerate_integration_axes}
        hNew = hh.rebinHist(h[sel], self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h[sel]

        # select sideband regions
        ha = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_dy}]
        hax = hNew[{self.name_x: self.sel_d2x, self.name_y: self.sel_dy}]
        hb = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_y}]
        hbx = hNew[{self.name_x: self.sel_d2x, self.name_y: self.sel_y}]

        a = ha.values(flow=flow)
        ax = hax.values(flow=flow)
        b = hb.values(flow=flow)
        bx = hbx.values(flow=flow)

        # fakerate factor
        y_num = ax*b**2
        y_den = bx*a**2
        # fakerate factor
        y = divide_arrays(y_num,y_den)

        if h.storage_type == hist.storage.Weight:
            # full variances
            avar = ha.variances(flow=flow)
            axvar = hax.variances(flow=flow)
            bvar = hb.variances(flow=flow)
            bxvar = hbx.variances(flow=flow)
            y_var = b**4/(bx**2*a**4)*axvar + ax**2*b**2/(bx**2*a**4)*4*bvar + 4*avar/a**2 + ax**2*b**4/(bx**4*a**4)*bxvar
            y_var[abs(y_den) <= 0]=0

        if self.hCorr:
            # apply QCD MC nonclosure correction and account for variance of correction
            cval = self.hCorr.values(flow=flow)
            y *= cval
            if h.storage_type == hist.storage.Weight:
                cvar = self.hCorr.variances(flow=flow)
                y_var = y**2 * cvar + cval**2 * y_var

        if self.throw_toys:
            logger.info("Throw toys")
            # throw toys for each parameter separately
            values = np.stack([a, ax, b, bx], axis=-1)
            variances = np.stack([avar, axvar, bvar, bxvar], axis=-1)

            nsamples=10000
            seed=42

            if self.throw_toys == "normal":
                # throw gaussian toys
                toy_shape = [*values.shape, nsamples]
                toy_size=np.product(toy_shape)
                np.random.seed(seed)  # For reproducibility
                toys = np.random.normal(0, 1, size=toy_size)
                toys = np.reshape(toys, toy_shape)
                toys = toys*np.sqrt(variances)[...,np.newaxis] + values[...,np.newaxis]
                # toy = toys[...,0,:] / toys[...,1,:]
                toy = (toys[...,1,:]*toys[...,2,:]**2) / (toys[...,0,:]**2 * toys[...,3,:])
                toy_mean = np.mean(toy, axis=-1)
                toy_var = np.var(toy, ddof=1, axis=-1) 
            elif self.throw_toys == "poisson":
                # throw posson toys
                toy_shape = [nsamples, *values.shape]
                rng = np.random.default_rng(seed)
                toys = rng.poisson(values, size=toy_shape)
                toy = (toys[...,1]*toys[...,2]**2) / (toys[...,0]**2 * toys[...,3])
                toy_mean = np.mean(toy, axis=0)
                toy_var = np.var(toy, ddof=1, axis=0) 

            y = toy_mean
            y_var = toy_var

            logger.info("Done with toys")


        if smoothing:
            x = self.get_bin_centers_smoothing(hNew, flow=True) # the bins where the smoothing is performed (can be different to the bin in h)
            y, y_var = self.smoothen(h, x, y, y_var, syst_variations=syst_variations, auxiliary_info=auxiliary_info, flow=flow)

        # broadcast abcd-x axis and application axes
        slices=[slice(None) if n in ha.axes.name else np.newaxis for n in h[{self.name_x: self.sel_x}].axes.name if n != self.name_y]
        y = y[*slices]
        y_var = y_var[*slices] if y_var is not None else None

        return y, y_var

class FakeSelector2DExtendedABCD(FakeSelector1DExtendedABCD):
    # extended ABCD method with 8 control regions as desribed in https://arxiv.org/abs/1906.10831 equation 15
    def __init__(self, h, *args, 
        interpolate_x=True, 
        interpolation_order=2,
        rebin_x="automatic", # can be a list of bin edges, "automatic", or None
        integrate_shapecorrection_x=False, smooth_shapecorrection=True, 
        smoothing_order_shapecorrection=[2,2,2],
        **kwargs
    ):
        super().__init__(h, *args, **kwargs)
        self.sel_d2y = None
        self.set_selections_y()
        self.set_selections_x(integrate_x=False)

        self.interpolate_x = interpolate_x
        self.interpolation_order = interpolation_order

        if rebin_x == "automatic":
            edges = h.axes[self.name_x].edges
            self.rebin_x = get_rebinning(edges, self.name_x)
        else:
            self.rebin_x = rebin_x

        # min and max (in application region) for transformation of bernstain polynomials into interval [0,1]
        self.axis_x_min = h[{self.name_x: self.sel_x}].axes[self.name_x].edges[0]
        if self.name_x == "mt":
            # mt does not have an upper bound, cap at 100
            edges = self.rebin_x if self.rebin_x is not None else h.axes[self.name_x].edges
            self.axis_x_max = extend_edges(h.axes[self.name_x].traits, edges)[-1]
        elif self.name_x in ["iso", "relIso", "relJetLeptonDiff", "dxy"]:
            # iso and dxy have a finite lower and upper bound in the application region
            self.axis_x_max = abcd_thresholds[self.name_x][1]
        else:
            self.axis_x_max = self.axis_x.edges[-1]

        # shape correction, can be interpolated in the abcd x-axis in 1D, in the x-axis and smoothing axis in 2D, or in the smoothing axis integrating out the abcd x-axis in 1D
        self.integrate_shapecorrection_x = integrate_shapecorrection_x # if the shape correction factor for the abcd x-axis should be inclusive or differential
        if self.integrate_shapecorrection_x and self.interpolate_x:
            raise RuntimeError("Can not integrate and interpolate x at the same time")
        self.smooth_shapecorrection = smooth_shapecorrection
        self.smoothing_order_shapecorrection = smoothing_order_shapecorrection

        # set smooth functions
        self.f_scf = None
        if self.interpolate_x and self.smooth_shapecorrection: 
            self.f_scf = get_regression_function(self.interpolation_order, self.smoothing_order_shapecorrection, pol=self.polynomial)
        elif self.interpolate_x:
            self.f_scf = get_regression_function(self.interpolation_order, pol=self.polynomial)
        elif self.smooth_shapecorrection:
            self.f_scf = get_regression_function(self.smoothing_order_shapecorrection, pol=self.polynomial)            

    # set slices object for selection of sideband regions
    def set_selections_y(self):
        y0, y1, y2, y3 = get_selection_edges(self.name_y, upper_bound=self.upper_bound_y)
        s = hist.tag.Slicer()
        self.sel_y = s[y0:y1:hist.sum] if y0 is not None and y1.imag > y0.imag else s[y1:y0:hist.sum]
        self.sel_dy = s[y1:y2:hist.sum] if y2.imag > y1.imag else s[y2:y1:hist.sum]
        self.sel_d2y = s[y2:y3:hist.sum] if y3 is None or y3.imag > y2.imag else s[y3:y2:hist.sum]

    def get_hist(self, h, is_nominal=False, variations_scf=False, variations_smoothing=False, variations_full=False, flow=True):
        if variations_scf and variations_smoothing:
            raise RuntimeError(f"Can only calculate vairances for fakerate factor or shape correction factor but not both")

        if self.smoothing_mode=="fakerate" or self.interpolate_x or self.smooth_shapecorrection:
            h = self.transfer_variances(h, set_nominal=is_nominal)

            y_frf, y_frf_var = self.compute_fakeratefactor(h, smoothing=True, syst_variations=variations_smoothing)
            y_scf, y_scf_var = self.compute_shapecorrection(h, syst_variations=variations_scf)
            c, cvar = self.get_yields_applicationregion(h)

            d = c * y_scf * y_frf

            if variations_scf and (self.interpolate_x or self.smooth_shapecorrection):
                dvar = c[..., np.newaxis,np.newaxis] * y_frf[...,np.newaxis,np.newaxis] * y_scf_var[...,:,:]
            elif variations_smoothing and self.smoothing_mode=="fakerate":
                dvar = c[..., np.newaxis,np.newaxis] * y_scf[..., np.newaxis,np.newaxis] * y_frf_var[...,:,:]
            elif self.smooth_shapecorrection or self.interpolate_x:
                # only take bin by bin uncertainty from c region
                dvar = (y_scf * y_frf)**2 * cvar
            else:
                # # take bin by bin uncertainty from c * c/cy 
                dvar = y_frf**2 *( 4*(c/cy)**2 * cvar + (c/cy**2)**2 * cyvar )

            if self.integrate_x:
                idx_x = [n for n in h.axes.name if n != self.name_y].index(self.name_x)
                d = d.sum(axis=idx_x)
                dvar = dvar.sum(axis=idx_x)
        elif self.smoothing_mode == "full":
            h = self.transfer_variances(h, set_nominal=is_nominal)
            d, dvar = self.calculate_fullABCD_smoothed(h, flow=flow, syst_variations=variations_smoothing)
        elif self.smoothing_mode == "binned":
            # no smoothing of rates
            d, dvar = self.calculate_fullABCD(h)
        else:
            raise ValueError("invalid choice of smoothing mode")

        # set histogram in signal region
        axes = [a for a in h[{self.name_x:self.sel_x if not self.integrate_x else hist.sum}].axes if a.name != self.name_y]
        hSignal = hist.Hist(*axes, storage=hist.storage.Double() if variations_smoothing else h.storage_type())
        hSignal.values(flow=flow)[...] = d
        if variations_scf or variations_smoothing or variations_full:
            hSignal = self.get_syst_hist(hSignal, d, dvar, flow=flow)
        elif hSignal.storage_type == hist.storage.Weight:
            hSignal.variances(flow=flow)[...] = dvar

        if self.global_scalefactor != 1:
            hSignal = hh.scaleHist(hSignal, self.global_scalefactor)

        return hSignal

    def compute_shapecorrection(self, h, syst_variations=False, apply=False, flow=True, auxiliary_info=False):
        # if apply=True, shape correction is multiplied to application region for correct statistical uncertainty, only allowed if not smoothing
        if apply and (self.interpolate_x or self.smooth_shapecorrection):
            raise NotImplementedError(f"Direct application of shapecorrection only supported when no smoothing is performed")
        # rebin in smoothing axis to have stable ratios
        sel = {n: hist.sum for n in self.fakerate_integration_axes}
        hNew = hh.rebinHist(h[sel], self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h[sel]
        if self.rebin_x is not None and (self.interpolate_x or self.smooth_shapecorrection):
            hNew = hh.rebinHist(hNew, self.name_x, self.rebin_x) # only rebin for regression
        hNew = hNew[{self.name_x: self.sel_x}]

        hc = hNew[{self.name_y: self.sel_dy}]
        hcy = hNew[{self.name_y: self.sel_d2y}]

        c = hc.values(flow=flow)
        cy = hcy.values(flow=flow)
        if h.storage_type == hist.storage.Weight:
            cvar = hc.variances(flow=flow)
            cyvar = hcy.variances(flow=flow)

        if self.integrate_shapecorrection_x:
            idx_x = hNew.axes.name.index(self.name_x)
            c = c.sum(axis=idx_x)
            cy = cy.sum(axis=idx_x)
            cvar = cvar.sum(axis=idx_x)
            cyvar = cyvar.sum(axis=idx_x)

        # shape correction factor
        y_num = c**2 if apply else c
        y_den = cy
        y = divide_arrays(y_num,y_den)

        if self.throw_toys:
            logger.info("Throw toys")
            # throw toys for nominator and denominator
            y_num_var = cvar
            y_den_var = cyvar
            values = np.stack([y_num, y_den], axis=-1)
            variances = np.stack([y_num_var, y_den_var], axis=-1)

            nsamples=10000
            seed=42

            if self.throw_toys == "normal":
                # throw gaussian toys
                toy_shape = [*values.shape, nsamples]
                toy_size=np.product(toy_shape)
                np.random.seed(seed)  # For reproducibility
                toys = np.random.normal(0, 1, size=toy_size)
                toys = np.reshape(toys, toy_shape)
                toys = toys*np.sqrt(variances)[...,np.newaxis] + values[...,np.newaxis]
                toy = toys[...,0,:] / toys[...,1,:]
                toy_mean = np.mean(toy, axis=-1)
                toy_var = np.var(toy, ddof=1, axis=-1) 
            elif self.throw_toys == "poisson":
                # throw posson toys
                toy_shape = [nsamples, *values.shape]
                rng = np.random.default_rng(seed)
                toys = rng.poisson(values, size=toy_shape)
                toy = toys[...,0] / toys[...,1]
                toy_mean = np.mean(toy, axis=0)
                toy_var = np.var(toy, ddof=1, axis=0) 

            y = toy_mean
            y_var = toy_var

            logger.info("Done with toys")

        if self.interpolate_x or self.smooth_shapecorrection:
            if h.storage_type == hist.storage.Weight:
                y_var = cvar/cy**2 + (c**2 * cyvar)/cy**4 # multiply out for better numerical stability (avoid NaNs from divisions)
                w = 1/np.sqrt(y_var)
            else:
                logger.warning("Smoothing extended ABCD on histogram without uncertainties, make an unweighted linear squared solution.")
                w = np.ones_like(y)
        else:
            if h.storage_type == hist.storage.Weight:
                if apply:
                    y_var = 4*c**2/cy**2*cvar + c**4/cy**4*cyvar # multiply out for better numerical stability (avoid NaNs from divisions)
                else:
                    y_var = 1/cy**2*cvar + c**2/cy**4*cyvar # multiply out for better numerical stability (avoid NaNs from divisions)

        if self.interpolate_x:
            axes = [n for n in h.axes.name if n not in [*self.fakerate_integration_axes, self.name_y] ]

            idx_ax_interpol = axes.index(self.name_x)

            # constract x points in application region
            x_interpol = self.get_bin_centers_interpolation(hNew, flow=flow)
            # x axis in original binning in passing region
            x_interpol_orig = self.get_bin_centers_interpolation(h[{self.name_x: self.sel_x}], flow=flow, cap=True)

            if self.smooth_shapecorrection:
                # interpolate scf in mT and smooth in pT (i.e. 2D)
                idx_ax_smoothing = hNew.axes.name.index(self.smoothing_axis_name)
                if idx_ax_smoothing != len(axes)-2 or idx_ax_interpol != len(axes)-1:
                    y = np.moveaxis(y, (idx_ax_smoothing, idx_ax_interpol), (-2, -1))
                    w = np.moveaxis(w, (idx_ax_smoothing, idx_ax_interpol), (-2, -1))

                x_smoothing = self.get_bin_centers_smoothing(hNew, flow=True)

                X, y, XTY = get_parameter_matrices_from2D(x_interpol, x_smoothing, y, w, 
                    self.interpolation_order, self.smoothing_order_shapecorrection, pol=self.polynomial, flatten=True)

                params, cov = self.solve(X, XTY)

                x_smooth_orig = self.get_bin_centers_smoothing(h, flow=True)
                y_smooth_orig = self.f_scf(x_interpol_orig, x_smooth_orig, params)

                if syst_variations:
                    y_smooth_var_orig = self.make_eigenvector_predictons_scf(params, cov, x_interpol_orig, x_smooth_orig)
                else: 
                    y_smooth_var_orig = None

                # move interpolation axis to original positon again
                if idx_ax_smoothing != len(axes)-2 or idx_ax_interpol != len(axes)-1:
                    y_smooth_orig = np.moveaxis(y_smooth_orig, (-2, -1), (idx_ax_smoothing, idx_ax_interpol))
                    y_smooth_var_orig = np.moveaxis(y_smooth_var_orig, (-4, -3), (idx_ax_smoothing, idx_ax_interpol)) if syst_variations else None

                if auxiliary_info:
                    y_pred = self.f_scf(x_interpol, x_smoothing, params)
            else:
                # interpolate scf in mT in 1D
                if idx_ax_interpol != len(axes)-1:
                    y = np.moveaxis(y, idx_ax_interpol, -1)
                    w = np.moveaxis(w, idx_ax_interpol, -1)

                X, XTY = get_parameter_matrices(x_interpol, y, w, self.smoothing_order_fakerate, pol=self.polynomial)
                params, cov = self.solve(X, XTY)

                y_smooth_orig = self.f_scf(x_interpol_orig, params)

                if syst_variations:
                    y_smooth_var_orig = self.make_eigenvector_predictons_scf(params, cov, x_interpol_orig)
                else: 
                    y_smooth_var_orig = None

                if auxiliary_info:
                    y_pred = self.f_scf(x_interpol, params)

                # move interpolation axis to original positon again
                if idx_ax_interpol != len(axes)-1:
                    y_smooth_orig = np.moveaxis(y_smooth_orig, -1, idx_ax_interpol)
                    y_smooth_var_orig = np.moveaxis(y_smooth_var_orig, -3, idx_ax_smoothing) if syst_variations else None
        
            # check for negative rates
            if np.sum(y_smooth_orig<0) > 0:
                logger.warning(f"Found {np.sum(y_smooth_orig<0)} bins with negative shape correction factors")
            if y_smooth_var_orig is not None and np.sum(y_smooth_var_orig<0) > 0:
                logger.warning(f"Found {np.sum(y_smooth_var_orig<0)} bins with negative shape correction factor variations")
                y_smooth_var_orig[y_smooth_var_orig<0] = 0

            y, y_var = y_smooth_orig, y_smooth_var_orig

            if auxiliary_info:
                # flatten
                y_pred = y_pred.reshape(y.shape) 
                w = w.reshape(y.shape)
                chi2, ndf = compute_chi2(y, y_pred, w, nparams=params.shape[-1])

        elif self.smooth_shapecorrection:
            # don't interpolate in mT, but smooth in pT in 1D

            # move smoothing axis to last
            axes = [n for n in h.axes.name if n not in [*self.fakerate_integration_axes, self.name_x, self.name_y] ]
            idx_ax_smoothing = axes.index(self.smoothing_axis_name)
            if idx_ax_smoothing != len(axes)-1:
                y = np.moveaxis(y, idx_ax_smoothing, -1)
                w = np.moveaxis(w, idx_ax_smoothing, -1)

            # smooth scf (e.g. in pT)
            x_smoothing = self.get_bin_centers_smoothing(hNew, flow=True)
            X, XTY = get_parameter_matrices(x_smoothing, y, w, self.smoothing_order_shapecorrection, pol=self.polynomial)
            params, cov = self.solve(X, XTY)

            # evaluate in range of original histogram
            x_smooth_orig = self.get_bin_centers_smoothing(h, flow=True)
            y_smooth_orig = self.f_scf(x_smooth_orig, params)

            if syst_variations:
                y_smooth_var_orig = self.make_eigenvector_predictons_scf(params, cov, x_smooth_orig)
            else: 
                y_smooth_var_orig = None

            # move smoothing axis to original positon again
            if idx_ax_smoothing != len(axes)-1:
                y_smooth_orig = np.moveaxis(y_smooth_orig, -1, idx_ax_smoothing)
                y_smooth_var_orig = np.moveaxis(y_smooth_var_orig, -3, idx_ax_smoothing) if syst_variations else None

            # check for negative rates
            if np.sum(y_smooth_orig<0) > 0:
                logger.warning(f"Found {np.sum(y_smooth_orig<0)} bins with negative shape correction factors")
            if y_smooth_var_orig is not None and np.sum(y_smooth_var_orig<0) > 0:
                logger.warning(f"Found {np.sum(y_smooth_var_orig<0)} bins with negative shape correction factor variations")

            y, y_var = y_smooth_orig, y_smooth_var_orig

            if auxiliary_info:
                y_pred = self.f_scf(x_smoothing, params)
                # flatten
                y_pred = y_pred.reshape(y.shape) 
                w = w.reshape(y.shape)
                chi2, ndf = compute_chi2(y, y_pred, w, nparams=params.shape[-1])

        # broadcast abcd-x axis and application axes
        if self.integrate_shapecorrection_x:
            slices=[slice(None) if n in hc.axes.name else np.newaxis for n in h.axes.name if n not in [self.name_y]]
            y = y[*slices]
            y_var = y_var[*new_axis] if y_var is not None else None
        else:
            slices=[slice(None) if n in hc.axes.name else np.newaxis for n in h.axes.name if n not in [self.name_x, self.name_y]]
            y = y[*slices]
            y_var = y_var[*slices] if y_var is not None else None

        if auxiliary_info:      
            return y, y_var, params, cov, chi2, ndf
        else:
            return y, y_var


    def compute_fakeratefactor(self, h, smoothing=False, syst_variations=False, flow=True, auxiliary_info=False):
        # rebin in smoothing axis to have stable ratios
        sel = {n: hist.sum for n in self.fakerate_integration_axes}
        hNew = hh.rebinHist(h[sel], self.smoothing_axis_name, self.rebin_smoothing_axis) if self.rebin_smoothing_axis is not None else h[sel]

        # select sideband regions
        ha = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_dy}]
        hax = hNew[{self.name_x: self.sel_d2x, self.name_y: self.sel_dy}]
        hay = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_d2y}]
        haxy = hNew[{self.name_x: self.sel_d2x, self.name_y: self.sel_d2y}]
        hb = hNew[{self.name_x: self.sel_dx, self.name_y: self.sel_y}]
        hbx = hNew[{self.name_x: self.sel_d2x, self.name_y: self.sel_y}]

        a = ha.values(flow=flow)
        ax = hax.values(flow=flow)
        ay = hay.values(flow=flow)
        axy = haxy.values(flow=flow)
        b = hb.values(flow=flow)
        bx = hbx.values(flow=flow)

        # fakerate factor
        y_num = (ax*ay*b)**2
        y_den = a**4 * axy * bx
        y = divide_arrays(y_num,y_den)

        if h.storage_type == hist.storage.Weight:
            # full variances
            avar = ha.variances(flow=flow)
            axvar = hax.variances(flow=flow)
            ayvar = hay.variances(flow=flow)
            axyvar = haxy.variances(flow=flow)
            bvar = hb.variances(flow=flow)
            bxvar = hbx.variances(flow=flow)
            y_var = y**2 * (4*axvar/ax**2 + 4*ayvar/ay**2 + axyvar/axy**2 + 4*bvar/b**2 + 16*avar/a**2 + bxvar/bx**2)
            y_var[abs(y_den) <= 0]=0

        if self.throw_toys:
            logger.info("Throw toys")
            # throw toys for each parameter separately
            values = np.stack([a, ax, ay, axy, b, bx], axis=-1)
            variances = np.stack([avar, axvar, ayvar, axyvar, bvar, bxvar], axis=-1)

            nsamples=10000
            seed=42 # For reproducibility

            if self.throw_toys == "normal":
                # throw gaussian toys
                toy_shape = [*values.shape, nsamples]
                toy_size=np.product(toy_shape)
                np.random.seed(seed)
                toys = np.random.normal(0, 1, size=toy_size)
                toys = np.reshape(toys, toy_shape)
                toys = toys*np.sqrt(variances)[...,np.newaxis] + values[...,np.newaxis]
                toy = (toys[...,1,:]*toys[...,2,:]*toys[...,4,:])**2 / (toys[...,0,:]**4 * toys[...,3,:] * toys[...,5,:])
                toy_mean = np.mean(toy, axis=-1)
                toy_var = np.var(toy, ddof=1, axis=-1) 
            elif self.throw_toys == "poisson":
                # throw posson toys
                toy_shape = [nsamples, *values.shape]
                rng = np.random.default_rng(seed)
                toys = rng.poisson(values, size=toy_shape)
                toys = toys.astype(np.double)
                toy = (toys[...,1]*toys[...,2]*toys[...,4])**2 / (toys[...,0]**4 * toys[...,3] * toys[...,5])
                toy_mean = np.mean(toy, axis=0)
                toy_var = np.var(toy, ddof=1, axis=0) 

            y = toy_mean
            y_var = toy_var

            logger.info("Done with toys")

        if smoothing:
            x = self.get_bin_centers_smoothing(hNew, flow=True) # the bins where the smoothing is performed (can be different to the bin in h)
            y, y_var = self.smoothen(h, x, y, y_var, syst_variations=syst_variations, auxiliary_info=auxiliary_info)

        # broadcast abcd-x axis and application axes
        slices=[slice(None) if n in ha.axes.name else np.newaxis for n in h[{self.name_x: self.sel_x}].axes.name if n != self.name_y]
        y = y[*slices]
        y_var = y_var[*slices] if y_var is not None else None

        return y, y_var

    def make_eigenvector_predictons_scf(self, params, cov, x1, x2=None):
        return make_eigenvector_predictons(params, cov, func=self.f_scf, x1=x1, x2=x2, force_positive=False)#self.polynomial=="bernstein")

    def get_bin_centers_interpolation(self, h, flow=True, cap=False):
        return self.get_bin_centers(h, self.name_x, self.axis_x_min, self.axis_x_max, flow=flow, cap=cap)

    def calculate_fullABCD(self, h, flow=True):
        if len(self.fakerate_integration_axes) > 0:
            logger.warning(f"Binned fake estimation is performed but fakerate integration axes {self.fakerate_integration_axes} are set, the bin-by-bin stat uncertainties are not correct along this axis.")
        if not self.integrate_x:
            logger.warning(f"Binned fake estimation is performed but ABCD x-axis is not integrated, the bin-by-bin stat uncertainties are not correct along this axis.")

        frf, frf_var = self.compute_fakeratefactor(h)
        c_scf, c_scf_var = self.compute_shapecorrection(h, apply=True)

        d = frf * c_scf

        if self.integrate_x:
            idx_x = [n for n in h.axes.name if n != self.name_y].index(self.name_x)
            d = d.sum(axis=idx_x)

            if h.storage_type == hist.storage.Weight:
                # the fakerate factor is independent of the abcd x-axis and treated fully correlated
                dvar = c_scf * frf_var**0.5
                dvar = np.moveaxis(dvar, idx_x, -1)
                dvar = np.einsum("...ni,...nj->...n", dvar, dvar)
                dvar += (frf**2 * c_scf_var).sum(axis=idx_x)

        elif h.storage_type == hist.storage.Weight:
            dvar = c_scf**2 * frf_var + frf**2 * c_scf_var

        return d, dvar
