# =============================================================================
# REDFIT MATHEMATICAL ENGINE (Schulz & Mudelsee, 2002)
# Exact Python port of redfit38e.f90:
# =============================================================================
import math
import warnings
import os as _os
import numpy as np
from scipy.special import gammainc
from scipy.stats import norm as _norm

# -----------------------------------------------------------------------------
# 0. Numerical Recipes RNG: ran1 + gasdev
#    Exact Python port of the Fortran NR routines included in redfit38e.f90
# -----------------------------------------------------------------------------
class NRRandom:
    """
    Exact Python port of Numerical Recipes ran1 (Park-Miller LCG with
    Bays-Durham shuffle, NTAB=32) and gasdev (Box-Muller transform).

    Matches GASDEV(idum) / RAN1(idum) in redfit38e Fortran source exactly.
    Reference: Press et al., Numerical Recipes in Fortran 77, 2nd ed., Â§7.1

    Convention (matching Fortran):
        idum must be NEGATIVE on construction â†’ triggers ran1 initialisation
        idum = -abs(seed)   (Fortran line: idum = -abs(idumini(1)))

    Usage:
        rng = NRRandom(idum)   # idum negative
        x   = rng.gasdev()     # one standard Gaussian deviate
    """
    # ran1 / Park-Miller parameters
    IA   = 16807
    IM   = 2147483647          # 2^31 âˆ’ 1 (Mersenne prime)
    AM   = 1.0 / 2147483647
    IQ   = 127773
    IR   = 2836
    NTAB = 32
    NDIV = 1 + (2147483647 - 1) // 32    # = 67108864
    EPS  = 1.2e-7
    RNMX = 1.0 - 1.2e-7

    def __init__(self, idum):
        """Initialise.  Pass a NEGATIVE integer (Fortran convention)."""
        self.idum  = int(idum)
        self.iv    = [0] * self.NTAB
        self.iy    = 0
        # gasdev saved state  (Fortran SAVE ISET, GSET)
        self._iset = 0
        self._gset = 0.0

    def _ran1(self):
        """NR ran1: Park-Miller LCG with Bays-Durham shuffle.
        Integer arithmetic is exact in Python (arbitrary precision), so the
        Schrage-method computation IA*(idum-k*IQ) âˆ’ IR*k is safe.
        """
        # â”€â”€ Initialise when idum â‰¤ 0 or iy == 0  (matches Fortran IF branch)
        if self.idum <= 0 or self.iy == 0:
            self.idum = max(-self.idum, 1)     # idum = max(âˆ’idum, 1)
            # Warm up shuffle table: NTAB+8 = 40 iterations  (j=40 â€¦ 1)
            for j in range(self.NTAB + 8, 0, -1):
                k = self.idum // self.IQ
                self.idum = (self.IA * (self.idum - k * self.IQ)
                             - self.IR * k)
                if self.idum < 0:
                    self.idum += self.IM
                if j <= self.NTAB:             # store into iv(j), 1-indexed
                    self.iv[j - 1] = self.idum
            self.iy = self.iv[0]               # iy = iv(1)

        # â”€â”€ Standard generation step (Schrage's method)
        k = self.idum // self.IQ
        self.idum = (self.IA * (self.idum - k * self.IQ) - self.IR * k)
        if self.idum < 0:
            self.idum += self.IM
        j       = self.iy // self.NDIV          # 0-indexed shuffle slot
        self.iy = self.iv[j]
        self.iv[j] = self.idum
        return min(self.AM * self.iy, self.RNMX)

    def gasdev(self):
        """NR gasdev: Box-Muller standard Gaussian deviate using ran1.
        Generates deviates in pairs; caches the second in _gset / _iset.
        Matches Fortran GASDEV(idum) algorithm and saved-state behaviour
        (SAVE ISET, GSET in the Fortran F77 source).
        """
        if self._iset == 0:
            # Rejection-sample inside the unit circle
            while True:
                v1 = 2.0 * self._ran1() - 1.0
                v2 = 2.0 * self._ran1() - 1.0
                r  = v1 * v1 + v2 * v2
                if r < 1.0 and r != 0.0:
                    break
            fac        = math.sqrt(-2.0 * math.log(r) / r)
            self._gset = v1 * fac          # cache partner deviate (GSET)
            self._iset = 1                 # flag that partner is ready (ISET=1)
            return v2 * fac
        else:
            self._iset = 0
            return self._gset             # return cached partner


# -----------------------------------------------------------------------------
# 1. Core Spectral Engine
#    rmtrend, winwgt, trig_arrays, ftfix, spectr, winbw
# -----------------------------------------------------------------------------
def rmtrend(t, x):
    """Remove linear trend by ordinary least squares.
    Matches Fortran SR rmtrend (lines 1414-1458 of redfit38e.f90) exactly.
    """
    t      = np.asarray(t, dtype=float)
    x      = np.asarray(x, dtype=float)
    n      = len(t)
    sx     = np.sum(t)
    sy     = np.sum(x)
    sxoss  = sx / n
    z      = t - sxoss
    st2    = np.sum(z * z)
    b      = np.sum(z * x) / st2
    a      = (sy - sx * b) / n
    return x - (a + b * t)


def winwgt(t, iwin):
    """Normalised window weights.
    Matches Fortran SR winwgt (lines 1328-1384) exactly for all five window
    types.  Fortran coefficients for Blackman-Harris preserved verbatim.

    iwin: 0=Rectangular, 1=Welch, 2=Hanning, 3=Triangular, 4=Blackman-Harris
    """
    nseg = len(t)
    rnp  = float(nseg)
    fac1 = rnp / 2.0 - 0.5
    fac2 = 1.0 / (rnp / 2.0 + 0.5)
    fac3 = rnp - 1.0
    fac4 = 2.0 * np.pi / (rnp - 1.0)
    tlen = t[-1] - t[0]
    jeff = rnp * (t - t[0]) / tlen          # fractional position (Fortran jeff)

    if   iwin == 0:
        w = np.ones(nseg)
    elif iwin == 1:
        w = 1.0 - ((jeff - fac1) * fac2) ** 2.0
    elif iwin == 2:
        w = 0.5 * (1.0 - np.cos(2.0 * np.pi * jeff / fac3))
    elif iwin == 3:
        w = 1.0 - np.abs((jeff - fac1) * fac2)
    elif iwin == 4:
        w = (0.4243801
             - 0.4973406 * np.cos(fac4 * jeff)
             + 0.0782793 * np.cos(fac4 * 2.0 * jeff))
    else:
        raise ValueError("iwin must be 0-4")

    sumw2 = np.sum(w * w)
    scal  = np.sqrt(rnp / sumw2)
    return w * scal


def trig_arrays(t, wz, nfreq):
    """Precompute Scargle tau-phase and cos/sin arrays.
    Matches Fortran SR trig (lines 1235-1300).
    Uses double-precision arithmetic for trig arguments, matching the v3.8e
    bugfix ('double :: arg' at line 1250 of redfit38e.f90).
    """
    nseg = len(t)
    tcos = np.zeros((nfreq, nseg))
    tsin = np.zeros((nfreq, nseg))
    wtau = np.zeros(nfreq)
    wrun = float(wz)                        # explicit float64 = Fortran double

    for ii in range(1, nfreq):              # ii=2..nfreq in Fortran (1-indexed)
        arg  = 2.0 * wrun * t               # double precision (Fortran: double arg)
        csum = np.sum(np.cos(arg))
        ssum = np.sum(np.sin(arg))
        # tau-phase estimation (Fortran lines 1279-1283)
        if abs(ssum) > 1.0e-4 or abs(csum) > 1.0e-4:
            watan = np.arctan2(ssum, csum)
        else:
            sumtc = np.sum(t * np.cos(arg))
            sumts = np.sum(t * np.sin(arg))
            watan = np.arctan2(-sumtc, sumts)
        wt        = 0.5 * watan
        wtau[ii]  = wt
        arg2      = wrun * t - wt
        tcos[ii]  = np.cos(arg2)
        tsin[ii]  = np.sin(arg2)
        wrun     += wz
    return tcos, tsin, wtau


def ftfix(x, t, tcos, tsin, wtau, nfreq):
    """Scargle (1989) Lomb-Scargle FT for unevenly spaced data.
    Matches Fortran SR ftfix (lines 1144-1231) exactly.
    Fortran constants: si=1.0, tzero=0.0  (lines 1012-1013)
    â†’ const2 = si*const1 = const1;  phase = wtnew (tzero term vanishes).
    """
    nn     = len(x)
    const1 = 1.0 / math.sqrt(2.0)          # Fortran parameter const1 = 1/âˆš2
    tol1   = 1.0e-4                        # Fortran parameter tol1
    tol2   = 1.0e-8                        # Fortran parameter tol2
    ftrx   = np.zeros(nfreq)
    ftix   = np.zeros(nfreq)

    sumx      = np.sum(x)
    # f=0 (DC) term  (Fortran: ftrx(1) = sumx / sqrt(fnn))
    ftrx[0]   = sumx / math.sqrt(float(nn))
    ftix[0]   = 0.0

    for ii in range(1, nfreq):
        tc    = tcos[ii]
        ts    = tsin[ii]
        cross = np.sum(t * tc * ts)
        scos2 = np.sum(tc * tc)
        ssin2 = np.sum(ts * ts)
        sumr  = np.sum(x * tc)
        sumi  = np.sum(x * ts)

        ftrd = const1 * sumr / math.sqrt(scos2)
        if ssin2 <= tol1:
            # Fortran: ftid = const2*sumx/sqrt(fnn)  where const2=si*const1=const1
            ftid = const1 * sumx / math.sqrt(float(nn))
            if abs(cross) > tol2:
                ftid = 0.0
        else:
            ftid = const1 * sumi / math.sqrt(ssin2)

        # phase = wtnew - wrun*tzero = wtau[ii]  (tzero=0 â†’ second term vanishes)
        phase    = wtau[ii]
        c        = math.cos(phase)
        s        = math.sin(phase)
        # Expand cmplx(ftrd,ftid)*cexp(i*phase) = (ftrd+i*ftid)*(c+i*s)
        ftrx[ii] = ftrd * c - ftid * s
        ftix[ii] = ftrd * s + ftid * c

    return ftrx, ftix


def spectr(t, x, ofac, hifac, n50, iwin, cache=None):
    """WOSA Lomb-Scargle autospectrum.
    Matches Fortran SR spectr (lines 1005-1111) exactly.
    Pass cache=None on the first call (ini=.true. path); pass the returned
    cache on subsequent calls to reuse trig arrays (ini=.false. path).
    Returns (freq, gxx, cache).
    """
    npts  = len(t)
    nseg  = int(2 * npts / (n50 + 1))             # points per segment
    avgdt = (t[-1] - t[0]) / float(npts - 1)      # average sampling interval
    tp    = avgdt * nseg                           # average segment period
    df    = 1.0 / (ofac * tp)                      # frequency spacing
    wz    = 2.0 * np.pi * df
    fnyq  = hifac * 1.0 / (2.0 * avgdt)           # average Nyquist frequency
    nfreq = int(fnyq / df) + 1                    # f[0]=0; f[nfreq-1]â‰ˆfNyq
    freq  = np.arange(nfreq) * df

    gxx         = np.zeros(nfreq)
    build_cache = (cache is None)
    if build_cache:
        cache = {'ww': [], 'tcos': [], 'tsin': [], 'wtau': []}

    for i in range(n50):
        # istart = (i-1)*nseg/2 in Fortran (1-based) = i*nseg//2 in Python (0-based)
        istart = (i * nseg) // 2
        twk    = t[istart:istart + nseg].copy()
        xwk    = x[istart:istart + nseg].copy()

        xwk = rmtrend(twk, xwk)                   # detrend (Fortran line 1077)

        if build_cache:
            ww         = winwgt(twk, iwin)
            tc, ts, wt = trig_arrays(twk, wz, nfreq)
            cache['ww'].append(ww)
            cache['tcos'].append(tc)
            cache['tsin'].append(ts)
            cache['wtau'].append(wt)
        else:
            ww = cache['ww'][i]
            tc = cache['tcos'][i]
            ts = cache['tsin'][i]
            wt = cache['wtau'][i]

        xwk        = ww * xwk                      # apply window (Fortran line 1082)
        ftrx, ftix = ftfix(xwk, twk, tc, ts, wt, nfreq)
        gxx       += ftrx ** 2 + ftix ** 2         # accumulate power

    # Scale (Fortran lines 1101-1105: scal = 2/(n50*nseg*df*ofac))
    scal = 2.0 / (n50 * nseg * df * ofac)
    gxx *= scal

    return freq, gxx, cache


def winbw(iwin, df, ofac):
    """6-dB bandwidth.
    Matches Fortran function winbw (lines 1388-1410).
    bw table (Fortran line 1406): Rect=1.21, Welch=1.59, Hanning=2.00,
                                  Triangular=1.78, Blackman-Harris=2.26
    """
    bw = [1.21, 1.59, 2.00, 1.78, 2.26]
    return df * ofac * bw[iwin]


# -----------------------------------------------------------------------------
# 2. Tau Estimation
#    brent_search, _ls_function, minls, rhoest, tauest, gettau
# -----------------------------------------------------------------------------
def brent_search(ax, bx, cx, func, tol=3.0e-8, itmax=100):
    """Brent's method line minimisation.
    Matches Fortran function brent (lines 1773-1862) exactly.
    """
    cgold = 0.3819660          # golden-ratio complement
    zeps  = 1.0e-18
    a     = min(ax, cx)
    b     = max(ax, cx)
    v = w = x = bx
    e  = 0.0
    d  = 0.0
    fx = fw = fv = func(x)

    for _ in range(itmax):
        xm   = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            break

        golden = True
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q      = abs(q)
            etemp  = e
            e      = d
            if not (abs(p) >= abs(0.5 * q * etemp)
                    or p <= q * (a - x)
                    or p >= q * (b - x)):
                d      = p / q
                u      = x + d
                if (u - a) < tol2 or (b - u) < tol2:
                    d = math.copysign(tol1, xm - x)
                golden = False

        if golden:
            e = (a - x) if x >= xm else (b - x)
            d = cgold * e

        u  = (x + d) if abs(d) >= tol1 else (x + math.copysign(tol1, d))
        fu = func(u)

        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu

    return x, fx


def _ls_function(a, t, x):
    """Least-squares misfit for AR(1) fit.
    Matches Fortran function ls (lines 1868-1879) exactly.
    """
    return np.sum(
        (x[1:] - x[:-1] * np.sign(a) * np.abs(a) ** (t[1:] - t[:-1])) ** 2.0
    )


def minls(t, x):
    """Minimise LS function over three brackets to guard against local minima.
    Matches Fortran SR minls (lines 1885-1922) exactly.
    """
    a_ar1 = 0.367879441     # 1/e
    tol   = 3.0e-8
    tol2  = 1.0e-6
    f     = lambda a: _ls_function(a, t, x)

    a1, dum1 = brent_search(-2.0,  a_ar1,                  2.0,   f, tol)
    a2, dum2 = brent_search(a_ar1, 0.5 * (a_ar1 + 1.0),   2.0,   f, tol)
    a3, dum3 = brent_search(-2.0,  0.5 * (a_ar1 - 1.0),   a_ar1, f, tol)

    mult = 0
    if ((abs(a2 - a1) > tol2 and abs(a2 - a_ar1) > tol2) or
            (abs(a3 - a1) > tol2 and abs(a3 - a_ar1) > tol2)):
        mult = 1

    dum_min = min(dum1, dum2, dum3)
    if   dum_min == dum2:
        amin = a2
    elif dum_min == dum3:
        amin = a3
    else:
        amin = a1
    return amin, mult


def rhoest(x):
    """Yule-Walker autocorrelation coefficient for equidistant data.
    Matches Fortran SR rhoest (lines 1928-1946): sum(x[i]*x[i-1]) / sum(x[i]Â²)
    """
    return np.sum(x[1:] * x[:-1]) / np.sum(x[1:] ** 2.0)


def tauest(t, x):
    """Estimate persistence time tau from an unevenly spaced segment.
    Matches Fortran SR tauest (lines 1684-1767) exactly:
      - Reverses the time axis (geological ages â†’ physical time direction).
      - Normalises x with ddof=1 variance (matching NR avevar).
      - Calls minls / brent for LS minimisation.
    """
    n     = len(t)
    # Reverse direction: geological ages â†’ increasing physical time (Fortran lines 1712-1715)
    tscal = -t[::-1].copy()
    xscal =  x[::-1].copy()

    # Normalise x to unit variance using ddof=1 (matches NR avevar)
    var   = np.var(xscal, ddof=1)
    fac   = math.sqrt(var)
    xscal = xscal / fac

    dt  = (tscal[-1] - tscal[0]) / float(n - 1)
    rho = rhoest(xscal)
    if rho <= 0.0:
        rho = 0.05
        warnings.warn("tauest: rho estimate <= 0, reset to 0.05")
    elif rho > 1.0:
        rho = 0.95
        warnings.warn("tauest: rho estimate > 1, reset to 0.95")

    scalt = -math.log(rho) / dt
    tscal = tscal * scalt

    amin, mult = minls(tscal, xscal)
    if mult == 1:
        warnings.warn("tauest: LS function has more than one minimum")
    if amin <= 0.0 or amin >= 1.0:
        warnings.warn(f"tauest: a_min out of (0,1) bounds (a_min={amin:.6f})")

    tau    = -1.0 / (scalt * math.log(amin))
    rhoavg = math.exp(-dt / tau)
    return tau, rhoavg


def gettau(t, x, n50):
    """Average tau over n50 WOSA segments with Kendall-Stuart bias correction.
    Matches Fortran SR gettau (lines 1601-1662) exactly.
    Bias correction formula (Fortran line 1643):
        rho = (rho*(nseg-1) + 1) / (nseg-4)
    avgdt is computed over the full time series (Fortran lines 1653-1654).
    """
    npts   = len(t)
    nseg   = int(2 * npts / (n50 + 1))
    rhosum = 0.0
    for i in range(n50):
        istart = (i * nseg) // 2
        twk    = t[istart:istart + nseg].copy()
        xwk    = x[istart:istart + nseg].copy()
        xwk    = rmtrend(twk, xwk)
        _, rho_i = tauest(twk, xwk)
        # Kendall & Stuart (1967) bias correction (Fortran line 1643)
        rho_i  = (rho_i * (float(nseg) - 1.0) + 1.0) / (float(nseg) - 4.0)
        rhosum += rho_i
    rho   = rhosum / float(n50)
    avgdt = np.sum(np.diff(t)) / float(npts - 1)   # = mean(diff(t))
    tau   = -avgdt / math.log(rho)
    return tau, rho


# -----------------------------------------------------------------------------
# 3. AR(1) Surrogate Generation â€” uses NR ran1 + gasdev
# -----------------------------------------------------------------------------
def makear1(t, tau, nrrng):
    """Generate one AR(1) surrogate at the given (uneven) sample times.
    Matches Fortran SR makear1 (lines 1115-1140) exactly.
    Uses NRRandom.gasdev() to replicate GASDEV(idum)/RAN1(idum) from Fortran.

    Parameters
    ----------
    t     : array of sample times (ages)
    tau   : persistence time
    nrrng : NRRandom instance (carries the mutable RNG state)
    """
    n      = len(t)
    red    = np.empty(n)
    red[0] = nrrng.gasdev()                          # Fortran: red(1) = gasdev(idum)
    for i in range(1, n):
        dt     = t[i] - t[i - 1]
        sigma  = math.sqrt(1.0 - math.exp(-2.0 * dt / tau))
        red[i] = math.exp(-dt / tau) * red[i - 1] + sigma * nrrng.gasdev()
    return red


# -----------------------------------------------------------------------------
# 4. Statistical Helpers
#    getdof, _getz, getchi2
# -----------------------------------------------------------------------------
def getdof(iwin, n50):
    """Effective degrees of freedom (Harris, 1978).
    Matches Fortran SR getdof (lines 1462-1481).
    c50 table (Fortran line 1472): [0.500, 0.344, 0.167, 0.250, 0.096]
    """
    c50   = [0.500, 0.344, 0.167, 0.250, 0.096]
    rn    = float(n50)
    c2    = 2.0 * c50[iwin] ** 2.0
    denom = 1.0 + c2 - c2 / rn
    neff  = rn / denom
    return 2.0 * neff


def _getz(alpha):
    """Normal-distribution percentile z such that P[Z <= z] = alpha.
    Matches Fortran function getz (lines 1537-1597).
    Fortran uses an erfcc Chebyshev approximation; scipy.stats.norm.ppf
    converges to the same mathematical value.
    """
    return float(_norm.ppf(alpha))


def getchi2(dof, alpha):
    """Chi-square value such that P[chi2 >= result] = alpha (upper tail).
    Matches Fortran function getchi2 (lines 1485-1533) algorithm exactly:
      - dof > 30: Sachs (1984) Eq. 1.132 cubic approximation
                  chi2 = dof*(1 - 2/(9*dof) + za*sqrt(2/(9*dof)))^3
                  where za = -getz(alpha)  [sign change as in Fortran comment]
      - dof <= 30: bisection on  ac = 1 - gammp(dof/2, chi2/2)
                   scipy.special.gammainc(a,x) = P(a,x) = gammp(a,x) in NR.
    """
    tol   = 1.0e-3
    itmax = 100
    if dof > 30.0:
        # Sachs (1984) approximation (Fortran lines 1500-1504)
        za   = -_getz(alpha)              # sign change noted in Fortran comment
        x    = 2.0 / 9.0 / dof
        chi2 = dof * (1.0 - x + za * math.sqrt(x)) ** 3.0
    else:
        # Bisection on the incomplete-gamma CDF (Fortran lines 1506-1529)
        lm  = 0.0
        rm  = 1000.0
        eps = (1.0 - alpha) * tol if alpha > 0.5 else alpha * tol
        for _ in range(itmax):
            chi2 = 0.5 * (lm + rm)
            # gammp(a,x) = regularised lower incomplete gamma P(a,x)
            ac   = 1.0 - float(gammainc(0.5 * dof, 0.5 * chi2))
            if abs(ac - alpha) <= eps:
                break
            if ac > alpha:
                lm = chi2
            else:
                rm = chi2
    return chi2


print("REDFIT Mathematical Engine v3.8e (Python) â€” successfully loaded into memory.")
