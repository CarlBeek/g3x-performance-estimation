"""Estimate aircraft performance tables from Garmin G3X CSV logs (data/*.csv), optionally
supplemented by an older G5/FDR export (*.zip containing high-rate "fdr_*.csv").

Model (all power in G3X "Engine Power %" units, so rated HP / prop efficiency / weight
are absorbed into the fitted constants):

    pwr = A*sig*V^3 + B*nz^2/(sig*V) + C*Ps

    V   = TAS / 100 kt, sig = density ratio from density altitude
    Ps  = specific excess power in fpm (VS + V/g * dV/dt)
    A   = parasite drag term, B = induced drag term, C = % power per 1000 fpm of climb

Power available (WOT, climb RPM) is fit linearly in sigma (Gagg-Ferrar form).
Fuel flow is fit linearly vs power, separately for climb / cruise / descent, which captures
how the pilot actually leans in each phase.

Usage: python perf_model.py [--data-dir data] [--fdr logs.zip ...] [--out performance_table.csv]
"""
import argparse
import glob
import os
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares

VNE_KTAS = 183.0
ROOT = os.path.dirname(os.path.abspath(__file__))

# G3X column name -> short name
COLS = {
    "Pressure Altitude (ft)": "palt", "Vertical Speed (ft/min)": "vs", "Indicated Airspeed (kt)": "ias",
    "True Airspeed (kt)": "tas", "Roll (deg)": "roll", "Density Altitude (ft)": "dalt",
    "Engine Power (%)": "pwr", "RPM": "rpm", "Fuel Flow (gal/hour)": "ff", "Flap Position": "flaps",
    "Height Above Ground (ft)": "agl", "Normal Acceleration (G)": "nz",
}

# G5/FDR export column name -> short name
FDR_COLS = {
    "Pressure Altitude (ft)": "palt", "Vertical Speed (fpm)": "vs", "Indicated Airspeed (kt)": "ias",
    "True Airspeed (kt)": "tas", "Roll (deg)": "roll", "Density Altitude (ft)": "dalt",
    "Engine Power (%)": "pwr", "RPM": "rpm", "Fuel Flow (gal/hr)": "ff", "Flap Position": "flaps",
    "Normal Accel (G)": "nz", "GPS Altitude (WGS84 ft)": "galt",
}


def sigma(da):
    return ((288.15 - 0.0019812 * np.asarray(da, float)) / 288.15) ** 4.2559


def load_g3x(data_dir):
    """One 1 Hz frame per G3X log file, plus the earliest UTC time covered by these logs."""
    frames, first_utc = [], []
    for f in sorted(glob.glob(os.path.join(data_dir, "*.csv"))):
        d = pd.read_csv(f, skiprows=[0, 2], low_memory=False, on_bad_lines="skip")
        d.columns = [c.strip() for c in d.columns]
        utc = pd.to_datetime(d["Date (yyyy-mm-dd)"] + " " + d["UTC Time (hh:mm:ss)"], errors="coerce")
        first_utc.append(utc.min())
        d = d[[c for c in COLS if c in d.columns]].rename(columns=COLS).apply(pd.to_numeric, errors="coerce")
        d["flight"] = os.path.basename(f)
        frames.append(d)
    return frames, min(t for t in first_utc if pd.notna(t))


def load_fdr(zip_path, before_utc):
    """G5/FDR export (~7-8 Hz, UTC timestamps): resample to 1 Hz, split into flights at power
    cycles, keep only data before `before_utc` (later flights are already in the G3X logs).
    It has no AGL column, so AGL is approximated as GPS altitude above the flight's first fix."""
    frames = []
    with zipfile.ZipFile(zip_path) as z:
        for name in z.namelist():
            if not name.endswith(".csv"):
                continue
            with z.open(name) as fh:
                d = pd.read_csv(fh, skiprows=[0], low_memory=False)
            d.columns = [c.strip() for c in d.columns]
            ts = pd.to_datetime(d["UTC Date (yyyy-mm-dd)"] + " " + d["UTC Time (hh:mm:ss)"], errors="coerce")
            ms = pd.to_numeric(d["Power Timestamp (msec)"], errors="coerce")
            d = d[list(FDR_COLS)].rename(columns=FDR_COLS).apply(pd.to_numeric, errors="coerce")
            d["seg"] = ((ms.diff() < 0) | (ms.diff() > 60000)).cumsum()
            d["sec"], d["ts"] = ms // 1000, ts
            d = d[d.ts < before_utc]
            for seg, x in d.groupby("seg"):
                x = x.groupby("sec").mean(numeric_only=True).reset_index(drop=True)
                if (x.ias > 60).sum() < 300:
                    continue
                x["agl"] = x.galt - x.galt.dropna().iloc[0]
                x["nz"] += 1.0  # FDR logs normal accel relative to 1 g (reads ~0 in level flight)
                x["flight"] = f"{name}#{seg}"
                frames.append(x[list(COLS.values()) + ["flight"]])
    return frames


def load(data_dir, fdr_zips):
    frames, g3x_start = load_g3x(data_dir)
    for z in fdr_zips:
        frames += load_fdr(z, g3x_start)
    return frames


def steady_segments(frames, win=31):
    """30 s rolling windows of steady, clean-configuration, wings-level flight (decimated to 10 s)."""
    out = []
    for d in frames:
        if len(d) < 200:
            continue
        d = d.copy()
        roll = lambda c: d[c].rolling(win, center=True)
        for c in ["ias", "tas", "vs", "pwr", "ff", "dalt", "rpm", "nz"]:
            d[c + "_m"] = roll(c).mean()
        d["ias_sd"], d["vs_sd"], d["pwr_sd"] = roll("ias").std(), roll("vs").std(), roll("pwr").std()
        d["roll_max"] = d.roll.abs().rolling(win, center=True).max()
        d["flap_max"] = roll("flaps").max()
        d["agl_min"] = roll("agl").min()
        d["dtas"] = (d.tas.shift(-15) - d.tas.shift(15)) / 30.0  # kt/s
        out.append(d.iloc[::10])
    s = pd.concat(out)
    s = s[(s.ias_m > 75) & (s.flap_max == 0) & (s.roll_max < 8) & (s.agl_min > 400) & (s.rpm_m > 1900)
          & (s.pwr_sd < 2.5) & (s.vs_sd < 250) & (s.ias_sd < 3)].copy()
    s["sig"] = sigma(s.dalt_m)
    s["Ps"] = s.vs_m + s.tas_m * 1.6878 * s.dtas * 1.6878 / 32.174 * 60
    return s.dropna(subset=["Ps", "pwr_m", "tas_m", "sig", "nz_m", "ff_m"])


def fit_polar(s):
    V = s.tas_m.values / 100
    X = np.c_[s.sig * V**3, s.nz_m**2 / (s.sig * V), s.Ps / 1000]
    return least_squares(lambda p: X @ p - s.pwr_m.values, [10, 20, 20], loss="soft_l1", f_scale=3).x


class Model:
    def __init__(self, s, raw):
        self.A, self.B, self.C = fit_polar(s)
        # climb power available: steady climbs above 6000 DA (below that the pilot often limits power)
        cl = s[(s.vs_m > 300) & (s.dalt_m > 6000)]
        self.pclimb = np.polyfit(cl.sig, cl.pwr_m, 1)
        # absolute WOT power (99th pct per 1000 ft DA band) - used to cap descent power
        a = raw[(raw.ias > 70) & (raw.agl > 300) & (raw.flaps == 0)].dropna(subset=["dalt", "pwr"])
        g = a.groupby(pd.cut(a.dalt, np.arange(1000, 17000, 1000)), observed=True)
        m = g.agg(n=("pwr", "size"), da=("dalt", "median"), p99=("pwr", lambda x: x.quantile(0.99)))
        m = m[m.n > 500]
        self.pmax = np.polyfit(sigma(m.da), m.p99, 1)
        # cruise: pilot flies ~75% down low, near-WOT at cruise RPM up high
        lv = s[(s.vs_m.abs() < 100)]
        hi = lv[lv.dalt_m > 9000]
        self.cruise_frac = float(np.median(hi.pwr_m / np.polyval(self.pmax, hi.sig)))
        # climb IAS schedule the pilot actually flies
        self.climb_ias = np.polyfit(cl.dalt_m, cl.ias_m, 1)
        de = s[s.vs_m < -300]
        self.desc_ias = float(de.ias_m.median())
        self.desc_vs = float(de.vs_m.median())
        # fuel flow vs power per phase
        self.ff = {k: np.polyfit(d.pwr_m, d.ff_m, 1) for k, d in
                   {"cruise": lv, "climb": s[s.vs_m > 300], "descent": de}.items()}

    def preq(self, tas, sig, ps=0.0):
        V = tas / 100
        return self.A * sig * V**3 + self.B / (sig * V) + self.C * ps / 1000

    def cruise_tas(self, pwr, sig):
        vmp = (self.B / (3 * self.A)) ** 0.25 * 100 / np.sqrt(sig)  # min-power TAS
        if self.preq(vmp, sig) >= pwr:
            return np.nan
        return brentq(lambda v: self.preq(v, sig) - pwr, vmp, 400)

    def row(self, da):
        sg = float(sigma(da))
        pmax, pcl = np.polyval(self.pmax, sg), np.polyval(self.pclimb, sg)
        # cruise
        pc = min(75.0, self.cruise_frac * pmax)
        tas_c = self.cruise_tas(pc, sg)
        if tas_c > VNE_KTAS:
            tas_c, pc = VNE_KTAS, self.preq(VNE_KTAS, sg)
        # climb at pilot's IAS schedule, full power
        ias_cl = np.polyval(self.climb_ias, da)
        tas_cl = ias_cl / np.sqrt(sg)
        roc = max(0.0, (pcl - self.preq(tas_cl, sg)) / self.C * 1000)
        # descent at pilot's typical IAS and rate, capped at VNE; never above cruise power
        # (if holding the rate would need more, the descent steepens instead)
        tas_d = min(self.desc_ias / np.sqrt(sg), VNE_KTAS)
        pd_ = self.preq(tas_d, sg, self.desc_vs)
        vs_d = self.desc_vs
        if pd_ > pc:
            pd_ = pc
            vs_d = (pc - self.preq(tas_d, sg)) / self.C * 1000
        ff = lambda k, p: np.polyval(self.ff[k], p)
        return dict(
            DA=da,
            cruise_pwr=pc, cruise_KTAS=tas_c, cruise_GPH=ff("cruise", pc),
            climb_KIAS=ias_cl, climb_KTAS=tas_cl, climb_GPH=ff("climb", pcl), climb_FPM=roc,
            desc_KIAS=tas_d * np.sqrt(sg), desc_KTAS=tas_d, desc_GPH=ff("descent", pd_), desc_FPM=vs_d,
        )


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data-dir", default=os.path.join(ROOT, "data"), help="directory of G3X log CSVs")
    p.add_argument("--fdr", nargs="*", default=None,
                   help="G5/FDR export zip(s); default: any *.zip in the project root")
    p.add_argument("--out", default=os.path.join(ROOT, "performance_table.csv"), help="output CSV path")
    p.add_argument("--bootstrap", type=int, default=200, help="bootstrap resamples for uncertainty bands")
    return p.parse_args()


def main():
    args = parse_args()
    fdr = sorted(glob.glob(os.path.join(ROOT, "*.zip"))) if args.fdr is None else args.fdr
    raw_frames = load(args.data_dir, fdr)
    raw = pd.concat(raw_frames)
    s = steady_segments(raw_frames)
    m = Model(s, raw)
    print(f"steady samples: {len(s)}  flights: {s.flight.nunique()}")
    print(f"polar A={m.A:.3f} B={m.B:.3f} C={m.C:.3f}  pclimb={m.pclimb.round(1)}  pmax={m.pmax.round(1)}")
    print(f"cruise frac of WOT={m.cruise_frac:.3f}  climb IAS={m.climb_ias.round(5)}  "
          f"descent {m.desc_ias:.0f} KIAS @ {m.desc_vs:.0f} fpm  ff={ {k: v.round(3) for k, v in m.ff.items()} }")

    das = np.arange(0, 23001, 1000)
    tab = pd.DataFrame([m.row(d) for d in das])

    # bootstrap by flight for uncertainty on the headline numbers
    flights = s.flight.unique()
    rng = np.random.default_rng(1)
    boots = []
    for _ in range(args.bootstrap):
        pick = rng.choice(flights, len(flights))
        sb = pd.concat([s[s.flight == f] for f in pick])
        try:
            boots.append(pd.DataFrame([Model(sb, raw).row(d) for d in das]))
        except Exception:  # degenerate resample (e.g. no qualifying climbs) - skip it
            pass
    for c in ["cruise_KTAS", "climb_FPM"]:
        q = np.array([b[c].values for b in boots])
        tab[c + "_lo"], tab[c + "_hi"] = np.nanpercentile(q, 10, 0), np.nanpercentile(q, 90, 0)

    tab = tab.round(1)
    tab.to_csv(args.out, index=False)
    pd.set_option("display.width", 250)
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
