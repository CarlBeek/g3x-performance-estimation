# G3X Performance Estimation

Build a POH-style performance table (cruise, climb and descent vs. density altitude) for a
light aircraft from nothing more than its Garmin G3X flight logs.

Rather than flying dedicated test points, the script mines ordinary flights for stretches of
steady, wings-level, clean-configuration flight, fits a simple energy/drag-polar model to them,
and then evaluates that model at every 1,000 ft of density altitude. It was written for an RV-6,
but nothing in it is type-specific apart from `VNE_KTAS`.

## The model

All power is expressed in the G3X's own *Engine Power %* units, so rated horsepower, propeller
efficiency and weight are absorbed into the fitted constants:

```
pwr = A·σ·V³ + B·nz² / (σ·V) + C·Ps
```

| Symbol | Meaning |
| --- | --- |
| `V`  | true airspeed / 100 kt |
| `σ`  | density ratio, from density altitude |
| `nz` | normal load factor |
| `Ps` | specific excess power in fpm: `VS + (V/g)·dV/dt` |
| `A`  | parasite drag term |
| `B`  | induced drag term |
| `C`  | % power per 1,000 fpm of climb |

Around this polar the script also fits:

- **Power available** – full-throttle climb power, linear in σ (Gagg-Ferrar form), plus the
  absolute WOT ceiling from the 99th-percentile power in each altitude band.
- **Pilot technique** – the climb IAS schedule, typical descent IAS and rate, and how close to
  WOT cruise is flown at altitude. The table therefore describes how the aircraft is *actually
  flown*, not a theoretical optimum.
- **Fuel flow** – linear in power, separately for climb, cruise and descent, which captures how
  the mixture is leaned in each phase.

Uncertainty bands on cruise speed and rate of climb (10th–90th percentile) come from a
bootstrap that resamples whole flights.

### Steady-state filter

Samples are taken from centred 30 s windows (decimated to one every 10 s) that satisfy:

- IAS > 75 kt, flaps up, |roll| < 8°, > 400 ft AGL, RPM > 1900
- standard deviation of power < 2.5 %, of vertical speed < 250 fpm, of IAS < 3 kt

## Usage

```sh
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# put G3X log CSVs (log_YYYYMMDD_HHMMSS_XXXX.csv) in ./data
python perf_model.py
```

Options:

| Flag | Default | |
| --- | --- | --- |
| `--data-dir` | `./data` | directory of G3X log CSVs |
| `--fdr` | any `*.zip` in the project root | optional G5/FDR export(s); only data older than the first G3X log is used |
| `--out` | `./performance_table.csv` | output path |
| `--bootstrap` | `200` | number of bootstrap resamples |

The script prints the fitted constants and writes the table to CSV.

## Output

`performance_table.csv` holds one row per 1,000 ft of density altitude from 0 to 23,000 ft:

| Column | |
| --- | --- |
| `cruise_pwr`, `cruise_KTAS`, `cruise_GPH` | level cruise at the pilot's usual power setting (≤ 75 %) |
| `climb_KIAS`, `climb_KTAS`, `climb_GPH`, `climb_FPM` | full-power climb at the pilot's usual climb speed |
| `desc_KIAS`, `desc_KTAS`, `desc_GPH`, `desc_FPM` | cruise descent at the pilot's usual speed and rate, capped at V<sub>NE</sub> |
| `*_lo`, `*_hi` | 10th / 90th percentile bootstrap bounds for cruise KTAS and climb FPM |

Example output from ~55 flights in an RV-6 (selected rows):

| DA (ft) | Cruise % | Cruise KTAS | Cruise GPH | Climb KIAS | Climb FPM | Climb GPH | Descent KTAS | Descent FPM |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0      | 75.0 | 153 | 9.0 | 121 | 1338 | 12.7 | 137 | −495 |
| 5,000  | 71.7 | 158 | 8.7 | 119 | 1004 | 11.1 | 148 | −495 |
| 10,000 | 63.8 | 157 | 8.0 | 117 |  690 |  9.6 | 159 | −495 |
| 15,000 | 56.7 | 156 | 7.3 | 115 |  394 |  8.4 | 173 | −553 |
| 20,000 | 50.4 | 153 | 6.8 | 113 |  110 |  7.3 | 183 | −847 |

## Caveats

This is an empirical estimate from unplanned flight data – weight, CG, temperature and
technique all vary from flight to flight and none of it is controlled. Values far outside the
altitudes you actually fly are extrapolations. It is **not** a substitute for your aircraft's
POH or your own flight testing; don't use it for flight planning where margins matter.
