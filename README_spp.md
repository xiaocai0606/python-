# GPS pseudorange single point positioning

This workspace contains a complete Python implementation for reading RINEX 2.x
observation/navigation files and solving GPS pseudorange single point positioning.

## Files

- `spp.py`: main program. It includes RINEX O/N readers, broadcast ephemeris
  satellite position calculation, Earth rotation correction, and iterative least
  squares positioning.
- `wuhn1660.12o`: copied test observation file.
- `brdc1660.12n`: copied test navigation file.
- `spp_first_epoch.csv`: first-epoch test result.
- `spp_20_epochs.csv`: first 20 epoch test results.

## Run

```powershell
python spp.py wuhn1660.12o brdc1660.12n --epoch 0
```

Solve from zero initial coordinates, matching the teaching flow:

```powershell
python spp.py wuhn1660.12o brdc1660.12n --epoch 0 --initial zero
```

Solve multiple epochs and save CSV:

```powershell
python spp.py wuhn1660.12o brdc1660.12n --all --max-epochs 20 --csv spp_20_epochs.csv
```

## Notes

- Pseudorange priority defaults to `C1,P1,P2`; change it with `--obs`.
- The program uses GPS broadcast ephemerides and applies satellite clock,
  relativistic clock, TGD, and Earth rotation corrections.
- Ionospheric and tropospheric corrections are not applied because the provided
  flow document focuses on the basic pseudorange SPP least-squares process.
