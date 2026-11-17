# Captions - moki-julio-circle-packing

Run: `python figures/make_figures.py .` from the repository root (Python 3.11, matplotlib 3.11; standard library otherwise).
The script re-checks every drawn packing in floating point (inside the container, no overlaps) and stops if one fails; it
also asserts the facts the triangle figure shows (11 shared frame points, 20 contacts, minimum distance).

## 1. `figures/packing_gallery.svg` (hero)

**Point.** Records are real, checkable pictures: six packings, one per container shape, that were still the best known.

**Caption.** Six of the record packings that were still the best known on 2026-10-18, one per container shape, drawn from
the certificates in this repository. In each, equal circles are as large as possible without overlapping each other or
the edge. The percentage is how much larger the radius is than the published record it beat. Most gains are small (here
0.014 % to 0.37 %); the pictures show new arrangements, not new kinds of packing.

**Alt text.** Six panels of tightly packed green circles inside a rectangle, a right triangle, a quarter circle, a half
circle, a regular pentagon and a regular 15-gon, each titled with its number of circles and its radius gain.

**Selection rule (in the script, not by hand).** Containers crc_800, crt, ccq, csc, cpt, cpd; in each, the smallest N
whose record is not in SUPERSEDED.csv and whose `relative_gain` is at least 1e-4.

**Sources.**
| item | file | field |
|---|---|---|
| circle centres and radius | `certificates/{crc_800,crt,ccq,csc,cpt,cpd}.zip` -> `<table>_<N>.txt` | line 1 `r`, then `x y` |
| N, gain, kind | `MANIFEST.csv` | `shelf`, `N`, `relative_gain`, `kind` (all six: NEW ARRANGEMENT) |
| "still best known" | `SUPERSEDED.csv` | (table, N) absent |
| date 2026-10-18 | `README.md` | "4,063 were still the best known on 2026-10-18" |
| container frames | `TECHNICAL.md` | "Same coordinate frame as Packomania's own files" |

Numbers in the panels: rectangle 1 x 0.8, N = 10, gain 3.710952e-03 (+0.371 %); right triangle N = 67, 1.389174e-04;
quarter circle N = 63, 1.679107e-04; half circle N = 41, 1.892882e-04; pentagon N = 107, 2.202269e-04; 15-gon N = 69,
1.240740e-03 (all `MANIFEST.csv`, `relative_gain`).

## 2. `figures/scoreboard.svg`

**Point.** Of 6,071 released records, 4,063 still stand; most of the ones others beat are in the square table.

**Caption.** Every released record, per container, split into those still the best known on 2026-10-18 (blue), those
others have since beaten (orange) and those since tied (green). The number at each bar is "still best known of
released". Most of the 1,974 beaten records are in the square table (1,727 of 2,213); 34 cube records were tied. We no
longer claim any beaten or tied record; each is listed with its source in SUPERSEDED.csv.

**Alt text.** Horizontal stacked bar chart of 21 containers. The square has 2,213 records of which 486 still stand; the
cube 455 with 232 standing and 34 tied; 6-D balls 77 with 23 standing; most other containers are entirely blue.

**Sources.**
| item | file | field |
|---|---|---|
| released per table (6,071) | `MANIFEST.csv` | rows per `shelf` |
| beaten (1,974) / tied (34) | `SUPERSEDED.csv` | `status` = BEATEN or BEATEN_BY_LARGER_N -> beaten; TIED -> tied |
| still best known (4,063) | computed | released - beaten - tied |
| container names | `RESULTS_TABLE.md` | summary table, column `container` |
| square 1,727 beaten / 486 standing; cube 34 tied | `SUPERSEDED.csv`, `MANIFEST.csv` | counts per table (match the README table) |

## 3. `figures/triangle16_two_optima.svg`

**Point.** For 16 circles in a triangle the best radius is now proved, and the proof covers both known arrangements
without choosing between them.

**Caption.** The two known best packings of 16 equal circles in an equilateral triangle (found in 1995 by Melissen &
Schuur and by Graham & Lubachevsky), drawn from this repository's proof certificate. Our computer-assisted proof shows
that no packing has a larger radius: r = 0.078655749492... for a triangle of side 1. The two arrangements share 11
circles (blue) and place 4 differently (orange); arrangement A is shown mirrored so the shared circles line up. The
dashed circle is loose and can move inside its hole; one valid position is shown. Lines mark the 20 circle-to-circle
contacts the proof uses. The proof does not say which arrangement is optimal, or whether both are.

**Alt text.** Two equilateral triangles, each holding 16 circles: 11 blue circles in the same places in both, 4 orange
circles placed differently, and one dashed loose circle near the centre; thin lines join touching circles.

**Sources.**
| item | file | field |
|---|---|---|
| frame points (15 per arrangement) | `proofs/triangle16/out/local_cert16.json` | `data[0].frame_uv` (A), `data[1].frame_uv` (B), barycentric (a, b) |
| loose circle position | same | `data[q].rattlers_uv` |
| 20 contacts | same | `data[q].pairs` |
| m (point form), r16 = m / (2 (1 + sqrt3 m)) | same, `data[0].m_lo_sq`; formula `proofs/triangle16/README.md` Section 1 | r16 = 0.078655749492... (README: 0.0786557494924822861883...) |
| 11 shared / 4 differ; A mirrored by a <-> b | `proofs/triangle16/README.md` Section 4 table, row q = 1; asserted by the script | |
| 1995 discoverers; "does not say which" | `README.md` "In plain words"; `proofs/triangle16/README.md` Section 2 | |

Geometry: points are mapped to the plane with vertices (0, 0), (1, 0), (1/2, sqrt3/2); circles of radius m/2 around the
points fit in the triangle pushed outward by m/2, which is the circle form scaled (README Section 1, "Circle form").
