# 🔥 FIREWATCH SA

**An Interactive Spatial Statistical Analysis of Satellite-Detected Fire Patterns in South Africa (2024–2025)**

🌐 Live dashboard: `https://<account>.shinyapps.io/firewatch-sa/`  ← add after deploying
💻 Code: `https://github.com/<username>/firewatch-sa`

An R Shiny dashboard that maps NASA FIRMS active-fire detections and tests whether their spatial
pattern is statistically meaningful, using three point-pattern methods (Chapter 4: Point Pattern Data):

| Method | Question it answers | Dashboard output |
|---|---|---|
| Kernel Density Estimation (KDE) | Where are detections concentrated? | Density layer on the map + concentration statistics |
| Nearest neighbour (Clark–Evans) | Clustered, random or dispersed? | R ratio, z-test, plain-language verdict, distance histogram |
| Ripley's K-function | At what distances does clustering occur? | L(r) − r / K(r) curve with CSR envelope + global DCLF test |

Plus a **Compare seasons** tab and **Download** buttons for results and filtered detections.

> Each point is a satellite-detected active-fire pixel, not a verified wildfire or burned area.

---

## Project structure

```
firewatch-sa/
├── app.R                      # the dashboard (UI + server)
├── R/                         # auto-loaded by Shiny
│   ├── constants.R            # projection, fire colours, helpers, flame logo
│   ├── analysis.R             # KDE, nearest neighbour, Ripley's K, seasons
│   └── about.R                # text for the About tab
├── scripts/
│   ├── 00_install_packages.R  # run once per computer
│   └── 01_prepare_data.R      # FIRMS file + GADM -> data/*.rds
├── tests/validate_methods.R   # checks the methods on patterns with known answers
├── www/styles.css             # fire styling
├── data/                      # cleaned .rds files (COMMITTED)
└── data-raw/                  # raw downloads (NOT committed)
```

---

## Step 0 — Choose the dataset

All options are from **NASA FIRMS**. Use **one sensor** for the whole project so 2024 and 2025 are comparable.

| Option | File name starts with | Resolution | Pros | Cons | Verdict |
|---|---|---|---|---|---|
| **VIIRS S-NPP, standard** | `fire_archive_SV-C2_` | 375 m | Science-quality, long record (2012+), has the `type` column used to remove industrial hot spots | Larger files | ✅ **Use this** |
| VIIRS NOAA-20, standard | `fire_archive_J1V-C2_` | 375 m | Same quality | Record only from 2018 | Good backup |
| VIIRS NOAA-21 | `fire_archive_J2V-C2_` | 375 m | Newest | Very short record | Avoid |
| MODIS C6.1 | `fire_archive_M-C61_` | 1 km | Record back to 2000, smaller files | Misses small fires | Only if VIIRS is too slow |
| Near-real-time | `fire_nrt_…` | — | Very recent | No `type` column; gets replaced by standard data | Avoid |

**Your file `fire_archive_SV-C2_828855` is already the recommended dataset** (SV = Suomi-NPP VIIRS,
C2 = Collection 2, *archive* = standard science-quality data). Just check it covers
**1 Jan 2024 – 31 Dec 2025** and **South Africa**. If the download zip also contains a `fire_nrt_…`
file, don't use it; the script ignores it.

### "It's an Excel spreadsheet"

FIRMS delivers **CSV** files. Windows shows CSVs with the Excel icon, so check the real extension
(File Explorer → *View* → tick *File name extensions*).

* If it is `fire_archive_SV-C2_828855.csv` → perfect, use it as it is.
* **Do not open and re-save it in Excel.** Excel stops at 1 048 576 rows (two years of South African
  VIIRS data can approach that) and it rewrites the dates.
* If you only have an `.xlsx` version, the script reads that too (and fixes Excel dates), but
  re-download the CSV if the row count is exactly 1 048 575.

Copy the file into `data-raw/`. Province boundaries are downloaded automatically (GADM).

---

## Step 1 — Install tools (everyone, ~20 min)

1. Install **R ≥ 4.3** (cran.r-project.org), **RStudio Desktop** (posit.co) and **Git** (git-scm.com).
2. Create free accounts on **github.com** and **shinyapps.io**.
3. RStudio → *Tools → Global Options → Git/SVN* → check that the Git executable is found (restart RStudio after installing Git).
4. Tell Git who you are (RStudio → *Terminal* tab):
   ```bash
   git config --global user.name  "Your Name"
   git config --global user.email "you@students.wits.ac.za"
   ```
5. GitHub needs a token instead of a password. In the R console:
   ```r
   install.packages(c("usethis", "gitcreds"))
   usethis::create_github_token()   # opens GitHub; tick "repo", generate, copy the token
   gitcreds::gitcreds_set()          # paste the token
   ```

---

## Step 2 — Create the GitHub repository (member 4)

1. Unzip this project folder somewhere simple, e.g. `Documents/firewatch-sa`.
2. Double-click `firewatch-sa.Rproj` to open it in RStudio.
3. On GitHub: **New repository** → name `firewatch-sa` → **Public** → do **not** add a README → *Create*.
4. In the RStudio *Terminal*:
   ```bash
   git init
   git add .
   git commit -m "Initial FIREWATCH SA scaffold"
   git branch -M main
   git remote add origin https://github.com/<username>/firewatch-sa.git
   git push -u origin main
   ```
5. GitHub → repo → *Settings → Collaborators* → invite the other three members.
6. **Everyone else:** accept the invite, then RStudio → *File → New Project → Version Control → Git* →
   paste `https://github.com/<username>/firewatch-sa.git` → *Create Project*.

### Working rules (avoid merge conflicts)

Each person works on their **own branch** and only edits **their own files**:

| Member | Owns | Branch | Files |
|---|---|---|---|
| 1 | Data download, cleaning, boundaries | `data-prep` | `scripts/01_prepare_data.R`, `data/` |
| 2 | Methods + validation | `methods` | `R/analysis.R`, `tests/validate_methods.R` |
| 3 | Interface, map, styling | `ui` | `app.R` (UI part), `www/styles.css`, `R/about.R` |
| 4 | Integration, testing, GitHub, deployment | `main` | `app.R` (server part), `README.md` |

Daily routine:

```bash
git checkout main && git pull           # 1. get everyone's latest work
git checkout -b methods                 # 2. first time only (use your branch name)
git checkout methods && git merge main  #    later times: switch to your branch, bring it up to date
# ... edit, save, test with shiny::runApp() ...
git add .
git commit -m "Add Ripley's K envelope"  # 3. small, descriptive commits
git push -u origin methods              # 4. push, then open a Pull Request on GitHub
```

Member 4 reviews and merges Pull Requests into `main`. (RStudio's *Git* tab does the same with buttons:
Pull ⬇, Commit, Push ⬆.)

---

## Step 3 — Install packages (everyone)

```r
source("scripts/00_install_packages.R")
```

---

## Step 4 — Prepare the data (member 1)

1. Put `fire_archive_SV-C2_828855.csv` (or `.xlsx`) in `data-raw/`.
2. Open `scripts/01_prepare_data.R`. Leave `TEST_MODE <- TRUE` and run:
   ```r
   source("scripts/01_prepare_data.R")
   ```
   This processes **August 2024** only and prints the Limpopo nearest-neighbour ratio. Nothing is saved.
3. If that works, set `TEST_MODE <- FALSE` and run it again (a few minutes for the full file).
4. Sanity-check the printed tables: Limpopo, Mpumalanga, KwaZulu-Natal and the Eastern Cape should
   have the most detections, peaking in winter/early spring (Jul–Sep).
5. Commit the four files in `data/` and push:
   ```bash
   git add data && git commit -m "Add cleaned FIRMS data" && git push
   ```

**What the script does:** keeps 2024–2025 → removes low-confidence detections → keeps only
`type = 0` (vegetation fires; removes e.g. industrial heat) → removes duplicates → projects to
South African Albers (km) → assigns each point to a province (GADM) → adds month and season →
builds `spatstat` windows for each province and the whole country (Prince Edward Islands removed)
→ saves a cleaning log that the dashboard shows on the About tab.

If GADM is down: download *South Africa → GeoPackage* from gadm.org, put `gadm41_ZAF.gpkg` in
`data-raw/`, and re-run.

---

## Step 5 — Validate the methods (member 2)

```r
source("tests/validate_methods.R")
```

Simulates a **random**, a **clustered** (Thomas process) and a **dispersed** (inhibition) pattern and
checks that both nearest neighbour and Ripley's K classify each one correctly. Screenshot the table and
the three plots for the presentation — it shows the methods work *before* they are used on real data.

---

## Step 6 — Run the dashboard (everyone)

```r
shiny::runApp()
```

Try:
1. *Limpopo · 2024 · August · Ripley's K* → **Run analysis**
2. *South Africa · Both years · All months · KDE* → **Run analysis** (density layer appears on the map)
3. *Compare seasons* tab → *South Africa* → **Compare seasons**

### How the app is built

* `app.R` loads the `.rds` files once. Filters update the **detection count** and map dots instantly.
* **Run analysis** builds a `spatstat` point pattern for the chosen region/period and always computes
  the nearest-neighbour ratio (for the value boxes), plus the chosen method.
* Ripley's K uses at most 2 000 randomly sampled points (adjustable). Random thinning does not change
  K(r), so the curve is valid and the app stays fast.
* The map shows at most 8 000 dots for speed; the statistics always use **all** detections.

---

## Step 7 — Deploy to shinyapps.io (member 4)

```r
install.packages("rsconnect")
# shinyapps.io -> Account -> Tokens -> Show -> Copy to clipboard, then paste & run, e.g.
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")

rsconnect::deployApp(
  appDir   = ".",
  appName  = "firewatch-sa",
  appFiles = c("app.R",
               list.files("R",    full.names = TRUE),
               list.files("www",  full.names = TRUE),
               list.files("data", full.names = TRUE))
)
```

* Your link: `https://<account>.shinyapps.io/firewatch-sa/` — put it at the top of this README.
* **Never commit the token** (`rsconnect/` is in `.gitignore`).
* If the app crashes online, check *shinyapps.io → Applications → firewatch-sa → Logs*. The free plan
  has 1 GB of memory; if it runs out, set *Instance size* larger or reduce the data to type = 0, nominal
  + high confidence (already done) and re-deploy.
* Re-deploy after every merge into `main` that changes the app.

---

## Step 8 — Timeline to 26 October

| Day | Goal |
|---|---|
| Thu 22 | Everyone installed; repo created and cloned; data prepared and committed |
| Fri 23 | Methods validated; app running locally for everyone; first deploy |
| Sat 24 | Polish UI, interpretation text, check numbers; test on a phone |
| Sun 25 | Final deploy; rehearse the demo; slides |
| Mon 26 | 14:00 presentation — share GitHub + shinyapps.io links |

### Final checklist
- [ ] App runs from a **fresh clone** (`git clone …` → `shiny::runApp()`)
- [ ] Deployed link works on a phone and another laptop
- [ ] README has both links
- [ ] Demo rehearsed (~2 min per method)
- [ ] Limitations slide: detections ≠ fires; CSR assumes constant intensity; pointwise envelopes; huge n makes p-values tiny

### Suggested demo script
1. **Problem** (30 s): fires are visible on the map, but is the pattern statistically meaningful?
2. **KDE** — South Africa, both years: hotspots in the north-east; change σ to show smoothing.
3. **Nearest neighbour** — Limpopo, August 2024: R well below 1 → clustered; read the explanation.
4. **Ripley's K** — same selection: clustering occurs from short distances up to tens of km.
5. **Compare seasons** — winter/spring dominate the north-east; the Western Cape behaves differently.
6. **Validation + limitations** (1 min).

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `data/fires.rds not found` | Run `scripts/01_prepare_data.R` with `TEST_MODE <- FALSE` |
| `More than 1% of dates could not be read` | The file was re-saved by Excel in an unusual date format — use the original CSV |
| `Missing columns: …` | Wrong file (e.g. a shapefile export) — download as CSV |
| `stopifnot(setequal(...))` fails | GADM province names changed; print `prov$province` and fix the spellings in the script |
| Ripley's K is slow | Use 19 simulations or fewer points; smaller max distance |
| Push rejected | `git pull` first, fix any conflict, then push |
| Git asks for a password | Use the token from Step 1 (`gitcreds::gitcreds_set()`) |

---

## Data acknowledgement

We acknowledge the use of data and imagery from NASA's Fire Information for Resource Management
System (FIRMS), part of NASA's Earth Science Data and Information System (ESDIS). Boundaries: GADM v4.1.
