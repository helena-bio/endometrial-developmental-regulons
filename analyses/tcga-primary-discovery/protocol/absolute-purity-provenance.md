ABSOLUTE PURITY-SOURCE PROVENANCE + CPE SAMPLE-FLOW RE-CONFIRMATION
================================================================================
TASK-028 -- v3 amendment support for the sealed Freeze B packet (experimenter)
STATUS: FACTS + PROVENANCE RECOMMENDATION. Nothing frozen (sealing is Vladimir-gated).

SCOPE / DISCIPLINE
--------------------------------------------------------------------------------
Tumour-blind. Read ONLY: the two purity tables and the frozen cohort
barcodes+subtype labels. Purity VALUE columns were read ONLY to determine
FINITE-value coverage (a covariate-coverage fact required for complete-case
sizing). NO tumour EXPRESSION value was read; NO expression scoring; NO model
fit; the 35 normals were not touched. Source selection is by PROVENANCE /
DEFINITION only -- NOT by coverage and NOT by any downstream result.

Git HEAD unchanged: 83503bad47b60193598b2b9ebe819c22c83e8ac1 (start == end).
ASCII only. All outputs under execution/absolute_provenance/.

INPUTS
--------------------------------------------------------------------------------
Cohort (READ-ONLY metadata): task027 freeze_a_redux/cohort_selected_primary.tsv
  507 primary-tumour patients; subtype totals POLE=49 / MMRd=148 / NSMP=147 /
  p53abn=163. Join key: patient-level TCGA-XX-XXXX.

Purity source A (Aran CPE table), Supplementary Data 1:
  experimenter_final/sources/aran_cpe.xlsx
  sha256 fce1ba13ff3a6432b7e4f260eb516e407c967d615ca87bc8dffa7765baf97529
  Sheet 'Supp Data 1'; HEADER ON ROW 4 (rows 1-2 title, row 3 blank).
  Columns row4: Sample ID | Cancer type | ESTIMATE | ABSOLUTE | LUMP | IHC | CPE.
  Aran, Sirota & Butte 2015, Nat Commun 6:8971, PMID 26634437. This one file
  carries BOTH an ABSOLUTE column AND the CPE (consensus) column. Missing values
  are the literal string 'NaN'. Sample IDs are aliquot-level (e.g. -01A).

Purity source B (dedicated PanCanAtlas ABSOLUTE):
  experimenter_final/sources/absolute_purity.txt
  = TCGA_mastercalls.abs_tables_JSedit.fixed.txt
  sha256 f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd
  GDC PanCanAtlas, https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
  Columns: array | sample | call status | purity | ploidy | Genome doublings |
  Coverage for 80% power | Cancer DNA fraction | Subclonal genome fraction |
  solution. 10,786 data rows (pan-cancer). 'array' is patient-primary
  (TCGA-XX-XXXX-01); 'sample' is aliquot.

Method note (finite value): a value counts as FINITE only if it is a real
number. Excluded as non-finite: None/empty, and the literal tokens
NaN / NA / NULL / '.'. When multiple aliquots map to one cohort patient, the
patient is finite if ANY of its aliquots is finite (union at patient level).

================================================================================
1. CPE RE-CONFIRMATION (for the v3 sample-flow)   -- ALL CONFIRMED (Y)
================================================================================
Aran CPE column, finite value over the 507 cohort, patient-level:
  POLE 49/49 | MMRd 148/148 | NSMP 146/147 | p53abn 163/163 | ALL 506/507

  - Finite CPE coverage = 506/507.                                        [Y]
  - The single missing patient is TCGA-BS-A0TG, subtype NSMP.             [Y]
  - CPE-adjusted complete-case subtype counts (drop the 1 missing NSMP):
      POLE=49 | MMRd=148 | NSMP=146 | p53abn=163
    Only NSMP drops by 1 (147 -> 146); POLE, MMRd, p53abn unchanged.      [Y]

(TCGA-BS-A0TG is present in the Aran table but its CPE cell is 'NaN', so it is
non-finite and correctly excluded from the CPE complete-case set.)

================================================================================
2. ABSOLUTE PROVENANCE -- EVERY ABSOLUTE SOURCE ON DISK
================================================================================
Two ABSOLUTE sources exist on disk. Facts for each (finite-value coverage,
patient-level, over the 507 cohort):

--- ABSOLUTE source (1): Aran-2015 ABSOLUTE COLUMN --------------------------------
  file    experimenter_final/sources/aran_cpe.xlsx
  sha256  fce1ba13ff3a6432b7e4f260eb516e407c967d615ca87bc8dffa7765baf97529
  column  'ABSOLUTE' (Supp Data 1, header row 4)
  release Aran, Sirota & Butte 2015 aggregation of ABSOLUTE (Carter et al. 2012,
          DNA copy-number ML purity) as available at that time; partial subset.
  FINITE coverage of the 507:
    POLE 27/49 | MMRd 100/148 | NSMP 105/147 | p53abn 134/163 | ALL 366/507
  complete-case subtype counts: POLE=27 | MMRd=100 | NSMP=105 | p53abn=134

--- ABSOLUTE source (2): TCGA PanCanAtlas ABSOLUTE (dedicated file) ---------------
  file    experimenter_final/sources/absolute_purity.txt
          (= TCGA_mastercalls.abs_tables_JSedit.fixed.txt)
  sha256  f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd
  column  'purity'
  release TCGA PanCanAtlas ABSOLUTE purity/ploidy master calls; GDC UUID
          4f277128-f793-4354-a13d-30cc7fe9f6b5; uniform pan-cancer ABSOLUTE run.
  FINITE coverage of the 507:
    POLE 48/49 | MMRd 148/148 | NSMP 144/147 | p53abn 162/163 | ALL 502/507
  complete-case subtype counts: POLE=48 | MMRd=148 | NSMP=144 | p53abn=162
  finite-missing (5): TCGA-A5-A0G1 (POLE), TCGA-A5-A1OH (p53abn),
    TCGA-BG-A0M9 (NSMP), TCGA-D1-A0ZU (NSMP), TCGA-EY-A1GT (NSMP)

RECONCILIATION WITH THE PRE-SEAL AUDIT (503/507 by ROW PRESENCE):
  The pre-seal COVARIATE_FINAL_AUDIT.md / task3 purity JSON reported ABSOLUTE
  503/507 by ROW PRESENCE (a barcode is present in the table), NOT by finite
  value. Re-running ROW PRESENCE here reproduces exactly 503/507
    (POLE 48/49 | MMRd 148/148 | NSMP 144/147 | p53abn 163/163 | ALL 503/507).
  The finite-value count is 502/507. The single-patient gap is TCGA-A5-A1OH
  (p53abn): its barcode IS present as a row (array=TCGA-A5-A1OH-01) but its
  'call status' and 'purity' cells are BLANK (uncalled sample), so it is present
  by row but NON-FINITE by value. Complete-case (OPTION 2) requires a finite
  value, hence 502 (not 503). This is a definition difference, not a mismatch:
  ROW PRESENCE 503 == pre-seal audit; FINITE complete-case 502 == the number the
  v3 complete-case ABSOLUTE sensitivity must use.

  Also note the pre-seal audit's 4-missing set (BS-A1OH not among its named,
  it listed POLE 1 + NSMP 3). The named row-presence misses here are
  A5-A0G1 (POLE), BG-A0M9 / D1-A0ZU / EY-A1GT (NSMP) = 1 POLE + 3 NSMP, i.e.
  identical to the pre-seal 4. A1OH is the 5th, added ONLY under the finite
  definition. All internally consistent.

================================================================================
3. WHY THE SOURCES DIFFER (366 vs 502) -- DEFINITIONAL / RELEASE REASON
================================================================================
The spread is a definitional/release difference between an OLDER PARTIAL
aggregation and the UNIFORM pan-cancer re-run:

  - Aran-2015 ABSOLUTE COLUMN is a 2015-vintage table that carried an ABSOLUTE
    value only where an ABSOLUTE call happened to be available at that time.
    Whole cancer types are BLANK: in the Aran table ACC has 0/92 finite ABSOLUTE
    and LIHC has 0/380 finite ABSOLUTE, while their CPE columns are ~fully
    populated. For UCEC specifically the Aran ABSOLUTE column is finite in only
    396/553 table rows (157 UCEC rows are 'NaN'). Over the 507 cohort this
    partial availability yields 366/507. The gaps are missing-DATA in the older
    aggregation, not "impure" tumours.

  - TCGA PanCanAtlas ABSOLUTE (Taylor-2018-era PanCanAtlas program; the master
    calls distributed on GDC) is the UNIFORM, re-run ABSOLUTE across all 33 TCGA
    types with a single called-solution table. It calls the large majority of
    UCEC samples, giving 502/507 finite over the cohort (503 present by row).

Provenance summary: same underlying METHOD name ("ABSOLUTE", Carter et al.
2012), but two DIFFERENT tabulations/releases -- Aran-2015's partial pan-cancer
subset vs the later uniform PanCanAtlas master calls. The ~136-patient gap
(502-366) is entirely explained by cancer-types/samples that were unpopulated in
the older Aran ABSOLUTE column but called in the PanCanAtlas run.

Provenance citations:
  - Aran D, Sirota M, Butte AJ. "Systematic pan-cancer analysis of tumour
    purity." Nat Commun 6:8971 (2015), PMID 26634437. Supp Data 1 tabulates
    ESTIMATE/ABSOLUTE/LUMP/IHC/CPE; ABSOLUTE is populated only where available.
  - TCGA PanCanAtlas, GDC publications page, ABSOLUTE purity/ploidy file
    TCGA_mastercalls.abs_tables_JSedit.fixed.txt (GDC UUID
    4f277128-f793-4354-a13d-30cc7fe9f6b5); uniform pan-cancer ABSOLUTE master
    calls used across the PanCanAtlas program.

================================================================================
4. PROVENANCE-RULE CHOICE (definition-only; NOT coverage-driven)
================================================================================
RULE (stated independent of any coverage or downstream number):
  Honour the sealed B2.6 spec text, which names the ABSOLUTE sensitivity as
  "TCGA PanCanAtlas ABSOLUTE". Bind that name to the DEDICATED PanCanAtlas
  master-calls file, NOT to the Aran-2015 ABSOLUTE column. Adopt it only if it
  is unambiguous and sufficiently documented on disk: one clear release, one
  clear purity column, one pinnable checksum, a documented source URL.

CHECK (is source (2) unambiguous + sufficiently documented?):
  release  YES -- single named GDC master-calls file, uniform pan-cancer run.
  column   YES -- single 'purity' column, unambiguous.
  checksum YES -- sha256 f430a975...763cd, pinnable on disk.
  source   YES -- GDC UUID 4f277128-... on the PanCanAtlas publications page.
  All four criteria met -> the ABSOLUTE sensitivity is PINNABLE. No need to
  fall back to REMOVING the ABSOLUTE sensitivity.

CHOICE: TCGA PanCanAtlas ABSOLUTE (source 2), the dedicated file.
  Rejected: the Aran-2015 ABSOLUTE column. Reason is PROVENANCE, not coverage:
  the sealed spec names "PanCanAtlas ABSOLUTE" and the Aran column is a
  different, older, partial tabulation of ABSOLUTE that does not match that
  name. (Its far lower 366/507 coverage is a consequence of, not the reason for,
  the choice.)

================================================================================
5. CHOSEN ABSOLUTE SOURCE -- EXACT PIN FOR v3 (complete-case sensitivity)
================================================================================
  source   TCGA PanCanAtlas ABSOLUTE purity/ploidy master calls
  file     experimenter_final/sources/absolute_purity.txt
           (= TCGA_mastercalls.abs_tables_JSedit.fixed.txt)
  sha256   f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd
  column   purity
  release  GDC PanCanAtlas, UUID 4f277128-f793-4354-a13d-30cc7fe9f6b5
  join     patient-level TCGA-XX-XXXX (from 'sample' aliquot; 'array' fallback)
  complete-case definition: FINITE purity value (a called ABSOLUTE solution);
           present-but-uncalled rows (e.g. TCGA-A5-A1OH) are NOT complete-case.
  complete-case n = 502/507
  complete-case subtype counts:
     POLE = 48   (drops 1: TCGA-A5-A0G1)
     MMRd = 148  (unchanged)
     NSMP = 144  (drops 3: TCGA-BG-A0M9, TCGA-D1-A0ZU, TCGA-EY-A1GT)
     p53abn = 162 (drops 1: TCGA-A5-A1OH, present-but-uncalled)

CONTRAST WITH THE CPE PRIMARY (for the v3 sample-flow table):
  CPE primary   complete-case 506/507  -> POLE 49 | MMRd 148 | NSMP 146 | p53abn 163
  ABSOLUTE sens complete-case 502/507  -> POLE 48 | MMRd 148 | NSMP 144 | p53abn 162

--------------------------------------------------------------------------------
Machine numbers: absolute_provenance.json (full) and absolute_counts.json
(compact). Manifest: MANIFEST.json. No scoring, no model, no tumour value read.
