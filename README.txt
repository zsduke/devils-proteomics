Devil's Proteomics
Version 1.0.0-public-rc1
Made by Zhenyao Sun
Department of Chemistry
Duke University

Local table analysis and a separate PDB coordinate viewer.

Included functions
- XLSX, XLS, CSV and TSV input; worksheet and header selection; column mapping.
- Supplied statistics or replicate-based differential-abundance analysis.
- Welch, Student, paired t and Mann-Whitney tests; Benjamini-Hochberg FDR.
- QC, missingness, sample correlation, coefficient of variation and PCA.
- Local GMT gene-set enrichment using the tested identifier population.
- Interactive volcano plots, reports, CSV results and MATLAB volcano export.
- Local PDB model/residue selection, user annotations and observed-sequence export.

Direct mass-spectrometry RAW input is not supported. Convert or process instrument
files with suitable software and supply a quantified table. The structure viewer
accepts a defined subset of PDB; it does not reconstruct missing coordinates.
See USER_GUIDE.txt for input limits and procedures.

Start
Use Python 3.12 with the packages listed in requirements-public.txt. No Python
environment, wheel or native library is bundled. Installation from a package
index requires network access. The listed versions record an existing runtime;
a fresh installation on another computer has not been verified.

  python -m pip install -r requirements-public.txt
  python public_app.py --check
  python public_app.py --port 8765 --data-dir ./analysis-data

Run these commands from the extracted source directory in a dedicated Python
environment. Open http://127.0.0.1:8765/ manually. Choose a different port if it
is occupied. The program does not stop another service or open a browser.
Use a new data directory. Press Ctrl+C in the terminal to stop the service.

For only the coordinate viewer, open static/structure.html in a compatible
browser. It uses bundled assets and does not require Python. Three-dimensional
display requires WebGL. Browser acceptance on a fresh installation is unverified.

Data handling
Table uploads and analysis outputs remain in the chosen local data directory.
PDB and annotation files are read by the browser and are not sent to the table
analysis API. No database connection is included. Download complete result
bundles before closing the session; the result list is not restored on restart.
External documentation links, if opened, are handled by the browser. Local
processing does not guarantee the security of the machine, browser or storage.

Verification scope
The recorded 36 offline checks comprise 19 PDB parser/selection/export checks,
9 mocked-DOM viewer checks, 3 bundled 3Dmol data-parser/GLModel checks, and 5
simulated HTTP/static-route boundary checks. These are not 36 statistical tests.
They do not establish real-browser rendering, WebGL behavior, live HTTP behavior,
MATLAB execution or fresh-install compatibility.

License
This is source-available software under the Devil's Proteomics Education and
Public Sharing License, LicenseRef-Devils-Proteomics-Education-Public-Sharing-1.0,
effective Nov 2, 2026. It is not an OSI open-source license. LICENSE.txt and the
accompanying RELEASE_MANIFEST.json govern the covered materials and release.

Permission is limited to eligible institutional users and noncommercial teaching,
learning or academic research. Covered improvements and validation materials must
be shared publicly by the earlier of 90 calendar days after completion or their
first external distribution, publication or service use. This summary does not
replace the full conditions, including legally mandatory disclosure limits.
Read the license before undertaking work that cannot meet these conditions.

Excluded private prediction technology is not licensed by this release. Existing
MIT rights in the seven identified source modules and third-party licenses remain
in force. See licenses/OWNER_MIT_SCOPE.txt and THIRD_PARTY_NOTICES.txt.
