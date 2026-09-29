# Notes for Claude Code

## Documentation layout

- `README.md` stays short: what protal is, a quick start, and links. Do not grow it.
- User documentation for running protal (installation, usage, map files, outputs) is on
  https://protal.earlham.ac.uk/main.php?site=documentation. `docs/` covers what the website does
  not: building from source, database files, building and training a database, simulation,
  qcmsa, testing. Update the page in `docs/` that matches a change in behaviour or options.
- When a change makes the website out of date, say so in your reply (the website is not in this
  repository).

## Audits, benchmarks and reviews

Store every audit, code review, benchmark or experiment report you write in `docs/claude/`,
following the conventions in `docs/claude/README.md`: one file per report named
`YYYY-MM-DD-<topic>.md`, with the commit, data and commands it is based on, and a line in the
table of that README. Do not put such reports in the repository root, next to the code, or only in
a chat reply or scratch directory.
