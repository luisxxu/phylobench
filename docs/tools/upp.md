# UPP (`align`, aln `upp`)

UPP (SEPP 4.5.6 `run_upp.py`) aligns the queries to the fixed backbone alignment and tree, as an alternative to `hmmalign` for the `bscampp` arms. It is an alignment lever (T29), not production parity: `aln=upp` cells exist only for `bscampp` (and the gated `bscampp-b5000`) on `gg2-wol` at B10k\*, are labelled by `aln`, and are never pooled with `hmmalign` cells (user-approved exception to the "same hmmalign alignment" invariant). Rules: `workflow/rules/upp.smk` (`upp_align`, `upp_queries`); no arm file.

## Command line

`{cell}` is `results/{dataset}/{backbone}/{region}/align`. `-d` and `-p` need absolute paths; the rule prefixes `$PWD`.

```sh
work=$PWD/{cell}/upp; rm -rf $work; mkdir -p $work/out $work/tmp
phylobench-measure --out {cell}/upp.run.json --threads {threads} --timeout 84000 \
    --tool-version 'sepp 4.5.6 (run_upp.py)' --env-lock envs/sepp.linux-64.pin.txt -- \
    run_upp.py -s results/{dataset}/regions/{region}/queries.fasta \
    -a data/{dataset}/{backbone}/aligned-dna-sequences.fasta \
    -t data/{dataset}/{backbone}/tree.nwk -A 100 -m dna -x {threads} \
    -o upp -d $work/out -p $work/tmp
```

Output goes to `{cell}/upp.log`. On `ok` the rule moves `upp_alignment_masked.fasta` to `{cell}/upp.masked.fasta` and `upp_insertion_columns.txt` to `{cell}/upp.insertion_columns.txt`, then removes `$work/tmp` (`$work/out` keeps the unmasked alignment). `ok` without the masked file hard-fails. On `status != ok` both outputs are empty files so the DNF row is recorded and the DAG continues; `upp_queries` then writes an empty `upp.aligned.fasta`.

`upp_queries` (unmeasured, phylobench env) runs `phylobench.regions.upp_query_rows` on the masked file: every row has the backbone width, every backbone row comes back identical after `upper()`, no unknown names, no lowercase or `.` in a query row. Any violation fails the job. Queries UPP dropped are omitted and scored as unplaced misses.

## Version and `-A`

`envs/sepp.yaml`, lock `envs/sepp.linux-64.pin.txt` (the `sepp-pinned` env). `pplacer --version` prints `dev` in this build; alpha23 is proven by the lock.

`-A 100` (decomposition subset size) is the single pinned setting. S0 on Barnacle: `-A 10` at B10k\* with 300 full-length queries timed out at 2,400 s (about 1,900 HMMs x 16 fragment chunks, about 30k `hmmsearch` launches); `-A 100` finished in 381.8 s.

## Threads and timeout

`-x {threads}` (16 in production). The per-rule harness `--timeout 84000` comes after the global `--timeout` in `{measure}` (the later flag wins), and equals the profiles' 1,440 min `upp_align` runtime minus 40 min so the harness writes a DNF row instead of SLURM killing the job silently. `tests/test_profile.py` ties the two.

## Hidden-binary check (S0 U1/U2)

- `run_upp.py -h` lists `-s -a -t -A -m -x -o -d -p`.
- `main.config` paths all sit under the env's `bin`; `~/.sepp` is absent; `sepp/config.py` mentions `~/.sepp` only in a docstring.
- The options dump lists `pasta path='run_pasta.py'` (unpinned), but PASTA runs only without a backbone. We always pass `-a` and `-t`, so it never runs.

## S0 numbers

| Probe | Result |
|---|---|
| Tiny fixture, `-A 10`, 2 threads | ok, 13.5 s, 3.2 GB; checker all True (width 1500, 1 insertion column, 10/10 queries) |
| B10k\*, 300 queries, `-A 10`, 16 threads | timeout at 2,400 s |
| B10k\*, 300 queries, `-A 100`, 16 threads, EPYC 7302 | ok, 381.8 s, 4,617 MB (cgroup); backbone width 1598, 10,000 rows; full width 1660; masked width 1598; 62 insertion columns; full minus insertion columns equals masked; backbone rows identical after `upper()`; 300/300 queries; 0 lowercase in masked queries; masked query characters `- A C G T` |
| Bfull\*, 300 queries, `-A 100`, `highmem` Xeon Gold 6242 (S0b) | timeout at 8,402 s, 31,102 MB: B2 unsized, no upp cell in T30 (`docs/runbook-t29.md`) |

The backbone input has IUPAC codes (`A B C G K M R S T W Y -`); UPP uppercases but does not change them.

## T29 S1 (B10k\*, EPYC 7302)

UPP align at full-length: 3,186 s and 4,694 MB for 3,000 queries (`hmmalign` 62 s); amplicon regions 312 to 835 s. `upp/bscampp` full-length e2e 3,248 s against 115 s for `hmmalign/bscampp`, for a within-3 gain of +0.0047 (lower bound -0.0007, not shown). Details: `docs/runbook-t29.md`.

## Smoke

UPP runs in the CI smoke for the `bscampp` leg (tiny fixture, 13.5 s) and passed the exact-placement check (`smoke_min_exact: 8`) in the T29 PR.

## Quirks

- UPP uppercases the backbone on read; all-gap backbone columns are kept, so the masked file has exactly the backbone width.
- Cost: UPP builds its HMM ensemble inside each region's `align` row, while `hmmalign` reuses one per-backbone `hmmbuild` that `end_to_end` excludes. Compare `upp` e2e against `hmmalign` e2e plus the `hmmbuild` row (`docs/runbook-t29.md`), and only within one `cpu_model`.
- A missing query (UPP dropped it) is a miss, counted and reported.
- Scratch is cell-local (`-p`): a tmpfs `TMPDIR` would count toward cgroup `memory.peak`.

## License

SEPP/UPP is GPL-3; the engine is HMMER 3.4 (BSD-3). `upp/bscampp` cells also run EPA-ng (AGPL-3). See `LICENSES.md`.
