# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make          # standard build
make asan=1   # build with AddressSanitizer
make clean    # remove objects and binary
```

Requires: g++/gcc, zlib (`-lz`), pthreads. No C++11 — the codebase intentionally uses C99/C++03 for portability.

The binary is `./hifiasm`.

## Testing

There is no unit test suite. CI (`.github/workflows/ci.yaml`) only verifies compilation. Functional testing is done by running the assembler on real read data and inspecting output GFA files.

## Architecture Overview

Hifiasm is a haplotype-resolved de novo genome assembler for PacBio HiFi and ONT reads. The main pipeline entry is `ha_assemble()` in [Assembly.cpp](Assembly.cpp), called from [main.cpp](main.cpp).

### Assembly Pipeline

1. **Read loading & error correction** — [Process_Read.cpp](Process_Read.cpp) loads FASTA/FASTQ input; [Correct.cpp](Correct.cpp) runs POA-based consensus correction, caching results to `.ec.bin`
2. **Overlap detection** — [Hash_Table.cpp](Hash_Table.cpp) builds k-mer index; [Overlaps.cpp](Overlaps.cpp) detects, validates, and chains read overlaps, caching to `.ovlp.*.bin`
3. **Graph construction & phasing** — [Assembly.cpp](Assembly.cpp) builds the overlap graph, separates haplotypes using heterozygous bubble detection
4. **Output** — [gfa_ut.cpp](gfa_ut.cpp) writes GFA format; primary contigs go to `.bp.p_ctg.gfa`, haplotype contigs to `.bp.hap1/2.p_ctg.gfa`

### Key Source Files

| File | Role |
|------|------|
| [Assembly.cpp](Assembly.cpp) | Core assembly graph logic, haplotig separation |
| [Overlaps.cpp](Overlaps.cpp) | Overlap detection, validation, chaining |
| [Correct.cpp](Correct.cpp) | Error correction pipeline |
| [Hash_Table.cpp](Hash_Table.cpp) | K-mer indexing for all-vs-all overlap |
| [POA.cpp](POA.cpp) | Partial Order Alignment consensus |
| [Trio.cpp](Trio.cpp) | Trio-binning for haplotype resolution |
| [HIC.cpp](HIC.cpp) | Hi-C read integration |
| [ecovlp.cpp](ecovlp.cpp) | Error-corrected overlap handling |
| [Purge_Dups.cpp](Purge_Dups.cpp) | Haplotig duplication purging |
| [CommandLines.cpp](CommandLines.cpp) | All CLI argument parsing |

### Key Data Structures

- **`hifiasm_opt_t`** (`CommandLines.h`) — central options struct; passed throughout the pipeline
- **`All_reads`** — compressed read storage with indexed access
- **`overlap_region`** — represents a pairwise read overlap with CIGAR alignment
- **`UC_Read`** — uncompressed single read

### Threading & Memory

- Threading via pthreads wrapped in [kthread.cpp](kthread.cpp)/[kthread.h](kthread.h)
- Memory: manual malloc/free; [kalloc.h](kalloc.h) provides NUMA-aware pool allocation
- Performance-critical paths use SSE4.2 and POPCNT intrinsics (especially in [Levenshtein_distance.h](Levenshtein_distance.h))

### Input Modes

| Flag | Mode |
|------|------|
| (default) | HiFi-only |
| `--ont` | ONT R10 |
| `--ul` | HiFi + ultra-long ONT |
| `--h1/--h2` | Hi-C phasing |
| `-1/-2/-3` | Trio binning with parental short reads |

### Caching

Intermediate results are cached to avoid recomputation across runs:
- `.ec.bin` — corrected reads
- `.ovlp.source.bin` / `.ovlp.reverse.bin` — computed overlaps

### EC-Round Overlaps vs Final Overlaps (two distinct PAF populations)

There are two separate overlap computation phases that both write to `R_INF.paf` / `R_INF.reverse_paf`, but they serve different purposes:

**EC-round overlaps** (`worker_hap_ec` → `push_ne_ovlp`, [ecovlp.cpp:2771](ecovlp.cpp#L2771)):
- Run during each error-correction round on *partially-corrected* reads
- Store the **trimmed exact-match sub-region** coordinates via `extract_max_exact` from alignment windows (`w_list`)
- `el=1` means the longest exact substring was found; `qns/qe/ts/te` are trimmed to that sub-region
- Purpose: guide the next POA/consensus correction round

**Final overlaps** (`worker_hap_dc_ec_gen_new_idx` / `worker_hap_dc_ec_gen` → `push_ff_ovlp`, [ecovlp.cpp:2827](ecovlp.cpp#L2827)):
- Run once after all EC rounds (`ha_ec_ff` in [Assembly.cpp:1979](Assembly.cpp#L1979)) on *fully-corrected* reads with a fresh k-mer index
- Store **full overlap chain coordinates** (not trimmed); `el = shared_seed` from `exact_ec_check` on real sequences
- `h_ec_lchain_fast_new` reconciles newly-found chains with EC PAF entries: if coordinates agree within 86.67%, quality metadata (`strong`/`no_l_indel`) is inherited from the EC entries
- Purpose: feed the assembly graph — these are the canonical overlaps used for graph construction

The EC PAF entries cannot feed the graph directly because: (1) read sequences change across EC rounds so coordinates go stale, (2) the graph needs full overlap extents not exact sub-region coordinates.

### Final Overlap Performance (realtime/streaming mode)

`ha_ec_ff` is the dominant bottleneck in streaming mode because it always iterates over **all** reads (`kt_for(..., n_a)`), even when only a small batch of new reads was added. Timing data shows it takes ~1.6 hr and is flat across batch sizes (20–100%), while EC rounds scale proportionally.

The planned incremental fix is in [ecovlp.cpp:6525](ecovlp.cpp#L6525) (currently commented out): for `continue_from_prev_state` runs, only process new reads (`kt_for_mod` over `[total_reads0, n_a)`) and dirty old reads (`kt_for_dirty` over `[0, total_reads0)`). Non-dirty old reads keep their existing PAF entries; their reverse-direction entries get updated as a side effect when new/dirty reads find overlaps to them. **`kt_for_mod` and `kt_for_dirty` do not exist yet** — they need to be implemented in [kthread.cpp](kthread.cpp)/[kthread.h](kthread.h).
