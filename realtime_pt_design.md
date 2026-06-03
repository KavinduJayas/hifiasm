---
name: realtime-minimizer-table
description: "Incremental minimizer table (ha_pt_t) for realtime/continue_from_prev_state mode — design decisions, data structures, and all code changes"
metadata: 
  node_type: memory
  type: project
  originSessionId: 98a88e0d-28a8-4a7f-91a7-ca8693c273d6
---

## Goal

`ha_ec_ff` (final overlap step) rebuilds `ha_idx` from ALL reads on every batch, making it O(total_reads) even when a small batch of new reads is added. The fix: persist `ha_idx` to disk, load it on subsequent batches, and only reprocess new reads + dirty reads (reads re-corrected by EC).

## Design Decisions

### Per-read round array (not per-position)
`ha_idxpos_t` is fully packed (64 bits: rid:28 + pos:27 + rev:1 + span:8). No room for a round field without breaking ABI. Instead, `ha_pt_s` gets `uint8_t *read_mz_round[n_reads]`. Since all positions for a given read are always inserted in the same round, per-read tracking is semantically equivalent to per-position tracking.

### Lazy deletion (no position array shrink)
The position arrays in `ha_pt1_t` are densely pre-allocated based on k-mer counts. There is no efficient in-place remove. Dirty reads' old positions remain in the primary table; they are masked at lookup time via `read_mz_round[rid] == 0xFF` (stale sentinel). Compaction (merge delta into primary, rebuild arrays) is deferred.

### Two-table lookup (primary + delta)
Dirty reads' CORRECTED positions cannot be appended to the primary's position arrays (no free slots). New reads' k-mers may not exist in the primary at all. Both cases are handled by a separate small `ha_idx_delta` table built from new + dirty reads only. At lookup, callers iterate both tables; stale entries from primary are skipped.

### `HAF_INCREMENTAL` pipeline flag
Rather than a new pipeline, the existing `sf_worker_for_mz` template worker gets an early-exit guard: if the flag is set, skip reads where `rid < total_reads0 && !(dirty_reads[rid] & 0x3F)`. This makes both the count and insertion phases process only new + dirty reads, producing the delta table at a fraction of the full cost.

### Opaque struct — accessor instead of header exposure
`ha_pt_s` is defined in `htab.cpp` (uses `yak_pt_t`, a static local type from `KHASHL_MAP_INIT`). Cannot move struct to header. Added `ha_pt_mz_round(const ha_pt_t*)` accessor; callers pre-fetch the raw pointer once per function, use it inline in the hot loop.

### Lightweight disk I/O (`ha_pt_table_save` / `ha_pt_table_load`)
`write_pt_index` / `load_pt_index` serialize All_reads + PAF data too — too heavy for `ha_ec_ff`. New functions write only the position table to `{prefix}.ff_pt`.

---

## Files Changed

### `htab.h`
- `extern ha_pt_t *ha_idx_delta;` — new global for delta table
- `ha_pt_gen(...)` signature: added `int extra_flags` parameter (callers pass 0 for normal builds)
- New declarations: `ha_pt_mark_stale`, `ha_pt_gen_delta`, `ha_pt_table_save`, `ha_pt_table_load`, `ha_pt_mz_round`

### `htab.cpp`
- **`ha_pt_s` struct**: added `uint8_t *read_mz_round`, `uint64_t n_reads`, `uint8_t mz_round`
- **`ha_pt_t *ha_idx_delta`**: global definition (alongside `ha_idx`)
- **Inner `ha_pt_gen(ct, n_thread, is_l)`**: initialise new fields to NULL/0
- **Outer `ha_pt_gen(asm_opt, ...)`**: added `extra_flags` param; ORs it into `extra_flag1/2`; allocates `read_mz_round` (`CALLOC`) after inner call, for both normal and fast (`ha_pt_gen_dp`) modes
- **`sf##_ha_pt_insert_list` (template)** and **`ha_pt_insert_list`**: after writing position fields, set `read_mz_round[rid] = mz_round` if array is present
- **`ha_pt_destroy`**: frees `read_mz_round`
- **`sf##_worker_for_mz` (template)**: HAF_INCREMENTAL guard at top — sets `mz[i].n=0` (skips minimizer computation + insertion) for non-dirty old reads
- **`#define HAF_INCREMENTAL 0x200`**: new pipeline flag
- **`write_pt_index`**: writes `n_reads`, `mz_round`, `read_mz_round` array after position arrays
- **`load_pt_index`**: reads those fields back; allocates and populates `read_mz_round`
- **`ha_pt_mz_round(pt)`**: returns `pt->read_mz_round` (or NULL); used by anchor.cpp to avoid struct access through incomplete type
- **`ha_pt_mark_stale(pt, rs)`**: iterates `rs->dirty_reads[0..total_reads0)`, sets `read_mz_round[rid] = 0xFF` for dirty reads
- **`ha_pt_gen_delta(asm_opt, flt_tab, rs, ...)`**: thin wrapper — calls `ha_pt_gen(..., HAF_INCREMENTAL)`, building a table covering only new + dirty reads
- **`ha_pt_table_save(pt, prefix)`**: writes position table only to `{prefix}.ff_pt` (k, pre, tot, tot_pos, bins, n_reads, mz_round, read_mz_round)
- **`ha_pt_table_load(prefix)`**: reads `{prefix}.ff_pt`, reconstructs `ha_pt_t` with all fields

### `anchor.cpp`
- **`seed1_t`**: added `const ha_idxpos_t *ad; int nd;` for delta table positions
- **All 6 `ha_pt_get` call sites** (`ha_get_new_candidates`, `ha_get_ug_candidates`, `minimizers_qgen0`, `minimizers_qgen0_amz`, `h_ec_lchain_re_gen`, `h_ec_lchain_re_gen3`, `h_ec_lchain_re_gen_srt`):
  - Pre-fetch: `const uint8_t *_ha_mz_rnd = ha_pt_mz_round(ha_idx);` at function top
  - Lookup loop: add `ha_idx_delta` branch → `seed[i].ad`, `seed[i].nd`; accumulate `n_a += n + nd`
  - Primary inner loop: add `if (_ha_mz_rnd && _ha_mz_rnd[y->rid] == 0xFF) continue;`
  - Delta inner loop: identical body over `s->ad[0..nd)` (always valid, no stale check)
  - After both loops: `ab->n_a = k;` to fix up count for downstream sort (stale skips reduce actual anchors)
- Long-read paths (`ha_ptl_get`, `seedl_t`) intentionally NOT modified — delta is built as short-read table only

### `Assembly.cpp`
- All 9 outer `ha_pt_gen` call sites: added trailing `0` argument for `extra_flags`
- **`ha_ec_ff`**:
  - Destroys `ha_idx_delta` at entry and exit
  - On `continue_from_prev_state`: tries `ha_pt_table_load`; if successful → `ha_pt_mark_stale` + `ha_pt_gen_delta` (incremental path); if file not found → falls through to full build
  - On full build: saves result with `ha_pt_table_save` for next batch

---

## Known Limitations / Future Work

1. **Stale entries accumulate**: primary table grows over batches (lazy deletion only). Fix: periodic merge of primary + delta (compact position arrays), or rebuild every N batches.
2. **2-bit `dirty_reads` round field** (bits 6–7) limits last-EC-round tracking to 4 values (0–3). Sufficient for default `number_of_round ≤ 3`.
3. **Delta misses singleton k-mers**: delta count phase applies the same `ha_ct_shrink` cutoff as full builds. New reads' k-mers appearing only once in the delta are dropped. Overlaps between a new unique read and old reads may be missed if the shared k-mer appears ≥2× in the full dataset but only once in the delta's subset.
4. **No delta for long-read paths** (`ha_ptl_get`): UL-mode overlap detection does not benefit from the delta table yet.
