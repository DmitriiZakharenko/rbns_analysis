# Sheet 2 – Merge sort improvements: brief answers

## (a) Insertion sort on short ranges
**Done correctly.** In `merge_sort_aux`, when `use_insertion_sort=True` and the range size is fewer than 8 elements (`(hi - lo + 1) < 8`), the code uses the provided `insertion_sort`. With `use_role_swap`, it copies from `array_from` to `array_to` and then sorts `array_to` in place so the result ends up in the right array.

---

## (b) Skip merge when already sorted
**Done correctly.** Before calling `merge`, the code checks `use_merge_improvement and array_from[mid] <= array_from[mid + 1]`. If true, the segment is already sorted, so merge is skipped. With `use_role_swap`, it copies `array_from[lo..hi]` to `array_to[lo..hi]` so the final result is still in `array_to`.

---

## (c) Role swap (no copy-back)
**Done correctly.**  
- **Outer call:** When `use_role_swap`, `merge_sort` calls `merge_sort_aux(tmp_array, array, ...)` so the “from” array is `tmp`, “to” is the original `array`.  
- **Recursion:** In `merge_sort_aux`, when `use_role_swap` the recursive calls use `(array_to, array_from)` so the roles alternate with depth.  
- **Merge:** When `use_role_swap`, `merge` does not copy from `array_to` back to `array_from`.  
- **(a) and (b):** With role swap, both the insertion-sort path and the “already sorted” path write the result into `array_to` (by copying then sorting, or by copying the segment), so the sorted data stays in the correct array at each level.

**Expert question – Why does this work?**  
At each recursion level we alternate which array is “from” and which is “to”. So after a merge step we write into the current “to” array. Deeper levels then read from that same buffer as their “from” and write into the other as “to”. Thus we never need to copy the merged result back into the “from” array: the next level reads from where we just wrote. At the top level we started with `tmp` as “from” and the original `array` as “to”, so the final sorted output is in the original `array`.

---

## (d) Experiments and plots

**What you did:** You ran `run_experiment` for three instance types: `random_unique`, `few_different_values`, and `almost_sorted`. Each run compares five configurations: no improvement, only merge improvement, only insertion sort, only role swap, and all three together.

**How to interpret (and what your numbers show):**

1. **Random unique**
   - **Without improvement** is slowest (e.g. ~0.009–0.019 s in your range).
   - **Role swap** and **insertion sort** each give a clear speedup (fewer copies / better small-block behavior).
   - **Merge improvement** helps a little (skips merge only when the two halves are already in order, which is rare in random data).
   - **All improvements** is fastest (~0.006–0.013 s): the gains combine.

2. **Few different values**
   - **Merge improvement** helps more than on random unique, because equal or near-equal values often make `array_from[mid] <= array_from[mid+1]`, so merges are skipped more often.
   - **Role swap** and **insertion sort** still help; **all improvements** stays best.

3. **Almost sorted**
   - **Merge improvement** helps the most here: many adjacent segments are already in order, so many merges are skipped.
   - **Insertion sort** is good on nearly sorted data.
   - Again **all improvements** is the best.

**Plots:** For each experiment, the plot has “input size” on the x-axis and “seconds” on the y-axis. The curve for “without improvement” is on top; “all improvements” is at the bottom. The other three curves sit between them. Your implementation is correct and the relative positions of the curves match the explanations above.

**Screenshots:** The assignment asks for a screenshot of each plot. Use the three figures produced by the three `run_experiment(...)` cells (random_unique, few_different_values, almost_sorted) and paste them into your submission as requested.

---

## Summary: is everything correct?

- **(a), (b), (c):** Implemented correctly; threshold 8 for insertion sort, merge skip when `array_from[mid] <= array_from[mid+1]`, and role swap (outer call, recursive swap, no copy-back, and (a)/(b) writing to `array_to`) are all correct.  
- **Test:** The example with all improvements gives `[2, 2, 2, 3, 4, 4, 5, 5, 5, 6, 8, 12]` — correct.  
- **Experiments:** The trends (no improvement slowest, all improvements fastest, merge improvement helping more on “few values” and “almost sorted”) are as expected. So yes, you did it right.
