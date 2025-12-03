# Applying the external-data changes locally

Follow these steps to run `hpg_apr.R` with the new external dataset loader:

1. **Update the working directory**
   - Open `hpg_apr.R` and set `setwd("<path-to-your-local-clone>")` so the script can find your data file and write outputs.

2. **Enable external data loading**
   - In `hpg_apr.R`, set `use_external_data <- TRUE`.
   - Point `external_data_path` to your file name if it differs from the default `SVR_input_data.RData`.

3. **Prepare your `.RData` file**
   - Ensure it contains an object named `sim` (time × units matrix with treated units first, then controls).
   - Optionally include `bands`, `t0`, and `num_controls` in the same file to override the defaults. Otherwise the script will derive `num_controls <- ncol(sim) - bands` and keep the existing `bands`/`t0` values.

4. **Run the script**
   - Execute from the repo root: `Rscript hpg_apr.R <sp_range> <t0> <errors_sp> <index>`.
   - When `use_external_data` is `TRUE`, only the `t0` argument affects downstream code; the other args remain for output paths/seeds.

5. **Outputs**
   - Results are saved under `Output/apr_sims/Results/ss<sp_range>/tt<t0>/ee<errors_sp>/res_<index>.RData`.

These steps apply the recent change that lets the master script load your saved `sim` matrix instead of generating simulated data.

## Bringing these changes into your own clone

If you want to pull this Codex-generated commit into your local fork, you can cherry-pick it from this branch:

1. Add a temporary remote that points to the Codex branch (replace `<codex_repo_url>` with the URL of this repository or your downloaded copy):

   ```bash
   git remote add codex <codex_repo_url>
   ```

2. Fetch the branch that contains the fix (named `work` here):

   ```bash
   git fetch codex work
   ```

3. Apply just the commit that introduced the external-data loader docs and wiring:

   ```bash
   git cherry-pick 7110576
   ```

4. Remove the temporary remote if you no longer need it:

   ```bash
   git remote remove codex
   ```

Alternatively, copy `hpg_apr.R` and `LOCAL_SETUP.md` from this repo into the same paths in your fork and commit them normally. Either approach will reproduce the Codex changes in your own repository.
