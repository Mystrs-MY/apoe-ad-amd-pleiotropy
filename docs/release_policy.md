# Public release policy

## Before submission

- Keep the GitHub repository private.
- Freeze a commit that matches the submitted manuscript and supplementary files.
- Run `python workflow/validate_release.py`.
- Confirm that no controlled data, credentials, manuscript drafts, or provider archives are tracked.

## At submission

- Create a versioned release, for example `v1.0.0-submission`.
- Record the commit SHA in the manuscript's Code Availability statement.
- If anonymous review is required, use the journal-approved private or anonymized code-access mechanism.
- Do not make the repository public until the target journal's policy and the authors' disclosure decision permit it.

## At public release

- Make the repository public.
- Archive the exact release with Zenodo or an equivalent service.
- Add the software DOI to `CITATION.cff`, the README, and the manuscript.
- Preserve future development on a new branch or release so the submitted snapshot remains immutable.

## Final security check

Review both tracked files and Git history. Deleting a secret from the latest commit is not enough if it appeared in earlier history. Rotate any credential that was ever exposed.
