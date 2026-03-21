# Contributing

Thanks for contributing to this open source research project.

## Scope

This repository supports coloc-confirmed TWAS analysis workflows for ankylosing spondylitis and related reproducible research outputs.

## License

By contributing, you agree your contributions are provided under the project MIT License.

## How To Contribute

1. Fork the repository and create a feature branch.
2. Keep changes focused and atomic.
3. Add or update tests for behavioral changes.
4. Update documentation for user-facing or workflow changes.
5. Open a pull request using the project PR template.

## Development Setup

1. Create and activate a Python environment.
2. Install dependencies from requirements.txt.
3. Run tests locally before opening a PR.

Example commands:

```bash
python -m venv .venv
# Windows Git Bash
source .venv/Scripts/activate
pip install -r requirements.txt
pytest -q
```

## Coding Guidelines

- Prefer small, readable functions and explicit naming.
- Keep public interfaces stable unless change is justified and documented.
- Avoid introducing silent fallbacks that hide real-data failures.
- Preserve reproducibility: document data sources, versions, and assumptions.

## Data And Security Guidelines

- Do not commit restricted or sensitive data.
- Use documented paths under data/ for local assets.
- Prefer checksums and provenance manifests when adding new datasets.

## Testing Expectations

- Add tests for new parsing, harmonization, validation, and pipeline logic.
- Run full test suite locally for substantial changes.
- Include commands and key outputs in PR description.

## Pull Request Expectations

Please include:

- Summary of changes
- Why the change is needed
- Files and modules affected
- Test evidence
- Risks and mitigation notes
- Follow-up tasks, if any

## Reporting Issues

Open an issue with:

- Problem description
- Steps to reproduce
- Expected vs actual behavior
- Relevant logs or stack traces
- Environment details (OS, Python version)

## Research Integrity Notes

- Clearly separate exploratory analyses from productionized pipeline steps.
- Prefer deterministic settings where possible.
- Record assumptions and thresholds in config and documentation.
