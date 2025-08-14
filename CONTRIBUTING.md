# Contributing to LikelihoodProfiler.jl

Thanks for your interest in improving **LikelihoodProfiler.jl**! This project welcomes small and simple contributions. If something here feels heavy — don’t worry, keep it simple.

## Questions & Support
- Please use **GitHub Issues** for questions and bug reports.
- Preferred language: **English**.

## Quick Start (for Pull Requests)
1. **Fork** the repo and create a branch from `master` (e.g., `feature/short-description`).
2. Install **Julia ≥ 1.10**.
3. In Julia REPL:
   ```julia
   import Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   Pkg.test()
   ```
4. Open a **Pull Request to `master`**. Describe what changed and why; link any related issues.

### What CI checks
- PRs run tests on **Ubuntu / macOS / Windows** for **Julia LTS** and **latest 1.x**.
- We merge when tests pass and a maintainer approves.

## Commit messages
- No strict format is required. Please use **clear, descriptive** messages (e.g., “Fix CI failure on Windows” rather than “fix”).

## Tests
- Add tests for new features and bug fixes.
- Keep tests reasonably fast and place them under `test/`.

## Dependencies
- Keep dependencies minimal. For new, non-trivial dependencies, please **open an issue** first.

## Versioning & Releases
- We follow **Semantic Versioning (SemVer)**.
- Releases and tagging are handled by maintainers.

## License
- By contributing, you agree your contributions are licensed under the **MIT License** of this repository.

## Code of conduct
- We do not maintain a separate code of conduct. Please be **respectful and constructive** in all interactions.

## Security
- If you believe you’ve found a security issue, please **open a private security advisory on GitHub** or contact the maintainers privately if possible. Please do not disclose publicly first.

## Good first issues
- Look for issues labeled **good first issue** or **help wanted** to get started.
