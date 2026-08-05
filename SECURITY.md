# Security policy

## Supported version

Until the first stable release, security fixes are applied to the latest code
on the `main` branch.

## Reporting a vulnerability

Please do not open a public issue for a suspected vulnerability. Use GitHub's
private vulnerability reporting feature for this repository. Include the
affected version or commit, reproduction steps, impact, and any suggested
mitigation.

For ordinary bugs that do not expose data or systems, open a regular GitHub
issue with a minimal reproducible example.

## Data handling

The R package analyzes local CSV data and does not transmit experimental data.
The optional FlowJo export scripts process local files. Do not add telemetry,
remote uploads, or network calls without an explicit design and security
review.
