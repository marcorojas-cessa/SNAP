## Summary
- What changed?
- Why was it needed?

## Scope
- [ ] Core change (`SNAP*.m`, `+snap_helpers`, `+snap_modules`)
- [ ] Contribution-only change (`+snap_contrib`, `examples`, `external_plugins`)
- [ ] Documentation-only

Label note:
- PRs touching core paths must include the `core` label.
- PRs labeled `contribution` must not touch core paths.

## Scientific/Behavioral Impact
- Expected effect on localization, fitting, segmentation, or classification behavior:
- Any channel-specific behavior introduced:

## Validation
- [ ] I ran this change locally
- [ ] I verified no UI/runtime errors in affected workflow(s)
- [ ] I verified channel-homogeneous behavior where applicable
- [ ] I updated docs for user-facing behavior

### Local checks performed
- MATLAB version:
- Commands run:
- Key results:

## Data/Schema Impact
- [ ] No export schema change
- [ ] Export schema changed (describe fields and why):

## Compatibility (Current SNAP State)
- [ ] Compatible with current SNAP modules/apps (`SNAP`, `SNAP_batch`, `SNAP_prepare`, `SNAP_classify`, `SNAP_train`)

## Notes for Reviewers
- Areas needing close review:
- Any follow-up work:
