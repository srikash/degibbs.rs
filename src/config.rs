//! Resolve CLI flags + BIDS sidecar into `Params` (explicit flags win over the sidecar).

use crate::{
    Params,
    pf::Pf,
    shift::ShiftParams,
    sidecar::{Sidecar, parse_axis},
};
use anyhow::{Result, anyhow, bail};

#[derive(Debug, Default, Clone)]
pub struct Requested {
    pub pe_axis: Option<String>,
    pub axes: Option<String>,
    pub pf: Option<String>,
    pub shift: ShiftParams,
    /// The input is complex (a complex NIfTI, or a magnitude image with `--phase`).
    pub complex: bool,
}

#[derive(Debug)]
pub struct Resolved {
    pub params: Params,
    pub notices: Vec<String>,
}

/// Parse `1`, `0.875`, `7/8`, ... into a supported PF mode.
pub fn parse_pf(s: &str) -> Result<Pf> {
    let v = match s.split_once('/') {
        Some((a, b)) => a.trim().parse::<f64>()? / b.trim().parse::<f64>()?,
        None => s.trim().parse::<f64>()?,
    };
    Pf::from_fraction(v)
}

fn parse_axes(s: &str) -> Result<(usize, usize)> {
    let (a, b) = s
        .split_once(',')
        .ok_or_else(|| anyhow!("--axes expects two comma-separated axes, e.g. 0,1"))?;
    let (a, b) = (parse_axis(a)?, parse_axis(b)?);
    if a == b {
        bail!("--axes must name two different axes");
    }
    Ok((a, b))
}

pub fn resolve(req: &Requested, sidecar: Option<&Sidecar>, ndim: usize) -> Result<Resolved> {
    let mut notices = Vec::new();
    let pe = match (&req.pe_axis, sidecar) {
        (Some(s), _) => parse_axis(s)?,
        (None, Some(sc)) => sc.pe_axis()?.ok_or_else(|| {
            anyhow!(
                "no phase-encoding axis: pass --pe-axis (the sidecar has no PhaseEncodingDirection)"
            )
        })?,
        (None, None) => {
            bail!("no phase-encoding axis: pass --pe-axis or provide a BIDS JSON sidecar")
        }
    };

    let plane = match (
        &req.axes,
        sidecar.map(|s| s.slice_axis()).transpose()?.flatten(),
    ) {
        (Some(s), _) => parse_axes(s)?,
        (None, Some(slice)) => {
            let rest: Vec<usize> = (0..3).filter(|&a| a != slice).collect();
            (rest[0], rest[1])
        }
        (None, None) => (0, 1),
    };
    let other = if plane.0 == pe {
        plane.1
    } else if plane.1 == pe {
        plane.0
    } else {
        bail!(
            "phase-encoding axis {pe} is not in the unringing plane (axes {},{}); pass --axes",
            plane.0,
            plane.1
        );
    };
    if pe >= ndim || other >= ndim {
        bail!("axes ({pe}, {other}) do not exist in a {ndim}D image");
    }

    let sidecar_pf = sidecar.and_then(|s| s.partial_fourier);
    let pf = if req.complex {
        // Complex data is taken to be POCS-reconstructed, so its k-space again ends symmetrically.
        if let Some(s) = &req.pf {
            let pf = parse_pf(s)?;
            if pf != Pf::Full {
                bail!(
                    "complex images are supported for full-Fourier data only: drop --pf or use --pf 1 \
                     (zero-filled partial-Fourier data should be processed as magnitude)"
                );
            }
        } else if let Some(v) = sidecar_pf.filter(|&v| v < 0.99) {
            notices.push(format!(
                "complex input is assumed to be POCS-reconstructed full-Fourier data: processing as --pf 1 \
                 and ignoring PartialFourier {v} in the sidecar (zero-filled data would keep its partial-Fourier ringing)"
            ));
        }
        Pf::Full
    } else {
        match (&req.pf, sidecar_pf) {
            (Some(s), _) => parse_pf(s)?,
            (None, Some(v)) => {
                let pf = Pf::from_fraction(v)?;
                if pf != Pf::Full {
                    notices.push(format!(
                        "partial Fourier factor {v} taken from the sidecar; if the data were POCS or homodyne reconstructed, pass --pf 1"
                    ));
                }
                pf
            }
            (None, None) => {
                notices.push("no partial Fourier factor given (no --pf, no PartialFourier in a sidecar): assuming full Fourier".into());
                Pf::Full
            }
        }
    };
    Ok(Resolved {
        params: Params {
            shift: req.shift,
            pf,
            pe_axis: pe,
            other_axis: other,
        },
        notices,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sc(json: &str) -> Sidecar {
        serde_json::from_str(json).unwrap()
    }

    #[test]
    fn sidecar_only() {
        let r = resolve(
            &Requested::default(),
            Some(&sc(
                r#"{"PhaseEncodingDirection":"j-","PartialFourier":0.75}"#,
            )),
            4,
        )
        .unwrap();
        assert_eq!(
            (r.params.pe_axis, r.params.other_axis, r.params.pf),
            (1, 0, Pf::P68)
        );
        assert_eq!(r.notices.len(), 1);
        assert!(r.notices[0].contains("pass --pf 1"));
    }

    #[test]
    fn flags_override_sidecar() {
        let req = Requested {
            pe_axis: Some("i".into()),
            pf: Some("7/8".into()),
            ..Default::default()
        };
        let r = resolve(
            &req,
            Some(&sc(
                r#"{"PhaseEncodingDirection":"j","PartialFourier":0.75}"#,
            )),
            4,
        )
        .unwrap();
        assert_eq!(
            (r.params.pe_axis, r.params.other_axis, r.params.pf),
            (0, 1, Pf::P78)
        );
    }

    fn complex_req(pf: Option<&str>) -> Requested {
        Requested {
            pe_axis: Some("j".into()),
            pf: pf.map(Into::into),
            complex: true,
            ..Default::default()
        }
    }

    #[test]
    fn complex_input_runs_as_full_fourier_and_ignores_the_sidecar() {
        let json = sc(r#"{"PartialFourier":0.75}"#);
        let r = resolve(&complex_req(None), Some(&json), 4).unwrap();
        assert_eq!(r.params.pf, Pf::Full);
        assert_eq!(r.notices.len(), 1);
        assert!(r.notices[0].contains("POCS"));
        // no notice when there is nothing to ignore
        let r = resolve(&complex_req(None), Some(&sc(r#"{"PartialFourier":1}"#)), 4).unwrap();
        assert!(r.notices.is_empty());
        let r = resolve(&complex_req(None), None, 4).unwrap();
        assert_eq!(r.params.pf, Pf::Full);
        assert!(r.notices.is_empty());
    }

    #[test]
    fn complex_input_accepts_explicit_full_fourier_only() {
        let json = sc(r#"{"PartialFourier":0.75}"#);
        let r = resolve(&complex_req(Some("1")), Some(&json), 4).unwrap();
        assert_eq!(r.params.pf, Pf::Full);
        assert!(r.notices.is_empty());
        for bad in ["6/8", "7/8", "5/8", "0.75"] {
            let e = resolve(&complex_req(Some(bad)), Some(&json), 4).unwrap_err();
            assert!(format!("{e}").contains("full-Fourier"), "{bad}: {e}");
        }
    }

    #[test]
    fn missing_pe_axis_is_an_error() {
        assert!(resolve(&Requested::default(), None, 4).is_err());
        assert!(resolve(&Requested::default(), Some(&sc("{}")), 4).is_err());
    }

    #[test]
    fn absent_pf_assumes_full_fourier_with_notice() {
        let req = Requested {
            pe_axis: Some("j".into()),
            ..Default::default()
        };
        let r = resolve(&req, None, 3).unwrap();
        assert_eq!(r.params.pf, Pf::Full);
        assert_eq!(r.notices.len(), 1);
    }

    #[test]
    fn slice_direction_picks_the_plane() {
        let r = resolve(
            &Requested::default(),
            Some(&sc(
                r#"{"PhaseEncodingDirection":"k","SliceEncodingDirection":"i"}"#,
            )),
            4,
        )
        .unwrap();
        assert_eq!((r.params.pe_axis, r.params.other_axis), (2, 1));
    }

    #[test]
    fn pe_axis_outside_default_plane_needs_axes() {
        let req = Requested {
            pe_axis: Some("k".into()),
            ..Default::default()
        };
        assert!(resolve(&req, None, 4).is_err());
        let req = Requested {
            pe_axis: Some("k".into()),
            axes: Some("1,2".into()),
            ..Default::default()
        };
        assert_eq!(resolve(&req, None, 4).unwrap().params.other_axis, 1);
    }

    #[test]
    fn pf_parsing_and_validation() {
        assert_eq!(parse_pf("1").unwrap(), Pf::Full);
        assert_eq!(parse_pf("0.875").unwrap(), Pf::P78);
        assert_eq!(parse_pf("5/8").unwrap(), Pf::P58);
        assert!(parse_pf("0.5").is_err());
        assert!(parse_pf("abc").is_err());
    }
}
