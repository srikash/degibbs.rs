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
}

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

    let pf = match (&req.pf, sidecar.and_then(|s| s.partial_fourier)) {
        (Some(s), _) => parse_pf(s)?,
        (None, Some(v)) => Pf::from_fraction(v)?,
        (None, None) => {
            notices.push("no partial Fourier factor given (no --pf, no PartialFourier in a sidecar): assuming full Fourier".into());
            Pf::Full
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
        assert!(r.notices.is_empty());
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
