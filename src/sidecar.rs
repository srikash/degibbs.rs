//! BIDS JSON sidecar parsing (keys per the BIDS specification).

use anyhow::{Context, Result, bail};
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// Axis written as `i|j|k` (or `0|1|2`) with an optional `+`/`-` sign.
pub fn parse_axis(s: &str) -> Result<usize> {
    let t = s.trim().trim_end_matches(['+', '-']);
    match t {
        "i" | "0" => Ok(0),
        "j" | "1" => Ok(1),
        "k" | "2" => Ok(2),
        _ => bail!("invalid axis '{s}': expected i, j, k (optionally signed) or 0, 1, 2"),
    }
}

#[derive(Debug, Default, Deserialize, PartialEq)]
#[serde(rename_all = "PascalCase")]
pub struct Sidecar {
    pub phase_encoding_direction: Option<String>,
    /// Written by some converters instead of `PhaseEncodingDirection` when the polarity is unknown.
    pub phase_encoding_axis: Option<String>,
    pub slice_encoding_direction: Option<String>,
    pub partial_fourier: Option<f64>,
}

impl Sidecar {
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading sidecar {}", path.display()))?;
        serde_json::from_str(&text).with_context(|| format!("parsing sidecar {}", path.display()))
    }

    pub fn pe_axis(&self) -> Result<Option<usize>> {
        self.phase_encoding_direction
            .as_deref()
            .or(self.phase_encoding_axis.as_deref())
            .map(parse_axis)
            .transpose()
    }

    pub fn slice_axis(&self) -> Result<Option<usize>> {
        self.slice_encoding_direction
            .as_deref()
            .map(parse_axis)
            .transpose()
    }
}

/// `dir/name.nii.gz` or `dir/name.nii` -> `dir/name.json`.
pub fn default_path(input: &Path) -> PathBuf {
    let name = input
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let stem = name
        .strip_suffix(".nii.gz")
        .or_else(|| name.strip_suffix(".nii"))
        .unwrap_or(&name);
    input.with_file_name(format!("{stem}.json"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axes_parse_with_and_without_sign() {
        for (s, want) in [("i", 0), ("j-", 1), ("k", 2), ("0", 0), ("1+", 1)] {
            assert_eq!(parse_axis(s).unwrap(), want, "{s}");
        }
        assert!(parse_axis("x").is_err());
        assert!(parse_axis("").is_err());
    }

    #[test]
    fn default_sidecar_path_strips_double_extension() {
        assert_eq!(
            default_path(Path::new("/a/b/in.nii.gz")),
            PathBuf::from("/a/b/in.json")
        );
        assert_eq!(default_path(Path::new("in.nii")), PathBuf::from("in.json"));
        assert_eq!(
            default_path(Path::new("weird")),
            PathBuf::from("weird.json")
        );
    }

    #[test]
    fn parses_bids_fields_and_ignores_the_rest() {
        let s: Sidecar = serde_json::from_str(
            r#"{"PhaseEncodingDirection":"j-","PartialFourier":0.75,"SliceEncodingDirection":"k","EchoTime":0.1}"#,
        )
        .unwrap();
        assert_eq!(s.pe_axis().unwrap(), Some(1));
        assert_eq!(s.slice_axis().unwrap(), Some(2));
        assert_eq!(s.partial_fourier, Some(0.75));
        let empty: Sidecar = serde_json::from_str("{}").unwrap();
        assert_eq!(empty.pe_axis().unwrap(), None);
    }

    #[test]
    fn phase_encoding_axis_is_a_fallback_and_direction_wins() {
        let s: Sidecar = serde_json::from_str(r#"{"PhaseEncodingAxis":"i"}"#).unwrap();
        assert_eq!(s.pe_axis().unwrap(), Some(0));
        let both: Sidecar =
            serde_json::from_str(r#"{"PhaseEncodingAxis":"i","PhaseEncodingDirection":"j-"}"#)
                .unwrap();
        assert_eq!(both.pe_axis().unwrap(), Some(1));
    }

    #[test]
    fn invalid_direction_is_an_error() {
        let s: Sidecar = serde_json::from_str(r#"{"PhaseEncodingDirection":"AP"}"#).unwrap();
        assert!(s.pe_axis().is_err());
    }
}
