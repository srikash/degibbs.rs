# syntax=docker/dockerfile:1
FROM rust:1-alpine AS build
RUN apk add --no-cache musl-dev
# The image is built from the crate published on crates.io, so it always matches `cargo install degibbs`.
ARG VERSION
RUN test -n "$VERSION" && cargo install degibbs --version "$VERSION" --locked --root /out

FROM alpine:3
LABEL org.opencontainers.image.title="degibbs" \
      org.opencontainers.image.description="Gibbs-ringing removal for MRI with partial-Fourier support" \
      org.opencontainers.image.source="https://github.com/srikash/degibbs.rs" \
      org.opencontainers.image.licenses="LicenseRef-NYU-Non-Commercial-Research"
RUN apk add --no-cache pigz
COPY --from=build /out/bin/degibbs /usr/local/bin/degibbs
COPY LICENCE /licenses/degibbs.rs-LICENCE
WORKDIR /data
ENTRYPOINT ["degibbs"]
