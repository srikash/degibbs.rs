# syntax=docker/dockerfile:1
FROM rust:1-alpine AS build
RUN apk add --no-cache musl-dev
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src src
RUN cargo build --release --locked

FROM alpine:3
LABEL org.opencontainers.image.title="degibbs" \
      org.opencontainers.image.description="Gibbs-ringing removal for magnitude MRI with partial-Fourier support" \
      org.opencontainers.image.source="https://github.com/srikash/degibbs.rs" \
      org.opencontainers.image.licenses="LicenseRef-NYU-Non-Commercial-Research"
COPY --from=build /src/target/release/degibbs /usr/local/bin/degibbs
COPY LICENCE /licenses/degibbs.rs-LICENCE
WORKDIR /data
ENTRYPOINT ["degibbs"]
