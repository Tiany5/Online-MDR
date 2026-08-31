#!/usr/bin/env bash
set -euo pipefail

if [[ "$(id -u)" -ne 0 ]]; then
  echo "This setup must run as root." >&2
  exit 1
fi

env_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
source_dir="$env_dir/sources"

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  r-base r-base-dev build-essential gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
  r-cran-biocmanager r-cran-matrix r-cran-foreach r-cran-doparallel \
  r-cran-ggplot2 r-cran-dplyr r-cran-patchwork r-cran-rcolorbrewer \
  r-cran-ggpubr r-cran-tidyr r-cran-iterators r-cran-glmnet \
  r-cran-mgcv r-cran-rstudioapi r-cran-scales r-cran-data.table \
  r-cran-kedd r-cran-rcppprogress

(
  cd "$source_dir"
  md5sum -c "$env_dir/checksums.md5"
)

install_if_needed() {
  local package="$1"
  local version="$2"
  local archive="$3"

  if Rscript -e 'args <- commandArgs(TRUE); ok <- requireNamespace(args[[1]], quietly=TRUE) && as.character(packageVersion(args[[1]])) == args[[2]]; quit(status=if (ok) 0L else 1L)' "$package" "$version"; then
    echo "$package $version is already installed"
  else
    MAKEFLAGS=-j4 R CMD INSTALL "$source_dir/$archive"
  fi
}

install_if_needed REBayes 2.60 REBayes_2.60.tar.gz
install_if_needed isotree 0.5.22 isotree_0.5.22.tar.gz
install_if_needed onlineFDR 2.2.0 onlineFDR_2.2.0.tar.gz

Rscript "$env_dir/audit_r_environment.R"

apt-get clean
rm -rf /var/lib/apt/lists/*

echo "R environment setup and verification completed successfully."
