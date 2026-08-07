#!/usr/bin/env bash

# Sets up environment and installs code linters.

set -euo pipefail
IFS=$'\n\t'

SCRIPT_NAME=$(basename "$0")

# Functions for logging and error handling.
log() { printf '%s: %s\n' "$SCRIPT_NAME" "$*" >&2; }
error() {
	log "ERROR: $*"
	exit 1
}
warn() { log "WARN: $*"; }

# Requires that a command is available in PATH, else exits with error.
require_cmd() {
	command -v "$1" >/dev/null 2>&1 || error "required command '$1' not found in PATH"
}

# Installs packages using apt-get, exits on failure.
install_packages() {
	local pkgs=("$@")
	if [ "${#pkgs[@]}" -eq 0 ]; then
		warn "No packages specified to install"
		return
	fi
	log "Installing packages: ${pkgs[*]}"
	apt-get install -yqq --no-install-recommends "${pkgs[@]}"
}

# Installs shellcheck if not already installed.
install_shellcheck() {
	if command -v shellcheck >/dev/null 2>&1; then
		log "shellcheck already installed."
		return
	fi
	install_packages shellcheck
	log "shellcheck installed."
}

# Adds the official LLVM apt repository and pins clang-format-<version> to it
# (apt otherwise picks whichever archive reports the higher version). The pin
# covers the clang-format package only, not its LLVM runtime dependencies.
# Arguments:
#   $1 : clang-format major version (e.g., 20)
add_llvm_apt_repo() {
	local clang_format_version="$1"
	local codename
	# shellcheck disable=SC1091 # /etc/os-release is a standard runtime file, not a repo input.
	codename=$(. /etc/os-release && echo "$VERSION_CODENAME")

	install_packages wget gnupg ca-certificates

	local keyring="/usr/share/keyrings/apt-llvm-org.gpg"
	wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | gpg --batch --yes --dearmor -o "$keyring"

	printf 'deb [signed-by=%s] http://apt.llvm.org/%s/ llvm-toolchain-%s-%s main\n' \
		"$keyring" "$codename" "$codename" "$clang_format_version" \
		>/etc/apt/sources.list.d/llvm.list

	printf 'Package: clang-format-%s\nPin: origin apt.llvm.org\nPin-Priority: 1001\n' \
		"$clang_format_version" >/etc/apt/preferences.d/clang-format-llvm.pref

	apt-get update -qq
}

# Installs specified version of clang-format from apt repositories.
# Arguments:
#   -c | --clang-format-version <version> : Version number (e.g., 20 for clang-format-20)
# Example:
#   install_clang_format_version -c 20
install_clang_format_version() {
	while [ "$#" -gt 0 ]; do
		case "$1" in
		-c | --clang-format-version)
			shift
			[ "$#" -gt 0 ] || error "missing value for $1"
			CLANG_FORMAT_VER="$1"
			;;
		--)
			shift
			break
			;;
		*)
			error "unknown argument to install_clang_format_from_llvm: $1"
			;;
		esac
		shift
	done

	if ! printf '%s' "$CLANG_FORMAT_VER" | grep -Eq '^[0-9]+$'; then
		error "invalid clang-format version: '$CLANG_FORMAT_VER'. Must be an integer like '20'"
	fi

	local pkg_name="clang-format-$CLANG_FORMAT_VER"

	installed_packages=$(apt list --installed 2>/dev/null | awk -F/ '{print $1}')
	if printf '%s' "$installed_packages" | grep -q "$pkg_name"; then
		log "$pkg_name already installed."
		return
	fi

	install_packages "$pkg_name"
	log "$pkg_name installed: $(dpkg-query -W -f='${Version}' "$pkg_name")"
}

apt-get update -qq
install_shellcheck

clang_format_version="20"
add_llvm_apt_repo "$clang_format_version"
install_clang_format_version --clang-format-version "$clang_format_version"

log "All requested linters installed"
